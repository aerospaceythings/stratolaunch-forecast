"""Offline test: write fake GFS GRIB subsets and run the whole pipeline (python -m tests.synthetic)."""
import os, sys
from datetime import datetime, timezone
import numpy as np
import eccodes as ec

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from forecast import gfs, run, site
import yaml

TMP = os.path.join(ROOT, '.cache', 'synthetic')
lats = np.arange(-34.75, -28.49, 0.25)
lons = np.arange(113.5, 120.76, 0.25)


def _grid(h, sample, lats=lats, lons=lons):
    g = ec.codes_grib_new_from_samples(sample)
    for k, v in [('Ni', len(lons)), ('Nj', len(lats)), ('latitudeOfFirstGridPointInDegrees', lats[-1]),
                 ('longitudeOfFirstGridPointInDegrees', lons[0]), ('latitudeOfLastGridPointInDegrees', lats[0]),
                 ('longitudeOfLastGridPointInDegrees', lons[-1]), ('iDirectionIncrementInDegrees', 0.25),
                 ('jDirectionIncrementInDegrees', 0.25), ('jScansPositively', 0)]:
        ec.codes_set(g, k, v)
    return g


def write(fh, path_up, path_sf, day_seed):
    rng = np.random.default_rng(fh)
    with open(path_up, 'wb') as f:
        for lev in gfs.PRES_LEVELS:
            z = 7200 * np.log(1013 / lev)
            jet = 35 * np.exp(-((z - 12000) / 4000) ** 2) + (8 if day_seed % 3 else 2)
            u = jet - (12 if z > 21000 else 0) + rng.normal(0, 1)
            v = 3 * np.sin(z / 4000 + day_seed)
            for sn, val in [('u', u), ('v', v), ('gh', z)]:
                g = _grid(fh, 'regular_ll_pl_grib2')
                ec.codes_set(g, 'shortName', sn)
                ec.codes_set(g, 'level', lev)
                ec.codes_set_values(g, np.full(len(lats) * len(lons), val, float))
                ec.codes_write(g, f); ec.codes_release(g)
    with open(path_sf, 'wb') as f:
        for sn, val, tl, lv in [('2t', 288.0, 'heightAboveGround', 2), ('2d', 286.5 if day_seed % 4 == 0 else 280, 'heightAboveGround', 2),
                                ('10u', 2.0 + day_seed % 5, 'heightAboveGround', 10), ('10v', 1.0, 'heightAboveGround', 10),
                                ('vis', 20000.0 if day_seed % 4 else 800.0, 'surface', 0), ('gust', 6.0, 'surface', 0),
                                ('cape', 50.0, 'surface', 0), ('prate', 0.0, 'surface', 0), ('orog', 280.0, 'surface', 0)]:
            g = _grid(fh, 'regular_ll_sfc_grib2')
            try:
                ec.codes_set(g, 'typeOfLevel', tl); ec.codes_set(g, 'level', lv)
                ec.codes_set(g, 'shortName', sn)
            except Exception as e:
                print('skip', sn, e); ec.codes_release(g); continue
            if sn == 'orog':
                for k, v in [('discipline', 0), ('parameterCategory', 3), ('parameterNumber', 5), ('typeOfFirstFixedSurface', 1)]:
                    ec.codes_set(g, k, v)
            ec.codes_set_values(g, np.full(len(lats) * len(lons), val, float))
            ec.codes_write(g, f); ec.codes_release(g)


def fake_download(cycle, hours, box, cache_dir, pause=0):
    os.makedirs(TMP, exist_ok=True)
    out = {}
    for fh in hours:
        up, sf = os.path.join(TMP, f'f{fh:03d}_up.grb2'), os.path.join(TMP, f'f{fh:03d}_sf.grb2')
        if not os.path.exists(up):
            write(fh, up, sf, fh // 24)
        out[fh] = (up, sf)
    return out


FL_LATS, FL_LONS = np.arange(-41, -21.99, 0.25), np.arange(105, 128.01, 0.25)


def fake_download_float(cycle, hours, box, cache_dir, pause=0):
    """Wide-area upper air: stratospheric easterlies that slowly change speed and meander north-south."""
    os.makedirs(TMP, exist_ok=True)
    out = {}
    for fh in hours:
        p = os.path.join(TMP, f'f{fh:03d}_float.grb2')
        if not os.path.exists(p):
            with open(p, 'wb') as f:
                for lev in gfs.FLOAT_LEVELS:
                    z = 7200 * np.log(1013 / lev)
                    if z < 17000:
                        u, v = 35 * np.exp(-((z - 12000) / 4000) ** 2) + 6, 2.0
                    else:
                        u = -(3 + 5 * np.sin(fh / 40)) - (z - 20000) / 1500
                        v = 4 * np.sin(fh / 17 + z / 3000)
                    for sn, val in [('u', u), ('v', v), ('gh', z)]:
                        g = _grid(fh, 'regular_ll_pl_grib2', FL_LATS, FL_LONS)
                        ec.codes_set(g, 'shortName', sn); ec.codes_set(g, 'level', lev)
                        ec.codes_set_values(g, np.full(len(FL_LATS) * len(FL_LONS), val, float))
                        ec.codes_write(g, f); ec.codes_release(g)
        out[fh] = p
    return out


if __name__ == '__main__':
    gfs.download = fake_download
    gfs.download_float = fake_download_float
    cfg = yaml.safe_load(open(os.path.join(ROOT, 'config.yaml')))
    cyc = datetime(2026, 11, 9, 0, tzinfo=timezone.utc)
    out = run.run(cfg, cyc, int(sys.argv[1]) if len(sys.argv) > 1 else 16)
    run.save(out)
    site.build(ROOT, cfg)
