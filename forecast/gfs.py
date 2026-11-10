"""Download a small subset of the NOAA GFS 0.25° forecast from NOMADS (grib filter) and read it."""
import os, time
from datetime import datetime, timedelta, timezone
import numpy as np
import requests

NOMADS = 'https://nomads.ncep.noaa.gov'
FILTER = NOMADS + '/cgi-bin/filter_gfs_0p25.pl'
PRES_LEVELS = [1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50, 40, 30, 20, 15, 10, 7, 5, 3]
FLOAT_LEVELS = [1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50, 40, 30, 20, 15, 10]   # floating: wide area, up to ~31 km
SESSION = requests.Session()
SESSION.headers['User-Agent'] = 'stratolaunch-forecast (github actions)'


def _get(url, params=None, tries=4):
    for k in range(tries):
        try:
            r = SESSION.get(url, params=params, timeout=90)
            if r.status_code == 200:
                return r
            if r.status_code == 404:
                return None
        except requests.RequestException:
            pass
        time.sleep(5 * (k + 1))
    return None


def latest_cycle(now=None, need_hour=384):
    """Newest GFS cycle whose last forecast hour is already published."""
    now = now or datetime.now(timezone.utc)
    t = now.replace(minute=0, second=0, microsecond=0, hour=now.hour // 6 * 6)
    for _ in range(8):
        d, h = t.strftime('%Y%m%d'), t.strftime('%H')
        url = f'{NOMADS}/pub/data/nccf/com/gfs/prod/gfs.{d}/{h}/atmos/gfs.t{h}z.pgrb2.0p25.f{need_hour:03d}.idx'
        try:
            if SESSION.head(url, timeout=30).status_code == 200:
                return t
        except requests.RequestException:
            pass
        t -= timedelta(hours=6)
    raise RuntimeError('no complete GFS cycle found on NOMADS in the last 48 h')


def _params(cycle, fh, box, kind):
    d, h = cycle.strftime('%Y%m%d'), cycle.strftime('%H')
    p = {'dir': f'/gfs.{d}/{h}/atmos', 'file': f'gfs.t{h}z.pgrb2.0p25.f{fh:03d}',
         'subregion': '', 'toplat': box['north'], 'bottomlat': box['south'],
         'leftlon': box['west'], 'rightlon': box['east']}
    if kind in ('upper', 'float'):
        p.update({'var_UGRD': 'on', 'var_VGRD': 'on', 'var_HGT': 'on'})
        p.update({f'lev_{L}_mb': 'on' for L in (PRES_LEVELS if kind == 'upper' else FLOAT_LEVELS)})
    else:
        p.update({v: 'on' for v in ['var_VIS', 'var_GUST', 'var_TMP', 'var_DPT', 'var_UGRD', 'var_VGRD',
                                    'var_PRATE', 'var_CAPE', 'var_LCDC', 'var_TCDC']})
        p.update({v: 'on' for v in ['lev_surface', 'lev_2_m_above_ground', 'lev_10_m_above_ground',
                                    'lev_low_cloud_layer', 'lev_entire_atmosphere']})
    return p


def download(cycle, hours, box, cache_dir, pause=0.7):
    """Fetch upper-air and surface subsets for each forecast hour. Returns {fh: (upper_path, sfc_path)}."""
    os.makedirs(cache_dir, exist_ok=True)
    out = {}
    for fh in sorted(set(hours)):
        paths = []
        for kind in ('upper', 'sfc'):
            path = os.path.join(cache_dir, f'{cycle:%Y%m%d%H}_f{fh:03d}_{kind}.grb2')
            if not (os.path.exists(path) and os.path.getsize(path) > 0):
                r = _get(FILTER, _params(cycle, fh, box, kind))
                if r is None or not r.content.startswith(b'GRIB'):
                    path = None
                else:
                    with open(path, 'wb') as f:
                        f.write(r.content)
                time.sleep(pause)
            paths.append(path)
        if paths[0]:
            out[fh] = tuple(paths)
        else:
            print(f'  f{fh:03d}: upper-air download failed, skipped')
    return out


def download_float(cycle, hours, box, cache_dir, pause=0.7):
    """Upper air over a wide area for the floating drift. Returns {fh: path}."""
    os.makedirs(cache_dir, exist_ok=True)
    out = {}
    for fh in sorted(set(hours)):
        path = os.path.join(cache_dir, f'{cycle:%Y%m%d%H}_f{fh:03d}_float.grb2')
        if not (os.path.exists(path) and os.path.getsize(path) > 0):
            r = _get(FILTER, _params(cycle, fh, box, 'float'))
            if r is None or not r.content.startswith(b'GRIB'):
                print(f'  f{fh:03d}: floating download failed, skipped')
                continue
            with open(path, 'wb') as f:
                f.write(r.content)
            time.sleep(pause)
        out[fh] = path
    return out


# ----------------------------------------------------------------- reading
SFC_MAP = {  # (shortName, typeOfLevel or None) -> our name
    ('vis', None): 'vis_m', ('gust', None): 'gust_ms', ('2t', None): 't2_k', ('2d', None): 'td2_k',
    ('10u', None): 'u10', ('10v', None): 'v10', ('prate', None): 'prate', ('cape', None): 'cape',
    ('lcc', None): 'low_cloud_pct', ('tcc', None): 'total_cloud_pct'}


def read_upper(path, stride=1):
    """-> lats (ascending), lons, levels_hPa (descending pressure), u, v, gh  each (L, ny, nx). stride: keep every n-th point."""
    import pygrib
    fields = {}
    lats = lons = None
    with pygrib.open(path) as g:
        for m in g:
            if m.typeOfLevel != 'isobaricInhPa' or m.shortName not in ('u', 'v', 'gh'):
                continue
            vals = np.asarray(m.values, dtype='f4')
            if lats is None:
                la, lo = m.latlons(); lats, lons = la[:, 0], lo[0, :]
            fields[(m.shortName, int(m.level))] = vals
    levels = sorted({lv for (_, lv) in fields}, reverse=True)
    levels = [lv for lv in levels if all((n, lv) in fields for n in ('u', 'v', 'gh'))]
    stack = {n: np.stack([fields[(n, lv)] for lv in levels]) for n in ('u', 'v', 'gh')}
    if lats[0] > lats[-1]:
        lats = lats[::-1]
        stack = {n: a[:, ::-1, :] for n, a in stack.items()}
    if stride > 1:
        lats, lons = lats[::stride], lons[::stride]
        stack = {n: a[:, ::stride, ::stride] for n, a in stack.items()}
    return lats, lons, np.array(levels), stack['u'], stack['v'], stack['gh']


def read_surface(path):
    """-> lats, lons, {name: (ny, nx)}"""
    import pygrib
    out, lats, lons = {}, None, None
    if not path:
        return lats, lons, out
    with pygrib.open(path) as g:
        for m in g:
            key = None
            for (sn, tl), name in SFC_MAP.items():
                if m.shortName == sn and (tl is None or m.typeOfLevel == tl):
                    key = name
                    break
            if key is None or key in out:
                continue
            if lats is None:
                la, lo = m.latlons(); lats, lons = la[:, 0], lo[0, :]
            out[key] = np.asarray(m.values, dtype='f4')
    if lats is not None and lats[0] > lats[-1]:
        lats = lats[::-1]
        out = {k: v[::-1, :] for k, v in out.items()}
    return lats, lons, out
