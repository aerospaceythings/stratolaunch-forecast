"""Historical outlook for days 17-30: the same launch, flown through ERA5 winds of the previous 5 years.

python -m forecast.history            # fill/refresh data/history_window/*.json
Needs a Copernicus CDS key: env CDSAPI_KEY (GitHub secret) or ~/.cdsapirc.
Only dates that are new to the window are downloaded, so a normal day fetches one date x 5 years.
"""
import glob, json, math, os
from datetime import datetime, timedelta, timezone, date
import numpy as np
import yaml

from . import model, sun

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'data', 'history_window')
CACHE = os.path.join(ROOT, '.cache', 'era5')
LEVELS = ['1000', '925', '850', '700', '600', '500', '400', '300', '250', '200', '150',
          '125', '100', '70', '50', '30', '20', '10', '7', '5']
TIMES = ['00:00', '03:00', '18:00', '21:00']
SFC_VARS = ['10m_u_component_of_wind', '10m_v_component_of_wind', '2m_temperature', '2m_dewpoint_temperature',
            'total_precipitation', 'convective_available_potential_energy', 'low_cloud_cover',
            'instantaneous_10m_wind_gust']
FIRST_DAY, LAST_DAY, NYEARS = 17, 30, 5


def client():
    import cdsapi
    key = os.environ.get('CDSAPI_KEY')
    if key:
        return cdsapi.Client(url=os.environ.get('CDSAPI_URL', 'https://cds.climate.copernicus.eu/api'),
                             key=key, quiet=True, progress=False)
    return cdsapi.Client(quiet=True, progress=False)


def _open_any(path):
    """CDS returns a zip (one .nc per stream/step type) for mixed single-level requests."""
    import xarray as xr, zipfile
    if zipfile.is_zipfile(path):
        d = path + '_unzipped'
        with zipfile.ZipFile(path) as z:
            z.extractall(d)
        parts = [xr.open_dataset(os.path.join(d, n), engine='netcdf4') for n in sorted(os.listdir(d)) if n.endswith('.nc')]
        return xr.merge(parts, compat='override', join='outer')
    return xr.open_dataset(path, engine='netcdf4')


def _day_path(kind, d):
    return os.path.join(CACHE, f'{kind}_{d:%Y%m%d}.nc')


def _split(ds, kind, keep_days):
    """Store a downloaded block as one file per UTC day, so every day is downloaded only once."""
    import pandas as pd
    for d in sorted({pd.Timestamp(t).date() for t in ds['time'].values}):
        part = ds.sel(time=slice(np.datetime64(d), np.datetime64(d + timedelta(days=1)) - np.timedelta64(1, 's')))
        if keep_days is None or d in keep_days:
            part.load().to_netcdf(_day_path(kind, d) + '.part', engine='netcdf4')
            os.replace(_day_path(kind, d) + '.part', _day_path(kind, d))


def migrate_cache():
    """Older runs cached whole requests (pl_YYYYMM_dd-dd.nc). Split them into per-day files, no re-download."""
    for f in glob.glob(os.path.join(CACHE, '*_??????_*.nc')):
        kind = os.path.basename(f).split('_')[0]
        try:
            ds = open_nc(f)
            _split(ds[['u', 'v', 'z']] if kind == 'pl' else ds, kind, None)
            print(f'  cache: split {os.path.basename(f)} into days')
        except Exception as e:
            print(f'  cache: could not split {os.path.basename(f)} ({e})')
        os.remove(f)
    for d in glob.glob(os.path.join(CACHE, '*_unzipped')):
        import shutil; shutil.rmtree(d, ignore_errors=True)


def load_days(days, box):
    """ERA5 for the given UTC days: from the per-day cache, downloading only the days that are missing
    (one request per month and dataset). Days no longer needed are removed from the cache."""
    import xarray as xr
    os.makedirs(CACHE, exist_ok=True)
    migrate_cache()
    want = set(days)
    for f in glob.glob(os.path.join(CACHE, '*_????????.nc')):            # keep the cache small
        try:
            if datetime.strptime(os.path.basename(f)[-11:-3], '%Y%m%d').date() not in want:
                os.remove(f)
        except ValueError:
            pass
    area = [box['north'], box['west'], box['south'], box['east']]
    missing = {}
    for d in days:
        if not (os.path.exists(_day_path('pl', d)) and os.path.exists(_day_path('sfc', d))):
            missing.setdefault((d.year, d.month), []).append(d.day)
    print(f'  ERA5: {len(days) - sum(map(len, missing.values()))} days cached, '
          f'{sum(map(len, missing.values()))} to download', flush=True)
    if missing:
        c = client()
    for (y, m), dd in sorted(missing.items()):
        base = {'product_type': ['reanalysis'], 'year': [str(y)], 'month': [f'{m:02d}'],
                'day': [f'{x:02d}' for x in dd], 'time': TIMES, 'area': area,
                'data_format': 'netcdf', 'download_format': 'unarchived'}
        print(f'  ERA5 download {y}-{m:02d} days {dd}', flush=True)
        tmp = os.path.join(CACHE, f'tmp_{y}{m:02d}')
        c.retrieve('reanalysis-era5-pressure-levels', dict(base, variable=['u_component_of_wind', 'v_component_of_wind',
                   'geopotential'], pressure_level=LEVELS), tmp + '_pl.nc')
        c.retrieve('reanalysis-era5-single-levels', dict(base, variable=SFC_VARS,
                   time=[f'{h:02d}:00' for h in range(24)]), tmp + '_sfc.nc')
        keep = {date(y, m, x) for x in dd}
        _split(open_nc(tmp + '_pl.nc')[['u', 'v', 'z']], 'pl', keep)
        _split(open_nc(tmp + '_sfc.nc'), 'sfc', keep)
        for f in glob.glob(tmp + '*'):
            if os.path.isdir(f):
                import shutil; shutil.rmtree(f, ignore_errors=True)
            else:
                os.remove(f)
    PL = xr.concat([xr.open_dataset(_day_path('pl', d), engine='netcdf4') for d in days], 'time').sortby('time')
    SF = xr.concat([xr.open_dataset(_day_path('sfc', d), engine='netcdf4') for d in days], 'time').sortby('time')
    return PL, SF


def open_nc(path):
    ds = _open_any(path)
    ren = {a: b for a, b in [('valid_time', 'time'), ('pressure_level', 'level'), ('latitude', 'lat'),
                             ('longitude', 'lon')] if a in ds.dims or a in ds.coords}
    ds = ds.rename(ren)
    ds = ds.drop_vars([v for v in ds.coords if v not in ('time', 'level', 'lat', 'lon')])
    return ds.sortby('lat').sortby('lon')


def flight_key(cfg, s):
    """Everything that changes a historical result: flight settings, site, limits and no-landing zones."""
    import hashlib
    f = cfg['flight']
    extra = json.dumps([s['lat'], s['lon'], s['ground_m'], cfg['burst_mode'], cfg.get('avoid_zones')], sort_keys=True)
    return f"{f['ascent_ms']}/{f['descent_sea_level_ms']}/{f['burst_km']}/" + hashlib.md5(extra.encode()).hexdigest()[:8]


def wanted(cfg, today):
    """Local calendar dates in the window (today + 17 .. today + 30)."""
    return [today + timedelta(days=k) for k in range(FIRST_DAY, LAST_DAY + 1)]


def fly_year(cfg, s, d_local, pl, sf):
    """One historical launch on the same calendar date in another year."""
    fl = cfg['flight']
    launch = sun.sunrise_utc(d_local, s['lat'], s['lon']) - timedelta(minutes=s['launch_minutes_before_sunrise'])
    t0 = np.datetime64(launch.replace(tzinfo=None))
    sel = pl.sel(time=slice(t0 - np.timedelta64(3, 'h'), t0 + np.timedelta64(6, 'h')))
    if sel.sizes['time'] < 2:
        return None
    ts = (sel['time'].values - t0) / np.timedelta64(1, 'h')          # hours relative to launch
    sel = sel.sortby('level', ascending=False)
    W = model.Winds(ts, sel['lat'].values, sel['lon'].values,
                    sel['u'].values, sel['v'].values, sel['z'].values / model.G0)
    (llat, llon), (blat, blon), dur, _ = model.fly(W, 0.0, s['lat'], s['lon'], s['ground_m'], fl)
    L, brg = model.dist_km(s['lat'], s['lon'], llat, llon)
    H, u, v = W.column(0.0, s['lat'], s['lon'])
    pm = model.profile_metrics(H, u, v, cfg['float_mode']['band_km'], cfg['float_mode']['min_layer_speed_ms'])
    # surface at the site: launch value, and maxima from launch to landing
    pt = sf.interp(lat=s['lat'], lon=s['lon'])
    near = pt.interp(time=t0)
    win = pt.sel(time=slice(t0 - np.timedelta64(1, 'h'), t0 + np.timedelta64(int(dur) + 60, 'm')))
    g = lambda v: float(near[v]) if v in near and np.isfinite(float(near[v])) else None
    gmax = lambda v: float(win[v].max()) if v in win and win.sizes.get('time', 0) else None
    w10 = None if g('u10') is None else math.hypot(g('u10'), g('v10'))
    sunset = sun.sunset_utc(d_local, s['lat'], s['lon'])
    r = dict(landing_km=round(L, 1), landing_bearing=round(brg), landing_dir=model.compass(brg),
             landing_lat=round(llat, 4), landing_lon=round(llon, 4), landing_over_ocean=model.over_ocean(llat, llon),
             wind10_kmh=None if w10 is None else round(w10 * 3.6, 1),
             gust_kmh=None if g('i10fg') is None else round(g('i10fg') * 3.6, 1),
             rain_mmh=None if gmax('tp') is None else round(gmax('tp') * 1000, 2),
             cape=None if gmax('cape') is None else round(gmax('cape')),
             daylight_after_landing_h=round((sunset - launch).total_seconds() / 3600 - dur / 60, 1),
             jet_kmh=round(pm['jet_kmh']), band_min_ms=round(pm['band_min_ms'], 2),
             steering_arc_deg=round(pm['steering_arc_deg']))
    r['fog'] = model.fog_risk(None, g('t2m'), g('d2m'), w10, None if g('lcc') is None else g('lcc') * 100)
    r['burst_status'], checks = model.burst_status(r, cfg['burst_mode'], cfg.get('avoid_zones'))
    r['reasons'] = [c['text'] for c in checks if c['level'] != 'go']
    r['float_status'], _ = model.float_status(pm, cfg['float_mode'])
    return r


def summarise(per_year):
    ys = [v for v in per_year.values() if v]
    if not ys:
        return {}
    L = np.array([v['landing_km'] for v in ys])
    return dict(n=len(ys), go=sum(v['burst_status'] == 'go' for v in ys),
                marginal=sum(v['burst_status'] == 'marginal' for v in ys),
                float_good=sum(v['float_status'] == 'good' for v in ys),
                landing_median_km=round(float(np.median(L))), landing_min_km=round(float(L.min())),
                landing_max_km=round(float(L.max())))


def run(cfg, today=None):
    SITES = model.sites(cfg)
    s = SITES[0]
    today = today or (datetime.now(timezone.utc) + timedelta(hours=s['timezone_hours'])).date()
    for f in glob.glob(os.path.join(OUT, '*.json')):                 # old single-site layout
        os.remove(f)
    dates = wanted(cfg, today)
    keep = {d.strftime('%m-%d') for d in dates}
    jobs = []                                                          # (site, date)
    for st in SITES:
        od = os.path.join(OUT, st['key'])
        os.makedirs(od, exist_ok=True)
        fk = flight_key(cfg, st)
        for f in glob.glob(os.path.join(od, '*.json')):               # drop dates that left the window
            if os.path.basename(f)[:5] not in keep:
                os.remove(f)
        for d in dates:
            p = os.path.join(od, d.strftime('%m-%d') + '.json')
            if os.path.exists(p):
                old = json.load(open(p))
                if old.get('flight') == fk and old.get('target_year') == d.year:
                    continue
            jobs.append((st, d))
    if not jobs:
        print('history window up to date')
        return
    todo = sorted({d for _, d in jobs})
    print(f'{len(jobs)} site-dates to compute: {", ".join(d.strftime("%d %b") for d in todo)}')
    # which UTC days to download, grouped per (year, month)
    need = {}
    for d in todo:
        for y in range(d.year - NYEARS, d.year):
            try:
                dl = date(y, d.month, d.day)
            except ValueError:                                         # 29 Feb
                continue
            launch = sun.sunrise_utc(dl, s['lat'], s['lon'])
            for dd in (launch.date() - timedelta(days=0), launch.date() + timedelta(days=1)):
                need.setdefault((dd.year, dd.month), set()).add(dd.day)
    box = cfg['forecast'].get('era5_box') or cfg['forecast']['box']
    PL, SF = load_days(sorted(date(y, m, x) for (y, m), ds in need.items() for x in ds), box)
    for st, d in jobs:
        per = {}
        for y in range(d.year - NYEARS, d.year):
            try:
                dl = date(y, d.month, d.day)
                per[str(y)] = fly_year(cfg, st, dl, PL, SF)
            except ValueError:
                continue
            except Exception as e:
                print(f'    {dl}: failed ({e})')
                per[str(y)] = None
        out = dict(date=d.strftime('%m-%d'), target_year=d.year, flight=flight_key(cfg, st), site=st['key'],
                   source='ERA5 reanalysis', years=per, summary=summarise(per))
        with open(os.path.join(OUT, st['key'], d.strftime('%m-%d') + '.json'), 'w') as f:
            json.dump(out, f, indent=0)
        sm = out['summary']
        print(f'  {st["key"]} {d:%d %b}: GO in {sm.get("go")}/{sm.get("n")} years, landing median {sm.get("landing_median_km")} km')


if __name__ == '__main__':
    if not os.environ.get('CDSAPI_KEY') and not os.path.exists(os.path.expanduser('~/.cdsapirc')):
        print('No Copernicus key (CDSAPI_KEY secret) - history window skipped.')
    else:
        run(yaml.safe_load(open(os.path.join(ROOT, 'config.yaml'))))
