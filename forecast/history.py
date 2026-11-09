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


def fetch(c, dataset, req, path):
    if os.path.exists(path):
        return path
    os.makedirs(os.path.dirname(path), exist_ok=True)
    c.retrieve(dataset, req, path + '.part')
    os.replace(path + '.part', path)
    return path


def open_nc(path):
    import xarray as xr
    ds = xr.open_dataset(path)
    ren = {a: b for a, b in [('valid_time', 'time'), ('pressure_level', 'level'), ('latitude', 'lat'),
                             ('longitude', 'lon')] if a in ds.dims or a in ds.coords}
    ds = ds.rename(ren)
    ds = ds.drop_vars([v for v in ds.coords if v not in ('time', 'level', 'lat', 'lon')])
    return ds.sortby('lat').sortby('lon')


def flight_key(cfg):
    f = cfg['flight']
    return f"{f['ascent_ms']}/{f['descent_sea_level_ms']}/{f['burst_km']}"


def wanted(cfg, today):
    """Local calendar dates in the window (today + 17 .. today + 30)."""
    return [today + timedelta(days=k) for k in range(FIRST_DAY, LAST_DAY + 1)]


def fly_year(cfg, d_local, pl, sf):
    """One historical launch on the same calendar date in another year."""
    s, fl = cfg['site'], cfg['flight']
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
    s = cfg['site']
    today = today or (datetime.now(timezone.utc) + timedelta(hours=s['timezone_hours'])).date()
    os.makedirs(OUT, exist_ok=True)
    fk = flight_key(cfg)
    dates = wanted(cfg, today)
    keep = {d.strftime('%m-%d') for d in dates}
    for f in glob.glob(os.path.join(OUT, '*.json')):                 # drop dates that left the window
        if os.path.basename(f)[:5] not in keep:
            os.remove(f)
    todo = []
    for d in dates:
        p = os.path.join(OUT, d.strftime('%m-%d') + '.json')
        if os.path.exists(p):
            old = json.load(open(p))
            if old.get('flight') == fk and old.get('target_year') == d.year:
                continue
        todo.append(d)
    if not todo:
        print('history window up to date')
        return
    print(f'{len(todo)} new dates: {", ".join(d.strftime("%d %b") for d in todo)}')
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
    box = cfg['forecast']['box']
    area = [box['north'], box['west'], box['south'], box['east']]
    c = client()
    files = {}
    for (y, m), days in sorted(need.items()):
        tag = f'{y}{m:02d}_' + '-'.join(f'{x:02d}' for x in sorted(days))
        base = {'product_type': ['reanalysis'], 'year': [str(y)], 'month': [f'{m:02d}'],
                'day': [f'{x:02d}' for x in sorted(days)], 'time': TIMES, 'area': area,
                'data_format': 'netcdf', 'download_format': 'unarchived'}
        print(f'  ERA5 {y}-{m:02d} days {sorted(days)}', flush=True)
        pl = fetch(c, 'reanalysis-era5-pressure-levels',
                   dict(base, variable=['u_component_of_wind', 'v_component_of_wind', 'geopotential'],
                        pressure_level=LEVELS), os.path.join(CACHE, f'pl_{tag}.nc'))
        sf = fetch(c, 'reanalysis-era5-single-levels', dict(base, variable=SFC_VARS, time=[f'{h:02d}:00' for h in range(24)]),
                   os.path.join(CACHE, f'sfc_{tag}.nc'))
        files[(y, m)] = (pl, sf)
    import xarray as xr
    PL = xr.concat([open_nc(p)[['u', 'v', 'z']] for p, _ in files.values()], 'time').sortby('time')
    SF = xr.concat([open_nc(f) for _, f in files.values()], 'time').sortby('time')
    PL = PL.isel(time=np.unique(PL['time'].values, return_index=True)[1])
    SF = SF.isel(time=np.unique(SF['time'].values, return_index=True)[1])
    for d in todo:
        per = {}
        for y in range(d.year - NYEARS, d.year):
            try:
                dl = date(y, d.month, d.day)
                per[str(y)] = fly_year(cfg, dl, PL, SF)
            except ValueError:
                continue
            except Exception as e:
                print(f'    {dl}: failed ({e})')
                per[str(y)] = None
        out = dict(date=d.strftime('%m-%d'), target_year=d.year, flight=fk, source='ERA5 reanalysis',
                   years=per, summary=summarise(per))
        with open(os.path.join(OUT, d.strftime('%m-%d') + '.json'), 'w') as f:
            json.dump(out, f, indent=0)
        sm = out['summary']
        print(f'  {d:%d %b}: GO in {sm.get("go")}/{sm.get("n")} years, landing median {sm.get("landing_median_km")} km')


if __name__ == '__main__':
    if not os.environ.get('CDSAPI_KEY') and not os.path.exists(os.path.expanduser('~/.cdsapirc')):
        print('No Copernicus key (CDSAPI_KEY secret) - history window skipped.')
    else:
        run(yaml.safe_load(open(os.path.join(ROOT, 'config.yaml'))))
