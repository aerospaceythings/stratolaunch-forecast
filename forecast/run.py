"""Daily forecast run: python -m forecast.run [--cycle YYYYMMDDHH] [--days N]"""
import argparse, json, math, os
from datetime import datetime, timedelta, timezone, date
import numpy as np
import yaml
from scipy.interpolate import RegularGridInterpolator as RGI

from . import gfs, model, sun, site

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'data')


def plan_days(cfg, cycle, ndays, last_fh=384):
    s = cfg['site']
    tz = timedelta(hours=s['timezone_hours'])
    first = (cycle + tz).date() + timedelta(days=1)
    days = []
    for k in range(ndays + 1):
        d = first + timedelta(days=k)
        launch = sun.sunrise_utc(d, s['lat'], s['lon']) - timedelta(minutes=s['launch_minutes_before_sunrise'])
        fh = (launch - cycle).total_seconds() / 3600
        if fh < 0 or fh + 3 > last_fh:
            continue
        if fh + 3 <= 120:
            hours = list(range(int(math.floor(fh)), int(math.ceil(fh + 2.6)) + 1))
        else:
            hours = list(range(int(fh // 3 * 3), int(math.ceil((fh + 2.6) / 3) * 3) + 1, 3))
        days.append(dict(date=d, launch_utc=launch, fh=fh, hours=hours,
                         sunset_utc=sun.sunset_utc(d, s['lat'], s['lon'])))
    return days[:ndays]


def surface_at(sfc, hours, fh, lat, lon):
    """Time-interpolated surface values at the site at forecast hour fh, plus window maxima."""
    vals = {}
    names = set().union(*[set(sfc[h][2].keys()) for h in hours if sfc.get(h)]) if any(sfc.get(h) for h in hours) else set()
    for n in names:
        ts, xs = [], []
        for h in hours:
            if sfc.get(h) and n in sfc[h][2]:
                la, lo, f = sfc[h]
                xs.append(float(RGI((la, lo), f[n], bounds_error=False, fill_value=None)([[lat, lon]])[0]))
                ts.append(h)
        if xs:
            vals[n] = (float(np.interp(fh, ts, xs)), max(xs))
    return vals


def run(cfg, cycle=None, ndays=None):
    s, fl = cfg['site'], cfg['flight']
    ndays = ndays or cfg['forecast']['max_days']
    cycle = cycle or gfs.latest_cycle()
    print(f'GFS cycle {cycle:%Y-%m-%d %HZ}')
    days = plan_days(cfg, cycle, ndays)
    hours = sorted({h for d in days for h in d['hours']})
    print(f'{len(days)} launch days, {len(hours)} forecast hours to fetch')
    files = gfs.download(cycle, hours, cfg['forecast']['box'], os.path.join(ROOT, '.cache'))
    upper, sfc = {}, {}
    for h, (up, sf) in files.items():
        upper[h] = gfs.read_upper(up)
        sfc[h] = gfs.read_surface(sf)
    tz = timedelta(hours=s['timezone_hours'])
    results = []
    for d in days:
        hs = [h for h in d['hours'] if h in upper]
        if len(hs) < 2 or hs[0] > d['fh'] or hs[-1] < d['fh'] + 2.3:
            print(f'  {d["date"]}: not enough forecast data, skipped')
            continue
        la, lo = upper[hs[0]][:2]
        common = sorted(set.intersection(*[set(upper[h][2].tolist()) for h in hs]), reverse=True)
        idx = {h: [upper[h][2].tolist().index(lv) for lv in common] for h in hs}
        W = model.Winds(hs, la, lo, *[np.stack([upper[h][k][idx[h]] for h in hs]) for k in (3, 4, 5)])
        (llat, llon), (blat, blon), dur, path = model.fly(W, d['fh'], s['lat'], s['lon'], s['ground_m'], fl)
        L, brg = model.dist_km(s['lat'], s['lon'], llat, llon)
        B, _ = model.dist_km(s['lat'], s['lon'], blat, blon)
        H, u, v = W.column(d['fh'], s['lat'], s['lon'])
        pm = model.profile_metrics(H, u, v, cfg['float_mode']['band_km'], cfg['float_mode']['min_layer_speed_ms'])
        sv = surface_at(sfc, hs, d['fh'], s['lat'], s['lon'])
        g = lambda n, i=0: sv[n][i] if n in sv else None
        wind10 = None if g('u10') is None else math.hypot(g('u10'), g('v10'))
        landing_utc = d['launch_utc'] + timedelta(minutes=dur)
        r = dict(date=d['date'].isoformat(), lead_days=round(d['fh'] / 24, 1), fh=round(d['fh'], 4),
                 sunset_fh=round((d['sunset_utc'] - cycle).total_seconds() / 3600, 3),
                 launch_local=(d['launch_utc'] + tz).strftime('%H:%M'),
                 launch_utc=d['launch_utc'].strftime('%Y-%m-%dT%H:%MZ'),
                 sunrise_local=(d['launch_utc'] + timedelta(minutes=s['launch_minutes_before_sunrise']) + tz).strftime('%H:%M'),
                 landing_lat=round(llat, 4), landing_lon=round(llon, 4), landing_km=round(L, 1), landing_bearing=round(brg),
                 landing_dir=model.compass(brg), landing_over_ocean=model.over_ocean(llat, llon),
                 burst_lat=round(blat, 4), burst_lon=round(blon, 4), burst_km=round(B, 1), flight_min=round(dur),
                 daylight_after_landing_h=round((d['sunset_utc'] - landing_utc).total_seconds() / 3600, 1),
                 wind10_kmh=None if wind10 is None else round(wind10 * 3.6, 1),
                 gust_kmh=None if g('gust_ms') is None else round(g('gust_ms') * 3.6, 1),
                 vis_km=None if g('vis_m') is None else round(g('vis_m') / 1000, 1),
                 t2_c=None if g('t2_k') is None else round(g('t2_k') - 273.15, 1),
                 dewpoint_gap_c=None if g('t2_k') is None or g('td2_k') is None else round(g('t2_k') - g('td2_k'), 1),
                 low_cloud_pct=None if g('low_cloud_pct') is None else round(g('low_cloud_pct')),
                 total_cloud_pct=None if g('total_cloud_pct') is None else round(g('total_cloud_pct')),
                 rain_mmh=None if g('prate', 1) is None else round(g('prate', 1) * 3600, 2),
                 cape=None if g('cape', 1) is None else round(g('cape', 1)),
                 path=path)
        r['fog'] = model.fog_risk(g('vis_m'), g('t2_k'), g('td2_k'), wind10, g('low_cloud_pct'))
        r.update({k: (round(v, 2) if isinstance(v, float) else v) for k, v in pm.items()})
        r['burst_status'], r['burst_checks'] = model.burst_status(r, cfg['burst_mode'], cfg.get('avoid_zones'))
        r['float_status'], r['float_checks'] = model.float_status(pm, cfg['float_mode'])
        r['confidence'] = 'high' if r['lead_days'] <= 3 else 'medium' if r['lead_days'] <= 7 else 'low'
        results.append(r)
        print(f'  {r["date"]} {r["launch_local"]}  {r["burst_status"]:8s} {L:5.0f} km {r["landing_dir"]:3s}  '
              f'wind {r["wind10_kmh"]}  fog {r["fog"]}  float {r["float_status"]}')
    export_winds(upper, cycle)
    return dict(cycle=cycle.strftime('%Y-%m-%dT%HZ'), generated=datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%MZ'),
                model='NOAA GFS 0.25°', site=s, flight=fl, days=results)


def export_winds(upper, cycle, step=4):
    """Coarse (1 deg) copy of the winds for the website's what-if sliders -> site/winds.json"""
    import base64
    hs = sorted(upper)
    if not hs:
        return
    common = sorted(set.intersection(*[set(upper[h][2].tolist()) for h in hs]), reverse=True)
    la, lo = upper[hs[0]][0], upper[hs[0]][1]
    yi, xi = np.arange(0, len(la), step), np.arange(0, len(lo), step)
    def pack(k, scale):
        arr = np.stack([upper[h][k][[upper[h][2].tolist().index(lv) for lv in common]][:, yi][:, :, xi] for h in hs])
        return base64.b64encode(np.round(arr / scale).astype('<i2').tobytes()).decode()
    out = dict(cycle=cycle.strftime('%Y-%m-%dT%HZ'), hours=hs, lats=la[yi].round(3).tolist(), lons=lo[xi].round(3).tolist(),
               levels=common, shape=[len(hs), len(common), len(yi), len(xi)],
               u=pack(3, 0.05), v=pack(4, 0.05), gh=pack(5, 2.0), scale=dict(u=0.05, v=0.05, gh=2.0))
    os.makedirs(os.path.join(ROOT, 'site'), exist_ok=True)
    with open(os.path.join(ROOT, 'site', 'winds.json'), 'w') as f:
        json.dump(out, f, separators=(',', ':'))


def save(out):
    os.makedirs(os.path.join(DATA, 'runs'), exist_ok=True)
    tag = out['cycle'].replace('-', '').replace(':', '').replace('T', '')
    with open(os.path.join(DATA, 'runs', f'{tag}.json'), 'w') as f:
        json.dump(out, f)
    with open(os.path.join(DATA, 'latest.json'), 'w') as f:
        json.dump(out, f)
    hp = os.path.join(DATA, 'history.json')
    hist = json.load(open(hp)) if os.path.exists(hp) else {}
    keep = ('lead_days', 'burst_status', 'landing_km', 'landing_dir', 'landing_lat', 'landing_lon', 'fog',
            'wind10_kmh', 'jet_kmh', 'float_status', 'band_min_ms', 'steering_arc_deg')
    for r in out['days']:
        lst = [e for e in hist.get(r['date'], []) if e['run'] != out['cycle']]
        lst.append(dict(run=out['cycle'], **{k: r[k] for k in keep}))
        hist[r['date']] = sorted(lst, key=lambda e: e['run'])
    with open(hp, 'w') as f:
        json.dump(hist, f, indent=0)
    return hist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cycle', help='YYYYMMDDHH (default: latest complete GFS run)')
    ap.add_argument('--days', type=int)
    ap.add_argument('--site-only', action='store_true', help='rebuild the website from saved data')
    a = ap.parse_args()
    cfg = yaml.safe_load(open(os.path.join(ROOT, 'config.yaml')))
    if not a.site_only:
        cyc = datetime.strptime(a.cycle, '%Y%m%d%H').replace(tzinfo=timezone.utc) if a.cycle else None
        save(run(cfg, cyc, a.days))
    site.build(ROOT, cfg)


if __name__ == '__main__':
    main()
