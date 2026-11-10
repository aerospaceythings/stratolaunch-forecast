"""Daily forecast run: python -m forecast.run [--cycle YYYYMMDDHH] [--days N]"""
import argparse, json, math, os
from datetime import datetime, timedelta, timezone, date
import numpy as np
import yaml
from scipy.interpolate import RegularGridInterpolator as RGI

from . import gfs, model, sun, site

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'data')


def plan_days(s, cycle, ndays, last_fh=384):
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
    fl = cfg['flight']
    ndays = ndays or cfg['forecast']['max_days']
    cycle = cycle or gfs.latest_cycle()
    print(f'GFS cycle {cycle:%Y-%m-%d %HZ}')
    SITES = model.sites(cfg)
    plans = {s['key']: plan_days(s, cycle, ndays) for s in SITES}
    hours = sorted({h for ds in plans.values() for d in ds for h in d['hours']})
    print(f'{len(SITES)} sites, {len(hours)} forecast hours to fetch')
    files = gfs.download(cycle, hours, cfg['forecast']['box'], os.path.join(ROOT, '.cache'))
    upper, sfc = {}, {}
    for h, (up, sf) in files.items():
        upper[h] = gfs.read_upper(up)
        sfc[h] = gfs.read_surface(sf)
    out_sites = {}
    for s in SITES:
        print(f'{s["name"]}:')
        out_sites[s['key']] = dict(site=s, days=fly_site(cfg, s, cycle, plans[s['key']], upper, sfc))
    export_winds(upper, cycle)
    floating(cfg, cycle, [(s, plans[s['key']], out_sites[s['key']]['days']) for s in SITES])
    return dict(cycle=cycle.strftime('%Y-%m-%dT%HZ'), generated=datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%MZ'),
                model='NOAA GFS 0.25°', flight=fl, sites=out_sites)


def fly_site(cfg, s, cycle, days, upper, sfc):
    fl = cfg['flight']
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
        r['confidence'] = 'high' if r['lead_days'] <= 3 else 'mid' if r['lead_days'] <= 7 else 'low'
        results.append(r)
        print(f'  {r["date"]} {r["launch_local"]}  {r["burst_status"]:8s} {L:5.0f} km {r["landing_dir"]:3s}  '
              f'wind {r["wind10_kmh"]}  fog {r["fog"]}')
    return results


def floating(cfg, cycle, per_site):
    """Free drift at the float altitude(s) from each dawn launch to the end of the forecast; when to cut down."""
    fc, fl = cfg['float_mode'], cfg['flight']
    per_site = [x for x in per_site if x[2]]
    if not per_site or not fc.get('float_km'):
        return
    fb = cfg['forecast'].get('float_box') or dict(north=-22, south=-41, west=105, east=128)
    first = int(min(d['fh'] for _, days, _ in per_site for d in days) // 3 * 3)
    hours = list(range(first, 385, 3))
    print(f'floating: {len(hours)} forecast hours over the wide area')
    files = gfs.download_float(cycle, hours, fb, os.path.join(ROOT, '.cache'))
    up = {h: gfs.read_upper(p, stride=cfg['forecast'].get('float_stride', 2)) for h, p in files.items()}
    hs = sorted(up)
    if len(hs) < 3:
        print('  floating: not enough data, skipped')
        return
    la, lo = up[hs[0]][:2]
    common = sorted(set.intersection(*[set(up[h][2].tolist()) for h in hs]), reverse=True)
    idx = {h: [up[h][2].tolist().index(lv) for lv in common] for h in hs}
    U, V, GH = [np.stack([up[h][k][idx[h]] for h in hs]) for k in (3, 4, 5)]
    W = model.Winds(hs, la, lo, U, V, GH)
    alts = [fc['float_km']] + [a for a in fc.get('sensitivity_km', []) if a != fc['float_km']]
    LVs = {a: model.Level(W, hs, la, lo, U, V, GH, a * 1000) for a in alts}
    for s, _, results in per_site:
        tz = timedelta(hours=s['timezone_hours'])
        for r in results:
            out = {a: model.float_flight(W, LVs[a], r['fh'], s, fl, dict(fc, _h=a * 1000), cfg.get('avoid_zones'))
                   for a in alts}
            launch = datetime.strptime(r['launch_utc'], '%Y-%m-%dT%H:%MZ').replace(tzinfo=timezone.utc)
            for f in out.values():
                if f['last_cut']:
                    f['last_cut']['local'] = (launch + timedelta(hours=f['last_cut']['h']) + tz).strftime('%a %d %b %H:%M')
            best = max(alts, key=lambda a: (out[a]['safe_h'], a == fc['float_km']))     # longest float; ties -> main
            m = dict(out[best], main_km=fc['float_km'])
            m['others'] = [dict(float_km=a, safe_h=out[a]['safe_h'], window_end=out[a]['window_end'],
                                open_end=out[a]['open_end'], path=out[a]['path'][::2], last_cut=out[a]['last_cut'])
                           for a in sorted(alts) if a != best]
            r['drift'] = m
            r['float_safe_h'] = m['safe_h']
            r['float_best_km'] = best
            wc = model.weather_checks(r, cfg['burst_mode'])                # same launch-weather limits as burst
            r['float_weather'] = model.worst(wc)
            r['float_weather_checks'] = wc
            r['float_status'] = model.float_rating(m['safe_h'], fc, r['float_weather'])
            print(f'  {s["key"]} {r["date"]} best {best} km: safe {m["safe_h"]} h ({m["window_end"]}), '
                  f'weather {r["float_weather"]}; ' + ', '.join(f'{a} km {out[a]["safe_h"]} h' for a in alts))

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
    """data/latest.json (all sites), data/runs/<cycle>.json (no drift paths), data/history.json {site: {date: [...]}}"""
    os.makedirs(os.path.join(DATA, 'runs'), exist_ok=True)
    tag = out['cycle'].replace('-', '').replace(':', '').replace('T', '')
    slim = {k: dict(v, days=[{a: b for a, b in r.items() if a != 'drift'} for r in v['days']]) for k, v in out['sites'].items()}
    with open(os.path.join(DATA, 'runs', f'{tag}.json'), 'w') as f:
        json.dump(dict(out, sites=slim), f)
    with open(os.path.join(DATA, 'latest.json'), 'w') as f:
        json.dump(out, f)
    hp = os.path.join(DATA, 'history.json')
    hist = json.load(open(hp)) if os.path.exists(hp) else {}
    if hist and all(isinstance(v, list) for v in hist.values()):      # old single-site format
        hist = {'northam': hist}
    keep = ('lead_days', 'burst_status', 'landing_km', 'landing_dir', 'landing_lat', 'landing_lon', 'fog',
            'wind10_kmh', 'jet_kmh', 'float_status', 'band_min_ms', 'steering_arc_deg', 'float_safe_h', 'float_best_km')
    for key, v in out['sites'].items():
        h = hist.setdefault(key, {})
        for r in v['days']:
            lst = [e for e in h.get(r['date'], []) if e['run'] != out['cycle']]
            lst.append(dict(run=out['cycle'], **{k: r.get(k) for k in keep}))
            h[r['date']] = sorted(lst, key=lambda e: e['run'])
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
