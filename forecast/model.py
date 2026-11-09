"""Trajectory, launch-condition checks and floating-band metrics."""
import numpy as np
from scipy.interpolate import RegularGridInterpolator as RGI

G0, KM_DEG = 9.80665, 111.32
# rough WA coastline (lat, lon); west of it (or south of the last point) = ocean
COAST = np.array([(-28.5, 114.6), (-29.0, 114.9), (-29.5, 114.97), (-30.0, 115.0), (-30.5, 115.06),
                  (-31.0, 115.33), (-31.5, 115.62), (-32.0, 115.75), (-32.5, 115.73), (-33.0, 115.67),
                  (-33.3, 115.64), (-33.55, 115.35), (-34.0, 115.05), (-34.4, 115.15)])


def over_ocean(lat, lon):
    if lat < COAST[-1, 0]:
        return True
    return bool(lon < np.interp(lat, COAST[::-1, 0], COAST[::-1, 1]))


def dist_km(lat0, lon0, lat, lon):
    dy = (lat - lat0) * KM_DEG
    dx = (lon - lon0) * KM_DEG * np.cos(np.radians(lat0))
    return float(np.hypot(dx, dy)), float(np.degrees(np.arctan2(dx, dy)) % 360)


def compass(brg):
    return ['N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW'][int((brg + 11.25) // 22.5) % 16]


class Winds:
    """u, v, geopotential height columns, interpolated in time and space. times in hours (any origin)."""

    def __init__(self, times_h, lats, lons, U, V, GH):
        # U, V, GH: (T, L, ny, nx) -> (T, ny, nx, L, 3)
        arr = np.moveaxis(np.stack([U, V, GH], -1), 1, 3)
        self.t = np.asarray(times_h, float)
        self.f = RGI((self.t, lats, lons), arr, bounds_error=False, fill_value=None)

    def column(self, t, lat, lon):
        c = self.f([[t, lat, lon]])[0]
        o = np.argsort(c[:, 2])
        return c[o, 2], c[o, 0], c[o, 1]          # height m (ascending), u, v

    def at(self, t, lat, lon, h):
        H, u, v = self.column(t, lat, lon)
        return np.interp(h, H, u), np.interp(h, H, v)


def rho(h):
    return 1.225 * np.exp(-h / 7200.0)


def fly(W, t0, lat0, lon0, ground_m, fl, dt=20.0):
    """Ascent -> burst -> descent. Returns landing (lat, lon), burst (lat, lon), duration min, path [[lat,lon],...]"""
    lat, lon, h, t = lat0, lon0, ground_m, t0
    path, burst, k = [[lat, lon]], None, 0
    while True:
        u, v = W.at(t, lat, lon, h)
        lat += v * dt / 1000 / KM_DEG
        lon += u * dt / 1000 / (KM_DEG * np.cos(np.radians(lat)))
        t += dt / 3600
        k += 1
        if burst is None:
            h += fl['ascent_ms'] * dt
            if h >= fl['burst_km'] * 1000:
                burst = (lat, lon)
        else:
            h -= fl['descent_sea_level_ms'] * np.sqrt(1.225 / rho(h)) * dt
        if k % 15 == 0:
            path.append([round(lat, 4), round(lon, 4)])
        if burst is not None and h <= ground_m:
            path.append([round(lat, 4), round(lon, 4)])
            return (lat, lon), burst, (t - t0) * 60, path


def profile_metrics(H, u, v, band_km, min_layer):
    z = np.arange(0.5, 35.01, 0.25) * 1000
    uu, vv = np.interp(z, H, u), np.interp(z, H, v)
    sp = np.hypot(uu, vv)
    zk = z / 1000
    band = (zk >= band_km[0]) & (zk <= band_km[1])
    trop = (zk >= 5) & (zk <= 16)
    strat = (zk >= 15) & (zk <= 30)
    zs, us = zk[strat], uu[strat]
    cross = np.where(np.diff(np.sign(us)) != 0)[0]
    zero = float(zs[cross[0]] + 0.25 * us[cross[0]] / (us[cross[0]] - us[cross[0] + 1])) if len(cross) else None
    hd = np.sort(np.degrees(np.arctan2(uu[band], vv[band]))[sp[band] > min_layer] % 360)
    arc = 0.0 if len(hd) < 2 else float(360 - np.diff(np.r_[hd, hd[0] + 360]).max())
    return dict(jet_kmh=float(sp[trop].max() * 3.6), jet_alt_km=float(zk[trop][sp[trop].argmax()]),
                zero_wind_km=zero, band_min_ms=float(sp[band].min()),
                band_min_alt_km=float(zk[band][sp[band].argmin()]), steering_arc_deg=arc,
                profile=dict(z_km=zk[::4].round(2).tolist(), u_kmh=(uu[::4] * 3.6).round(1).tolist(),
                             v_kmh=(vv[::4] * 3.6).round(1).tolist()))


def fog_risk(vis_m, t2_k, td2_k, wind10_ms, low_cloud_pct):
    """Heuristic fog risk at dawn: model visibility, or saturated + calm air near the ground."""
    dd = None if t2_k is None or td2_k is None else t2_k - td2_k
    if vis_m is not None and vis_m < 1000:
        return 'high'
    if dd is not None and dd <= 1.0 and wind10_ms is not None and wind10_ms < 2.5:
        return 'high'
    if vis_m is not None and vis_m < 5000:
        return 'moderate'
    if dd is not None and dd <= 2.5 and wind10_ms is not None and wind10_ms < 4:
        return 'moderate'
    return 'low'


def burst_status(r, c, avoid):
    """r: day result dict; c: burst_mode config. Returns status, list of (check, level, text)."""
    checks = []

    def add(name, level, text):
        checks.append(dict(check=name, level=level, text=text))

    L = r['landing_km']
    if r['landing_over_ocean']:
        add('landing', 'nogo', 'lands over the ocean')
    elif L > c['marginal_landing_km']:
        add('landing', 'nogo', f'lands {L:.0f} km away')
    elif L > c['max_landing_km']:
        add('landing', 'marginal', f'lands {L:.0f} km away')
    else:
        add('landing', 'go', f'{L:.0f} km {compass(r["landing_bearing"])}')
    for z in avoid or []:
        if dist_km(z['lat'], z['lon'], r['landing_lat'], r['landing_lon'])[0] < z['radius_km']:
            add('zone', 'nogo', f'lands in {z["name"]}')
    w = r['wind10_kmh']
    if w is None:
        add('wind', 'marginal', 'no surface wind data')
    else:
        add('wind', 'go' if w < c['max_launch_wind_kmh'] else 'marginal' if w < c['marginal_launch_wind_kmh'] else 'nogo',
            f'launch wind {w:.0f} km/h')
    g = r.get('gust_kmh')
    if g is not None and g > c['max_gust_kmh']:
        add('gust', 'marginal', f'gusts {g:.0f} km/h')
    f = r['fog']
    allowed = ['low'] if c['fog_allowed'] == 'low' else ['low', 'moderate']
    add('fog', 'go' if f in allowed else 'marginal' if f == 'moderate' else 'nogo', f'fog risk {f}')
    rain = r['rain_mmh']
    if rain is not None:
        add('rain', 'go' if rain < c['max_rain_mmh'] else 'marginal' if rain < c['marginal_rain_mmh'] else 'nogo',
            'no rain' if rain < c['max_rain_mmh'] else f'rain {rain:.1f} mm/h')
    cape = r['cape']
    if cape is not None:
        add('storm', 'go' if cape < c['max_cape'] else 'marginal' if cape < c['marginal_cape'] else 'nogo',
            'no storm risk' if cape < c['max_cape'] else f'storm risk (CAPE {cape:.0f})')
    if r['daylight_after_landing_h'] < c['min_daylight_after_landing_h']:
        add('daylight', 'nogo', 'lands too close to sunset')
    lv = [x['level'] for x in checks]
    status = 'nogo' if 'nogo' in lv else 'marginal' if 'marginal' in lv else 'go'
    return status, checks


def float_status(m, c):
    calm = m['band_min_ms'] <= c['max_band_min_wind_ms']
    arc = m['steering_arc_deg'] >= c['min_steering_arc_deg']
    return ('good' if calm and arc else 'fair' if calm or arc else 'poor'), dict(calm_layer=calm, steering=arc)
