"""Sunrise / sunset (NOAA solar calculator equations), UTC."""
import math
from datetime import datetime, timedelta, timezone


def _sun_event(date, lat, lon, rising, zenith=90.833):
    n = date.timetuple().tm_yday
    g = 2 * math.pi / 365 * (n - 1)
    eqt = 229.18 * (0.000075 + 0.001868 * math.cos(g) - 0.032077 * math.sin(g)
                    - 0.014615 * math.cos(2 * g) - 0.040849 * math.sin(2 * g))
    dec = (0.006918 - 0.399912 * math.cos(g) + 0.070257 * math.sin(g) - 0.006758 * math.cos(2 * g)
           + 0.000907 * math.sin(2 * g) - 0.002697 * math.cos(3 * g) + 0.00148 * math.sin(3 * g))
    la = math.radians(lat)
    ha = math.degrees(math.acos(math.cos(math.radians(zenith)) / (math.cos(la) * math.cos(dec))
                                - math.tan(la) * math.tan(dec)))
    minutes = 720 - 4 * (lon + (ha if rising else -ha)) - eqt
    return datetime(date.year, date.month, date.day, tzinfo=timezone.utc) + timedelta(minutes=minutes)


def sunrise_utc(local_date, lat, lon):
    """Sunrise (UTC datetime) on the given local calendar date (east of Greenwich, so it may fall on the previous UTC day)."""
    return _sun_event(local_date, lat, lon, True)


def sunset_utc(local_date, lat, lon):
    return _sun_event(local_date, lat, lon, False)
