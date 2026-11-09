# Stratolaunch forecast

A daily launch forecast for the balloon launch site near Northam, WA.

Every day at **13:20 AWST**, GitHub runs the forecast on its own servers and updates the site. For each of the next ~16 mornings it shows:

- **🎈 Burst flight:** a **GO / MARGINAL / NO-GO** rating for a launch 15 minutes before sunrise, with:
  - predicted landing spot and flight path
  - launch wind and gusts
  - fog risk
  - rain and storm risk
- **🛰️ Floating flight:** whether the 18–23.5 km band has a calm layer and a wide choice of wind directions. This is what a future floating platform would need.
- **📈 How each day's forecast has changed:** every run is saved, so you can watch a date's prediction settle as the day gets closer.

## How confident is it?

| Days ahead | Confidence | Use it for |
|---|---|---|
| 1–3 | high | final go/no-go, recovery planning |
| 4–7 | medium | choosing the launch day |
| 8–16 | low (shown faded) | spotting promising periods |

## Changing the rules

Edit **`config.yaml`** on GitHub (pencil icon) and commit. The site rebuilds automatically within a few minutes. You can change:
- the landing, wind, fog, rain and storm limits
- the climb rate, burst height and fall speed
- the launch time relative to sunrise
- the floating-band criteria
- places the balloon should not land (`avoid_zones`)

## Running it by hand

- **On GitHub:** **Actions → Daily launch forecast → Run workflow**.
- **Locally:**
  ```
  pip install -r requirements.txt
  python -m forecast.run            # latest GFS run, 16 days
  python -m forecast.run --site-only   # rebuild the page from saved data
  python -m tests.synthetic         # offline test with fake forecast data (needs: pip install eccodes)
  ```
  The page is written to `site/index.html`.

## Where the numbers come from

- **Forecast:** NOAA GFS 0.25°, the newest complete run, downloaded as a small area around the site from NOMADS. Pressure levels from 1000 to 3 hPa (surface to ~40 km), plus the surface fields (visibility, 2 m temperature and dew point, 10 m wind, gusts, precipitation rate, CAPE, cloud).
- **Trajectory:** ascent at a constant rate to burst height, then a parachute descent whose fall speed scales with 1/√(air density). Settings come from the 13 Oct 2026 flight.
- **Fog risk:** a rule of thumb. It's "high" when model visibility is under 1 km, or the air is saturated (dew-point gap ≤ 1 °C) and calm. It's "moderate" when visibility is under 5 km or the air is nearly saturated with light wind.
- **Landing over the ocean:** checked against a simplified WA coastline.

**Limits:**
- Single-model forecast with no ensemble yet.
- Simple balloon model.
- The fog rule hasn't yet been checked against observations.
- The website is public, so keep internal information out of this repo.

## Setup (one time)

**Settings → Pages → Build and deployment → Source: GitHub Actions.** Then run the workflow once from the Actions tab.
