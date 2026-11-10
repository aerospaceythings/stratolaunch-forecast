# Stratolaunch forecast

A daily launch forecast for two balloon launch sites in WA: **Northam Airport** and **Merredin Airport**. Switch between them at the top of the page.

Every day at **13:20 AWST**, GitHub runs the forecast on its own servers and updates the site. For each of the next ~16 mornings it shows:

- **🎈 Burst flight:** a **GO / MARGINAL / NO-GO** rating for a launch 15 minutes before sunrise, with:
  - predicted landing spot and flight path
  - launch wind and gusts
  - fog risk
  - rain and storm risk
- **🛰️ Floating flight:** the same dawn launch, but the balloon climbs to 20 km and drifts with the wind for up to the full 16-day forecast. There is no station-keeping. Every hour along the drift, the model asks: *if we cut down now, where would the payload land?* The site shows:
  - how long you can float before the **last safe cut-down**, meaning the payload still lands on land within 300 km
  - when that cut-down is (local time)
  - an hour-by-hour green/red cut-down timeline
  - the same at 18 and 22 km, because the real float height is never exact

  GOOD means at least 24 h of float, FAIR 6–24 h. Floating is forecast-only (no history cards).
- **📈 How each day's forecast has changed:** every run is saved, so you can watch a date's prediction settle as the day gets closer.

## No-landing zones

The landing (burst) or cut-down landing (floating) is **NO-GO** inside any of these zones, which are drawn in red on the map:
- **Perth metro and coastal suburbs:** from Guilderton to Mandurah/Dawesville, from the coast to the Darling Scarp.
- **Darling Range forest:** the Perth Hills, Mundaring, Avon Valley NP, Dwellingup and Collie.
- **Towns:** circles around the wheatbelt and regional towns, plus Bunbury and Kalgoorlie.

They are rough outlines. Check them on the map and adjust `avoid_zones` in `config.yaml`.

## Days 17–30: history window

Past the 16-day forecast, the site shows **history cards**: the same dawn launch flown through ERA5 reanalysis weather for the same calendar date in each of the previous 5 years. Each card shows how many years it would have been GO, plus the typical and range of landing distances. Click a card for the year-by-year list.

A separate daily job (`.github/workflows/history.yml`, 10:30 Perth time) downloads only the date that newly enters the window, so the window always reaches 30 days ahead. **It needs your Copernicus key:** in GitHub, **Settings → Secrets and variables → Actions → New repository secret**, name `CDSAPI_KEY`, value = the key from your `.cdsapirc` file (the part after `key:`). Without the secret the job skips quietly.

## Flight settings and burst calculator (on the website)

- **Flight settings:** sliders for ascent rate, landing speed and burst altitude. Every landing is recalculated in your browser from a coarser (1°) copy of the same forecast winds. Defaults come from `config.yaml` (5 m/s up, 5 m/s landing, 33,500 m burst), and nothing you try on the page is saved.
- **Burst altitude calculator:** pick a Kaymont balloon (1000–3000 g) or enter a custom one, enter the payload mass and gas, and fill by target ascent rate or neck lift. It gives burst altitude, ascent rate, neck lift and time to burst, using the same method as the CUSF burst calculator. **Use this…** pushes the result into the flight settings.
- **Map background:** switch between Dark / Streets / Satellite (top-right). If one provider is blocked, the next one loads automatically. To use your own tile service with an API key, set `map.custom_tiles_url` in `config.yaml`.

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
- the floating altitude, the extra heights to compare, the 300 km limit and the GOOD/FAIR hours (`float_mode`)
- places the balloon should not land (`avoid_zones`: polygons or circles)
- the launch sites (`site` and `other_sites`)

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
