# FiMoS

Minimal instructions to install and run the simulation.

## Install

From the repo root:

```bash
npm install
```

## Run Simulation (No GUI)

Example (writes CSV):

```bash
node js/adapters/cli.mjs \
  --geojson data/Most_shoreline_polygon_UTM33.geojson \
  --steps 5 \
  --fish 10 \
  --out out_tracks.csv
```

Output format: `time,fish_id,x,y,z,fish_state`.

## Run GUI (Local Server)

Start a static server from the repo root:

```bash
npx http-server .
```

Then open:

```
http://127.0.0.1:8080/Behaviour_simulater.html
```

If the GeoJSON fails to load, hard-reload the page to clear any cached script files.
