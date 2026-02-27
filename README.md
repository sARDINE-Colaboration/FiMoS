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

### Optional State Config

You can provide a JSON file with arbitrary states:

```json
{
  "state_vectors": [
    {
      "name": "resting",
      "means": [0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
      "sojourn_time": 20,
      "transition_probs": [0, 0.5, 0.5],
      "state_id": 0
    },
    {
      "name": "foraging",
      "means": [0.1, 0.4, 0.9, 0.9, 0.8, 0.0, 0.0, 0.0],
      "sojourn_time": 20,
      "transition_probs": [0.5, 0, 0.5],
      "state_id": 1
    },
    {
      "name": "active",
      "means": [0.8, 0.8, 0.1, 0.1, 0.2, 0.0, 0.0, 0.0],
      "sojourn_time": 20,
      "transition_probs": [0.5, 0.5, 0],
      "state_id": 2
    }
  ],
  "substates": [2, 1, 1]
}
```

Use it like:

```bash
node js/adapters/cli.mjs \
  --geojson data/Most_shoreline_polygon_UTM33.geojson \
  --steps 5 \
  --fish 10 \
  --state-config input_state.json \
  --out out_tracks.csv
```

## Run GUI (Local Server)

Start a static server from the repo root:

```bash
npx http-server .
```

Then open:

```
http://127.0.0.1:8080/Behaviour_simulater.html
```

If the GeoJSON fails to load, hard-reload (firefox: cmd+shift+r) the page to clear any cached script files.
