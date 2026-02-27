import fs from "fs";
import path from "path";
import { SimulationCore } from "../core/sim.mjs";

function parseArgs(argv) {
  const args = new Map();
  for (let i = 2; i < argv.length; i++) {
    const arg = argv[i];
    if (arg.startsWith("--")) {
      const key = arg.slice(2);
      const next = argv[i + 1];
      if (next && !next.startsWith("--")) {
        args.set(key, next);
        i++;
      } else {
        args.set(key, true);
      }
    }
  }
  return args;
}

function usage() {
  return [
    "Usage:",
    "  node js/adapters/cli.mjs --geojson data/Most_shoreline_polygon_UTM33.geojson --steps 10 --fish 5 --out out_tracks.csv",
    "",
    "Options:",
    "  --geojson <path>   Input GeoJSON file",
    "  --steps <n>        Number of dt_output steps to run (default: 1)",
    "  --fish <n>         Number of fish (default: 0)",
    "  --width <n>        Viewport width (default: 1200)",
    "  --height <n>       Viewport height (default: 800)",
    "  --out <path>       Output CSV file (default: stdout)",
    "  --help             Show this help",
  ].join("\n");
}

const args = parseArgs(process.argv);
if (args.has("help")) {
  console.log(usage());
  process.exit(0);
}

const geojsonPath = args.get("geojson");
if (!geojsonPath) {
  console.error("Missing --geojson path\n");
  console.error(usage());
  process.exit(1);
}

const steps = Number(args.get("steps") ?? 1);
const fishCount = Number(args.get("fish") ?? 0);
const width = Number(args.get("width") ?? 1200);
const height = Number(args.get("height") ?? 800);
const outPath = args.get("out");

const resolvedPath = path.resolve(process.cwd(), geojsonPath);
const geojson = JSON.parse(fs.readFileSync(resolvedPath, "utf-8"));

const sim = new SimulationCore();
sim.loadGeoJSON(geojson, { width, height });
if (fishCount > 0) sim.resetFish(fishCount);

const rows = [];
rows.push(["time", "fish_id", "x", "y", "z", "fish_state"].join(","));
for (let i = 0; i < steps; i++) {
  const snapshot = sim.step();
  snapshot.fish.forEach(fish => {
    const [x, y, z] = fish.position;
    rows.push([snapshot.time, fish.id, x, y, z, fish.state].join(","));
  });
}
const serialized = rows.join("\n");
if (outPath) {
  fs.writeFileSync(path.resolve(process.cwd(), outPath), serialized);
} else {
  console.log(serialized);
}
