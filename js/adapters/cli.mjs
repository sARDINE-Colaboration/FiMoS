import fs from "fs";
import path from "path";
import { SimulationCore } from "../core/sim.mjs";
import { Fish, FoodPatch } from "../core/physics.mjs";
import {
  calculateDistanceToPolygon,
  getMatrixIndicesAndPoints,
  sampleDepthMap,
  mapDepth,
  triangulateGround,
  mapDepthPointIdxToTriangleStarts,
} from "../core/depth.mjs";
import {
  calculateBbox,
  normalizeToLineString,
  rescaleCoordinates2D,
  rescaleCoordinates3D,
} from "../core/geo.mjs";
import {
  calculateTransitionMatrix,
  drawStateParameterArray,
  drawStateParameterFromArray,
  stateSwitch,
} from "../core/state.mjs";
import { createDefaultStateConfig } from "../core/defaults.mjs";
import { Delaunay } from "d3-delaunay";

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
    "  --depth-res <n>    Depth resolution (default: 4)",
    "  --max-depth <n>    Max depth for shore-distance mode (default: -8)",
    "  --state-config <path>  JSON state config file with state_vectors (optional)",
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
const depthResolution = Number(args.get("depth-res") ?? 4);
const maxDepthArg = Number(args.get("max-depth") ?? -8);
const stateConfigPath = args.get("state-config");
const useStateConfig = Boolean(stateConfigPath);

const resolvedPath = path.resolve(process.cwd(), geojsonPath);
const geojson = JSON.parse(fs.readFileSync(resolvedPath, "utf-8"));

const sim = new SimulationCore({ dt: 0.05, dt_output: 1 });
sim.loadGeoJSON(geojson, { width, height });

if (geojson.features.length > 0) {
  normalizeToLineString(geojson.features[0]);
}

const bbox = calculateBbox(geojson);
const geojsonLimits = {
  minX: bbox[0],
  minY: bbox[1],
  maxX: bbox[2],
  maxY: bbox[3],
};

const rescaled2D = rescaleCoordinates2D(
  geojson.features[0].geometry.coordinates,
  { width, height },
  geojsonLimits
);
geojson.features[0].geometry.coordinates = rescaled2D.coords;

const rescaleState = {
  pixel_per_meter: rescaled2D.pixel_per_meter,
  x_as_max_extent: rescaled2D.x_as_max_extent,
  rescale_offset: rescaled2D.rescale_offset,
};

let maxDepth = maxDepthArg;
let depth_map = [];
let depth_map_points_idxs = [];
let depth_points = [];

if (geojson.features.length > 1 && geojson.features[1].geometry.type === "MultiPoint") {
  geojson.features[1].geometry.coordinates = rescaleCoordinates3D(
    geojson.features[1].geometry.coordinates,
    { width, height },
    rescaleState,
    geojsonLimits
  );
  depth_map = sampleDepthMap(geojson, {
    worldWidth: width,
    worldHeight: height,
    depthResolution,
  });
  [depth_map_points_idxs, depth_points] = getMatrixIndicesAndPoints(depth_map, depthResolution);
  const flatDepths = depth_map.flat().filter(d => !Number.isNaN(d));
  maxDepth = Math.min(...flatDepths);
} else {
  const depthMatrix = calculateDistanceToPolygon(geojson, {
    worldWidth: width,
    worldHeight: height,
    depthResolution,
  });
  depth_map = mapDepth(depthMatrix, maxDepth);
  [depth_map_points_idxs, depth_points] = getMatrixIndicesAndPoints(depth_map, depthResolution);
}

const shorelineCoords = geojson.features[0].geometry.coordinates;
const [triangulated_points, triangles] = triangulateGround(
  depth_points,
  shorelineCoords,
  depthResolution,
  Delaunay
);
const depth_point_idx_to_triangle_starts = mapDepthPointIdxToTriangleStarts(depth_points.length, triangles);

const stateDefaults = createDefaultStateConfig();
let state_means_matrix = stateDefaults.state_means_matrix;
let states_present = stateDefaults.states_present;
let substates = stateDefaults.substates;
let transition_probs = stateDefaults.transition_probs;
let sojourn_times = stateDefaults.sojourn_times;

if (stateConfigPath) {
  const resolvedStatePath = path.resolve(process.cwd(), stateConfigPath);
  const stateConfig = JSON.parse(fs.readFileSync(resolvedStatePath, "utf-8"));

  if (!Array.isArray(stateConfig.state_vectors) || stateConfig.state_vectors.length === 0) {
    throw new Error("state-config: state_vectors must be a non-empty array.");
  }

  const maxStateId = Math.max(...stateConfig.state_vectors.map(s => s.state_id));
  const expectedCount = maxStateId + 1;
  const idSet = new Set(stateConfig.state_vectors.map(s => s.state_id));
  if (idSet.size !== stateConfig.state_vectors.length) {
    throw new Error("state-config: state_id values must be unique.");
  }
  for (let i = 0; i < expectedCount; i++) {
    if (!idSet.has(i)) {
      throw new Error("state-config: state_id values must be contiguous starting at 0.");
    }
  }

  const ordered = stateConfig.state_vectors.slice().sort((a, b) => a.state_id - b.state_id);
  states_present = ordered.map(s => s.name);
  state_means_matrix = ordered.map(s => s.means);
  sojourn_times = ordered.map(s => s.sojourn_time);
  transition_probs = ordered.map(s => s.transition_probs);

  if (!transition_probs.every(row => Array.isArray(row) && row.length === ordered.length)) {
    throw new Error("state-config: each transition_probs must be an array sized to number of states.");
  }

  if (Array.isArray(stateConfig.substates)) {
    substates = stateConfig.substates;
  } else {
    substates = states_present.map(() => 1);
  }
  if (substates.length !== states_present.length) {
    throw new Error("state-config: substates length must match number of states.");
  }
}
let transition_matrix = calculateTransitionMatrix(transition_probs, sojourn_times, sim.config.dt_output);

const param = { globalStateStdev: 0.05 };

const fishes = [];
const food_patches = [];

const env = {
  viewport: { width, height },
  pixel_per_meter: rescaleState.pixel_per_meter,
  depth_map,
  depth_map_points_idxs,
  depth_points: triangulated_points,
  triangles,
  depth_point_idx_to_triangle_starts,
  triangles_to_draw: [],
  points_to_draw: [],
  fishes,
  food_patches,
  dist_matrix: [],
  depthResolution,
  geojson,
  food_patch_update_time: 1000,
  state_means_matrix,
  lower_limits: [0.1, 0.01, 0.01, 0.01, 0.01, 0, 0, 0],
  upper_limits: [3, 3, 1, 1, 1, 10, 5, 5],
  substates,
  param,
  states_present,
  transition_matrix,
  draw_state_parameter_array: () => drawStateParameterArray(states_present, env),
  draw_state_parameter_from_array: (arr, state) => drawStateParameterFromArray(arr, state, env),
  zero_vector: stateDefaults.zero_vector,
};

for (let i = 0; i < 5; i++) {
  food_patches.push(new FoodPatch(i + 1, env, geojson));
}

for (let i = 0; i < fishCount; i++) {
  fishes.push(new Fish(i, env, geojson));
}

if (useStateConfig && fishes.length > 0) {
  fishes.forEach(fish => {
    const randomState = Math.floor(Math.random() * states_present.length);
    fish.parameter_array = drawStateParameterArray(states_present, env);
    fish.state = randomState;
    const newParams = drawStateParameterFromArray(fish.parameter_array, fish.state, env);
    fish.beta = newParams[0];
    fish.v0 = newParams[1];
    fish.D_phi = newParams[2];
    fish.D_theta = newParams[3];
    fish.D_v = newParams[4];
    fish.patch_strength = newParams[5];
    fish.strength_att = newParams[6];
    fish.strength_align = newParams[7];
  });
}

const rows = [];
rows.push(["time", "fish_id", "x", "y", "z", "fish_state"].join(","));
for (let i = 0; i < steps; i++) {
  env.food_patches.forEach(patch => patch.updatePatch());
  const snapshot = sim.stepLegacy(env, geojson.features[0].geometry.coordinates);
  env.fishes.forEach(fish => {
    const [x, y, z] = fish.position;
    rows.push([snapshot.time, fish.id, x, y, z, fish.state].join(","));
  });
  env.fishes.forEach(fish => {
    stateSwitch(fish, env);
  });
}
const serialized = rows.join("\n");
if (outPath) {
  fs.writeFileSync(path.resolve(process.cwd(), outPath), serialized);
} else {
  console.log(serialized);
}
