import { Fish } from "./fish.mjs";
import { calculateBbox } from "./geo.mjs";
import { computeDistMatrixSparse } from "./legacy_physics.mjs";

export class SimulationCore {
  constructor(config = {}, rng = Math.random, hooks = {}) {
    this.config = {
      dt: 0.05,
      dt_output: 1,
      maxIter: 10000,
      depthResolution: 4,
      ...config,
    };
    this.rng = rng;
    this.hooks = hooks;
    this.time = 0;
    this.geojson = null;
    this.viewport = { width: 0, height: 0 };
    this.bbox = null;
    this.fishes = [];
    this.food_patches = [];
  }

  loadGeoJSON(geojson, viewport) {
    this.geojson = geojson;
    this.viewport = { ...this.viewport, ...viewport };
    this.bbox = calculateBbox(geojson);
  }

  resetFish(n) {
    this.fishes = [];
    for (let i = 0; i < n; i++) {
      const fish = new Fish(i);
      fish.initialize([0, 0, 0], 0);
      this.fishes.push(fish);
    }
  }

  step() {
    if (this.hooks.step) {
      return this.hooks.step(this);
    }
    // TODO: integrate physics + depth + social interactions.
    this.time += this.config.dt_output;
    return this.getState();
  }

  stepLegacy(env, polygonCoordinates) {
    const sim_steps = Math.floor(this.config.dt_output / this.config.dt);
    for (let i = 0; i < sim_steps; i++) {
      env.dist_matrix = computeDistMatrixSparse(env.fishes);
      env.fishes.forEach(fish => {
        fish.updatePosition(polygonCoordinates, this.config.dt);
      });
    }
    this.time += this.config.dt_output;
    return this.getState();
  }

  getState() {
    if (this.hooks.getState) {
      return this.hooks.getState(this);
    }
    return {
      time: this.time,
      fish: this.fishes.map(fish => ({
        id: fish.id,
        state: fish.state,
        position: fish.positions[fish.positions.length - 1] || [0, 0, 0],
      })),
    };
  }
}
