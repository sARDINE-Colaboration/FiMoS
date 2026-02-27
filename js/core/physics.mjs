import { Vector } from "./vector.mjs";
import {
  crossProduct,
  dotProduct,
  depthIndicesWithConflicts,
  getTrianglesFromPointIndices,
  getTriangleCoords,
  lineIntersectsTriangle_ppk,
} from "./collision.mjs";

export function randomNormal(mean, sd) {
  const u = 1 - Math.random();
  const v = Math.random();
  const num = Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
  return num * sd + mean;
}

export function pointInsidePolygon(point, polygon) {
  const x = point[0];
  const y = point[1];
  let inside = false;
  for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
    const xi = polygon[i][0], yi = polygon[i][1];
    const xj = polygon[j][0], yj = polygon[j][1];
    const intersect = ((yi > y) !== (yj > y)) &&
      (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
    if (intersect) inside = !inside;
  }
  return inside;
}

export function lineSegmentIntersection(p1, p2, p3, p4) {
  const x1 = p1[0], y1 = p1[1];
  const x2 = p2[0], y2 = p2[1];
  const x3 = p3[0], y3 = p3[1];
  const x4 = p4[0], y4 = p4[1];
  let nominator = (x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4);
  const denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
  const t = nominator / denominator;
  if (t >= 0 && t <= 1) {
    nominator = (x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3);
    const u = -nominator / denominator;
    if (u >= 0 && u <= 1) {
      return [t, x1 + t * (x2 - x1), y1 + t * (y2 - y1)];
    }
  }
  return false;
}

export function computeDistMatrixSparse(fishes) {
  const dist_matrix = [];
  for (let i = 0; i < fishes.length; i++) {
    const row = [];
    for (let j = i + 1; j < fishes.length; j++) {
      const distance = fishes[i].position.sub(fishes[j].position).norm();
      row.push(distance);
    }
    let row_first = [];
    for (let k = 0; k < i; k++) {
      row_first.push(dist_matrix[k][i]);
    }
    row_first = row_first.concat(0);
    dist_matrix.push(row_first.concat(row));
  }
  return dist_matrix;
}

export function closestCollisionWithTriangle(p, p_next, triangles_as_indices, depth_points, triangles_to_draw) {
  let t_min = 1.2;
  let norm_min = null;
  let triangle_min = null;
  for (let i = 0; i < triangles_as_indices.length; i++) {
    const triangle = getTriangleCoords(triangles_as_indices[i], depth_points);
    const [t, norm] = lineIntersectsTriangle_ppk(p, p_next, triangle);
    if (t < t_min) {
      t_min = t;
      norm_min = norm;
      triangle_min = triangle;
    }
  }
  if (t_min <= 1) {
    if (triangles_to_draw) triangles_to_draw.push(triangle_min);
    norm_min = new Vector(norm_min[0], norm_min[1], norm_min[2]);
    norm_min = norm_min.div_scalar(norm_min.norm());
    if (norm_min[2] < 0) {
      norm_min = norm_min.mul_scalar(-1);
    }
  }
  return [t_min, norm_min, triangle_min];
}

export function checkIfBelowGround(p, env) {
  const indices2D = depthIndicesWithConflicts(p, p, env.depth_map, env.depthResolution);
  const indices = [...new Set(
    indices2D.map(index => env.depth_map_points_idxs[index[0]][index[1]])
  )];
  const triangles_as_indices = getTrianglesFromPointIndices(indices, env.triangles, env.depth_point_idx_to_triangle_starts);
  const p_surface = [p[0], p[1], 0];
  const [t_min, norm_min] = closestCollisionWithTriangle(p, p_surface, triangles_as_indices, env.depth_points, env.triangles_to_draw);
  return [t_min, norm_min];
}

export function surfaceRepulsionForce(fish, response_time, strength) {
  const p = fish.position, v = fish.velocity, phi = fish.phi, speed_ms = fish.speed_ms;
  const p_after_respTime = p.add(v.mul_scalar(response_time));
  const jump = p_after_respTime[2];
  let force_theta = crossProduct([
    Math.cos(phi + Math.PI / 2),
    Math.sin(phi + Math.PI / 2),
    0,
  ], v);
  force_theta = new Vector(force_theta[0], force_theta[1], force_theta[2]);
  if (speed_ms !== 0) {
    force_theta = force_theta.div_scalar(force_theta.norm());
    if (force_theta[2] > 0) force_theta = force_theta.mul_scalar(-1);
  }
  if (jump < 0) {
    force_theta[0] = 0; force_theta[1] = 0; force_theta[2] = 0;
  } else {
    force_theta[2] -= 1;
    force_theta = force_theta.div_scalar(force_theta.norm());
    let strength_scaled = strength;
    if (p[2] < 0) {
      const depth = Math.abs(p[2]);
      strength_scaled = strength * (1 - depth / (jump + depth));
    }
    force_theta = force_theta.mul_scalar(strength_scaled);
  }
  return force_theta;
}

export function groundRepulsionForce(fish, response_time, strength, env) {
  const p = fish.position, v = fish.velocity, phi = fish.phi, speed_ms = fish.speed_ms;
  let force_theta = crossProduct([
    Math.cos(phi + Math.PI / 2),
    Math.sin(phi + Math.PI / 2),
    0,
  ], v);
  force_theta = new Vector(force_theta[0], force_theta[1], force_theta[2]);
  if (speed_ms !== 0) {
    force_theta = force_theta.div_scalar(force_theta.norm());
    if (force_theta[2] < 0) force_theta = force_theta.mul_scalar(-1);
  }
  const [t_below_ground, norm_below_ground] = checkIfBelowGround(p, env);
  if (-0.01 < t_below_ground && t_below_ground < 1.01) {
    force_theta = force_theta.add(norm_below_ground);
    force_theta = force_theta.div_scalar(force_theta.norm());
    return force_theta.mul_scalar(strength);
  }

  const p_after_respTime = p.add(v.mul_scalar(response_time));
  const indices2D = depthIndicesWithConflicts(p, p_after_respTime, env.depth_map, env.depthResolution);
  let force = new Vector(0, 0, 0);
  if (indices2D.length === 0 || fish.speed_ms === 0) {
    return force;
  }
  const indices = [...new Set(
    indices2D.map(index => env.depth_map_points_idxs[index[0]][index[1]])
  )];
  const triangles_as_indices = getTrianglesFromPointIndices(indices, env.triangles, env.depth_point_idx_to_triangle_starts);
  const [t_min, norm_min] = closestCollisionWithTriangle(p, p_after_respTime, triangles_as_indices, env.depth_points, env.triangles_to_draw);
  if (t_min > 1) {
    fish.ground_avoidance_sign = 0;
    return force;
  }
  const phi_ground = norm_min.phi();
  let phi_force = 0;
  if (fish.ground_avoidance_sign !== 0 || fish.shore_avoidance_sign !== 0) {
    if (fish.ground_avoidance_sign === 0) {
      fish.ground_avoidance_sign = fish.shore_avoidance_sign;
    }
    phi_force = phi + fish.ground_avoidance_sign * Math.PI / 2;
  } else {
    phi_force = phi + Math.PI / 2;
    fish.ground_avoidance_sign = 1;
    if (Math.cos(phi_force) * Math.cos(phi_ground) + Math.sin(phi_force) * Math.sin(phi_ground) < 0) {
      phi_force += Math.PI;
      fish.ground_avoidance_sign = -1;
    }
  }
  void phi_force; // force direction is computed via force_theta below
  force_theta = force_theta.add(norm_min);
  force_theta = force_theta.div_scalar(force_theta.norm());
  force = force_theta.mul_scalar(strength * (1 - t_min));
  return force;
}

export function shoreRepulsionForce(fish, polygon, response_time, strength) {
  const p = fish.position, v = fish.velocity, phi = fish.phi;
  const p_after_respTime = p.add(v.mul_scalar(response_time));
  let force = new Vector(0, 0, 0);
  let intersection_distance_in_t = 1.1;
  let index_line_segment = null;
  for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
    const intersection = lineSegmentIntersection(p, p_after_respTime, polygon[i], polygon[j]);
    if (intersection !== false) {
      if (intersection[0] < intersection_distance_in_t) {
        intersection_distance_in_t = intersection[0];
        index_line_segment = j;
      }
    }
  }
  if (index_line_segment !== null) {
    const j = index_line_segment;
    const p1 = polygon[j], p2 = polygon[j + 1];
    const vel_seg = new Vector(p1[0] - p2[0], p1[1] - p2[1]);
    let phi_force = vel_seg.phi() + Math.PI / 2;
    if (Math.cos(phi_force) * Math.cos(phi) + Math.sin(phi_force) * Math.sin(phi) > 0) {
      phi_force += Math.PI;
    }
    let phi_force2;
    if (fish.shore_avoidance_sign !== 0 || fish.ground_avoidance_sign !== 0) {
      if (fish.shore_avoidance_sign === 0) {
        fish.shore_avoidance_sign = fish.ground_avoidance_sign;
      }
      phi_force2 = fish.phi + fish.shore_avoidance_sign * Math.PI / 2;
    } else {
      phi_force2 = fish.phi + Math.PI / 2;
      fish.shore_avoidance_sign = 1;
      if (Math.cos(phi_force) * Math.cos(phi_force2) + Math.sin(phi_force) * Math.sin(phi_force2) < 0) {
        phi_force2 += Math.PI;
        fish.shore_avoidance_sign = -1;
      }
    }
    force[0] = Math.cos(phi_force2);
    force[1] = Math.sin(phi_force2);
    force = force.mul_scalar(strength * (1 - intersection_distance_in_t));
  } else {
    fish.shore_avoidance_sign = 0;
  }
  return force;
}

export function patchAttractionForce(fish, patch_sensing_length, strength, env) {
  const p = fish.position;
  let force = new Vector(0, 0, 0);
  let dist_threshold = patch_sensing_length * env.pixel_per_meter;
  let index_patch = null;

  env.food_patches.forEach((patch, index) => {
    const patch_dist = p.sub(patch.position).norm();
    if (patch_dist < dist_threshold) {
      dist_threshold = patch_dist;
      index_patch = index;
    }
  });

  if (index_patch !== null) {
    const closest_patch = env.food_patches[index_patch];
    let direction_to_patch = closest_patch.position.sub(p);
    const magnitude = direction_to_patch.norm();
    if (magnitude > 0) {
      direction_to_patch = direction_to_patch.div_scalar(magnitude);
    }
    force = direction_to_patch.mul_scalar(strength);
  }
  return force;
}

export function socialForce(fish, strength_att, strength_align, env) {
  const p = fish.position, v_meter = fish.velocity.div_scalar(env.pixel_per_meter);
  const t_min = 2;
  const t_max = 10;
  const d_min = t_min * fish.speed_ms * env.pixel_per_meter;
  const d_max = t_max * fish.speed_ms * env.pixel_per_meter;
  const dist_to_others = env.dist_matrix[fish.id] || [];
  let force_rep = new Vector(0, 0, 0);
  let force_att = new Vector(0, 0, 0);
  let force_ali = new Vector(0, 0, 0);
  let counter_rep = 0;
  fish.social_nn = [];
  for (let i = 0; i < dist_to_others.length; i++) {
    const d = dist_to_others[i];
    if (d < d_max && i !== fish.id) {
      const dist_weight = 1 - (d - d_min) / (d_max - d_min);
      fish.social_nn.push(i);
      let direction_to_other = env.fishes[i].position.sub(p);
      if (d > 0) {
        direction_to_other = direction_to_other.div_scalar(d);
      }
      direction_to_other = force_rep.add(direction_to_other.mul_scalar(strength_att * dist_weight));
      if (d < d_min) {
        counter_rep += 1;
        force_rep = force_rep.add(direction_to_other.mul_scalar(-1));
      } else {
        force_att = force_att.add(direction_to_other);
        const v_meter_other = env.fishes[i].velocity.div_scalar(env.pixel_per_meter);
        const v_diff = v_meter_other.sub(v_meter);
        force_ali = force_ali.add(v_diff.mul_scalar(strength_align * dist_weight));
      }
    }
  }
  if (counter_rep > 0) {
    return force_rep.div_scalar(counter_rep);
  }
  return force_att.add(force_ali);
}

export class Fish {
  constructor(id, env, geojson) {
    this.env = env;
    this.id = id;
    this.timestamp = [];
    this.timestamp_states = [];
    this.positions = [];
    this.states = [];
    this.parameters = [];
    this.position = new Vector(0, 0, 0);
    this.velocity = new Vector(0, 0, 0);
    this.speed_ms = 0.0;
    this.phi = 0.0;
    this.theta = Math.PI / 2;
    this.turn_phi = 2;
    this.turn_theta = 5;
    this.state = 0;
    this.parameter_array = env.draw_state_parameter_array ? env.draw_state_parameter_array(env) : [];
    const initial_parameters = env.draw_state_parameter_from_array
      ? env.draw_state_parameter_from_array(this.parameter_array, this.state, env)
      : env.zero_vector ?? [0, 0, 0, 0, 0, 0, 0, 0];
    this.beta = initial_parameters[0];
    this.v0 = initial_parameters[1];
    this.D_phi = initial_parameters[2];
    this.D_theta = initial_parameters[3];
    this.D_v = initial_parameters[4];
    this.patch_strength = initial_parameters[5];
    this.strength_att = initial_parameters[6];
    this.strength_align = initial_parameters[7];
    this.social_nn = [];
    this.patch_sensing_length = 10;
    this.shore_time = 10;
    this.shore_strength = 5;
    this.depth_time = 10;
    this.depth_stength = 10;
    this.shore_avoidance_sign = 0;
    this.ground_avoidance_sign = 0;
    if (geojson) {
      this.createFish(id, geojson);
    }
    this.debug_sign = 1;
  }

  createFish(id, geojson) {
    if (geojson) {
      const polygonCoordinates = geojson.features[0].geometry.coordinates;
      let tempx, tempy;
      do {
        tempx = Math.random() * this.env.viewport.width;
        tempy = Math.random() * this.env.viewport.height;
      } while (!pointInsidePolygon([tempx, tempy], polygonCoordinates));
      this.position[0] = tempx;
      this.position[1] = tempy;
      this.position[2] = 0;
    }
    this.phi = Math.random() * Math.PI * 2;
    this.theta = Math.PI / 2 + Math.random() * Math.PI / 8;
    this.speed_ms = Math.max(0, randomNormal(this.v0, this.v0 / 4));
    this._computeVelocity_inPixelPerSec();
    this.positions.push(this.position);
    this.timestamp.push(0);
    this.timestamp_states.push(0);
    this.states.push(this.state);
    this.parameters.push([this.beta, this.v0, this.D_phi, this.D_theta, this.D_v, this.patch_strength, this.patch_sensing_length, this.strength_att, this.strength_align]);
  }

  updatePosition(polygonCoordinates, dt) {
    const force_social = socialForce(this, this.strength_att, this.strength_align, this.env);
    const force_patch_attraction = patchAttractionForce(this, this.patch_sensing_length, this.patch_strength, this.env);
    const force_shore_repulsion = shoreRepulsionForce(this, polygonCoordinates, this.shore_time, this.shore_strength);
    const force_surface_repulsion = surfaceRepulsionForce(this, this.depth_time, this.depth_stength);
    const force_ground_repulsion = groundRepulsionForce(this, this.depth_time, this.depth_stength, this.env);
    let force = new Vector(0, 0, 0)
      .add(force_social)
      .add(force_patch_attraction)
      .add(force_shore_repulsion)
      .add(force_surface_repulsion)
      .add(force_ground_repulsion);

    if (Number.isNaN(force[0])) {
      console.log("force_social:", force_social, "force_patch_attraction:", force_patch_attraction, "force_shore_repulsion:", force_shore_repulsion, "force_surface_repulsion:", force_surface_repulsion, "force_ground_repulsion:", force_ground_repulsion);
      console.log("Fish id:", this.id, "position", this.position, "velocity:", this.velocity, "force:", force, "speed_ms:", this.speed_ms, "phi:", this.phi, "theta:", this.theta);
      throw new Error("NaN force detected for fish id: " + this.id);
    }

    let speed_ms = this.speed_ms, phi = this.phi, theta = this.theta;
    const cos_phi = Math.cos(phi), sin_phi = Math.sin(phi);
    const cos_theta = Math.cos(theta), sin_theta = Math.sin(theta);
    const u_v = new Vector(cos_phi * sin_theta, sin_phi * sin_theta, cos_theta);
    const u_phi = new Vector(-sin_phi * sin_theta, cos_phi * sin_theta, 0);
    const u_theta = new Vector(cos_phi * cos_theta, sin_phi * cos_theta, -sin_theta);

    const force_v = force.dot(u_v);
    speed_ms += (this.beta * (this.v0 - speed_ms) + force_v) * dt;
    speed_ms += Math.sqrt(this.D_v * dt) * randomNormal(0, 1);
    speed_ms = Math.max(0, Math.min(speed_ms, 20 * this.v0));

    const force_phi = force.dot(u_phi);
    const rand_phi = Math.sqrt(this.D_phi * dt) * randomNormal(0, 1);
    phi += (force_phi * dt + rand_phi) / (speed_ms + this.turn_phi);

    const force_theta = force.dot(u_theta);
    const rand_theta = Math.sqrt(this.D_phi * dt) * randomNormal(0, 1);
    theta += (force_theta * dt + rand_theta) / (speed_ms + this.turn_theta);

    if (theta < 0) {
      phi += Math.PI;
      theta *= -1;
    }
    if (theta > Math.PI) {
      phi += Math.PI;
      theta = Math.PI - (theta - Math.PI);
    }
    this.phi = phi;
    this.theta = theta;
    this.speed_ms = speed_ms;

    const speed_ms_before = this.speed_ms;
    this._computeVelocity_inPixelPerSec();
    if (!this._preventGroundCollision(dt)) {
      if (!this._preventGroundCollision(dt)) {
        this._preventGroundCollision(dt);
      }
    }
    this.position = this.position.add(this.velocity.mul_scalar(dt));
    this.stayInPolygon(polygonCoordinates, dt);
    this.speed_ms = speed_ms_before;
    this._computeVelocity_inPixelPerSec();
  }

  stayInPolygon(polygonCoordinates, dt) {
    if (!pointInsidePolygon([this.position[0], this.position[1]], polygonCoordinates)) {
      this.phi += Math.PI;
      this.velocity = this.velocity.mul_scalar(-1);
      this.velocity[2] = -this.velocity[2];
      this.position = this.position.add(this.velocity.mul_scalar(3 * dt));
    }
  }

  _preventGroundCollision(dt) {
    const p = this.position, v = this.velocity, speed_ms = this.speed_ms;
    const p_next = p.add(v.mul_scalar(dt));
    const indices2D = depthIndicesWithConflicts(p, p_next, this.env.depth_map, this.env.depthResolution);
    if (indices2D.length === 0 || speed_ms === 0) {
      return true;
    }
    const indices = [...new Set(
      indices2D.map(index => this.env.depth_map_points_idxs[index[0]][index[1]])
    )];
    const triangles_as_indices = getTrianglesFromPointIndices(indices, this.env.triangles, this.env.depth_point_idx_to_triangle_starts);
    const [t_min, norm_min, tria_min] = closestCollisionWithTriangle(p, p_next, triangles_as_indices, this.env.depth_points, this.env.triangles_to_draw);
    if (t_min < 0 || 1 < t_min) {
      return true;
    }
    console.warn("Fish id:", this.id, "is colliding with ground");
    const projection_v = norm_min.dot(v);
    if (t_min === 1) {
      if (projection_v > 0) {
        this.velocity = v.mul_scalar(1.05);
        console.warn("Collision-type A: t=1 AND triangle-norm-projection > 0 --> accelerate (1.05) to above the ground");
        return true;
      } else {
        this.velocity = v.mul_scalar(0.95);
        return true;
      }
    } else if (projection_v > 0) {
      console.warn("Collision-type B: t<1 AND triangle-norm-projection > 0 --> next step will be above the ground");
      const pc = p.add(v.mul_scalar(dt * t_min));
      monitorBug(this, p, pc, tria_min, this.env.points_to_draw);
      return false;
    } else {
      console.warn("Collision-type C: t<1 AND triangle-norm-projection < 0 --> reflect at ground");
      const pc = p.add(v.mul_scalar(dt * t_min));
      monitorBug(this, p, pc, tria_min, this.env.points_to_draw);
      const almost = 0.9999;
      this.position = p.add(v.mul_scalar(almost * t_min * dt));
      this.velocity = this.velocity.add(norm_min.mul_scalar(2 * Math.abs(projection_v)));
      this.velocity = this.velocity.mul_scalar(1 - t_min * almost);
      this.phi = this.velocity.phi();
      this.theta = this.velocity.theta();
      this.speed_ms = this.velocity.norm() / this.env.pixel_per_meter;
      return false;
    }
  }

  _computeVelocity_inPixelPerSec() {
    this.velocity[0] = Math.cos(this.phi) * Math.sin(this.theta);
    this.velocity[1] = Math.sin(this.phi) * Math.sin(this.theta);
    this.velocity[2] = Math.cos(this.theta);
    this.velocity = this.velocity.mul_scalar(this.speed_ms * this.env.pixel_per_meter);
  }

  recordState(max_iter) {
    const timestamp_now = Math.abs(this.timestamp[this.timestamp.length - 1]) + 1;
    this.positions.push(this.position);
    this.timestamp.push(timestamp_now * this.debug_sign);
    this.debug_sign = 1;

    if (this.positions.length > max_iter) {
      this.positions.shift();
      this.timestamp.shift();
    }

    if (this.state !== this.states[this.states.length - 1]) {
      this.states.push(this.state);
      this.timestamp_states.push(timestamp_now);
      this.parameters.push([this.beta, this.v0, this.D_phi, this.D_theta, this.D_v, this.patch_strength, this.patch_sensing_length, this.strength_att, this.strength_align]);

      if ((timestamp_now - this.timestamp_states[0]) > max_iter) {
        this.states.shift();
        this.timestamp_states.shift();
        this.parameters.shift();
      }
    }
  }
}

export class FoodPatch {
  constructor(id, env, geojson) {
    this.env = env;
    this.id = id;
    this.position = new Vector(0, 0, 0);
    if (geojson) {
      this.createPatch(id, geojson);
    }
  }

  createPatch(id, geojson) {
    if (geojson) {
      const polygonCoordinates = geojson.features[0].geometry.coordinates;
      let tempx, tempy;
      do {
        tempx = Math.random() * this.env.viewport.width;
        tempy = Math.random() * this.env.viewport.height;
      } while (!pointInsidePolygon([tempx, tempy], polygonCoordinates));
      this.position[0] = tempx;
      this.position[1] = tempy;
      this.position[2] = 0;
    }
  }

  updatePatch() {
    if (this.env.geojson) {
      const polygonCoordinates = this.env.geojson.features[0].geometry.coordinates;
      const updateChance = 1 / this.env.food_patch_update_time;
      if (Math.random() < updateChance) {
        let tempx, tempy;
        do {
          tempx = Math.random() * this.env.viewport.width;
          tempy = Math.random() * this.env.viewport.height;
        } while (!pointInsidePolygon([tempx, tempy], polygonCoordinates));
        this.position[0] = tempx;
        this.position[1] = tempy;
        this.position[2] = 0;
      }
    }
  }
}

export function monitorBug(fish, p1, p2, triangle, debug_points) {
  if (!debug_points) return;
  debug_points.push([p1[0], p1[1], p1[2], 1]);
  debug_points.push([p2[0], p2[1], p2[2], -1]);
  triangle.forEach(tria_point => {
    debug_points.push([tria_point[0], tria_point[1], tria_point[2], -1]);
  });
  fish.debug_sign = -1;
}

export function calculateDistance(x1, y1, x2, y2) {
  return Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2));
}

export { crossProduct, dotProduct };
