function randomNormal(mean, sd) {
  const u = 1 - Math.random();
  const v = Math.random();
  const num = Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
  return num * sd + mean;
}

export function generateTransitionProbs(dimension) {
  const transition_probs = new Array(dimension).fill(0).map(() => new Array(dimension).fill(0));
  const prob = 1 / (dimension - 1);

  for (let i = 0; i < dimension; i++) {
    for (let j = 0; j < dimension; j++) {
      if (i !== j) {
        transition_probs[i][j] = prob;
      }
    }
  }

  return transition_probs;
}

export function calculateTransitionMatrix(transition_probs, sojourn_times, dt_output) {
  const prob_staying = [];
  for (const state in transition_probs) {
    prob_staying[state] = 1 - (1 / (sojourn_times[state] * 60 / dt_output));
  }
  const transition_matrix = [];
  for (const state in transition_probs) {
    const probs = transition_probs[state];
    const row = [];
    for (const next_state in probs) {
      if (next_state == state) {
        row.push(prob_staying[state]);
      } else {
        row.push((1 - prob_staying[state]) * probs[next_state]);
      }
    }
    transition_matrix.push(row);
  }
  return transition_matrix;
}

export function drawStateParameters(state, env) {
  if (state in env.state_means_matrix) {
    let parameterVector = env.state_means_matrix[state].map(mean => {
      const result = randomNormal(mean, env.param.globalStateStdev);
      return result;
    });
    parameterVector = parameterVector.map(value => Math.max(0, Math.min(1, value)));
    parameterVector = parameterVector.map((value, idx) => {
      const scaled = value * (env.upper_limits[idx] - env.lower_limits[idx]) + env.lower_limits[idx];
      return scaled;
    });
    return parameterVector;
  }
  console.log("state not in dictionary");
  return false;
}

export function drawStateParameterArray(states_present, env) {
  const state_param_array = [];
  for (const state in states_present) {
    state_param_array[state] = [];
    for (let i = 0; i < env.substates[state]; i++) {
      const state_param_vector = drawStateParameters(state, env);
      state_param_array[state].push(state_param_vector);
    }
  }
  return state_param_array;
}

export function drawStateParameterFromArray(state_param_array, state, env) {
  if (state_param_array.length === 0) {
    return env.zero_vector;
  }
  const index = Math.floor(Math.random() * env.substates[state]);
  return state_param_array[state][index];
}

export function updateStateParameters(fish, new_state, env) {
  if (fish.parameter_array.length > 0) {
    const new_parameters = drawStateParameterFromArray(fish.parameter_array, new_state, env);
    fish.state = new_state;
    fish.beta = new_parameters[0];
    fish.v0 = new_parameters[1];
    fish.D_phi = new_parameters[2];
    fish.D_theta = new_parameters[3];
    fish.D_v = new_parameters[4];
    fish.patch_strength = new_parameters[5];
    fish.strength_att = new_parameters[6];
    fish.strength_align = new_parameters[7];
  } else {
    fish.v0 = 0;
    fish.D_phi = 0;
    fish.D_theta = 0;
    fish.D_v = 0;
    fish.patch_strength = 0;
    fish.strength_att = 0;
    fish.strength_align = 0;
  }
}

export function stateSwitch(fish, env) {
  if (env.states_present.length <= 1) return;
  const state = fish.state;
  let new_state = state;
  const probs = env.transition_matrix[state];
  const ran = Math.random();
  let cumulative_prob = 0;

  for (let i = 0; i < probs.length; i++) {
    cumulative_prob += probs[i];
    if (ran < cumulative_prob) {
      new_state = i;
      break;
    }
  }

  if (new_state !== fish.state) {
    updateStateParameters(fish, new_state, env);
  }
}
