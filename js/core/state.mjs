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

export function stateSwitch(fish, transition_matrix, states_present) {
  if (states_present.length <= 1) return;

  const state = fish.state;
  let new_state = state;
  const probs = transition_matrix[state];
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
    fish.setState(new_state);
  }
}
