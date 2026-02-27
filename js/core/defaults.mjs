export const parameterNames = ['beta', 'v0', 'D_phi', 'D_theta', 'D_v', 'patch_strength', 'strength_att', 'strength_align'];

export const zero_vector = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
export const resting_vector = [0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
export const foraging_vector = [0.1, 0.4, 0.9, 0.9, 0.8, 0.0, 0.0, 0.0];
export const active_vector = [0.8, 0.8, 0.1, 0.1, 0.2, 0.0, 0.0, 0.0];

export const resting_sojourn = [20];
export const foraging_sojourn = [20];
export const active_sojourn = [20];

export function createDefaultStateConfig() {
  const all_states = ['resting', 'foraging', 'active'];
  const all_state_means = [
    [...resting_vector],
    [...foraging_vector],
    [...active_vector],
  ];
  const all_sojourn_times = [
    [...resting_sojourn],
    [...foraging_sojourn],
    [...active_sojourn],
  ];

  return {
    parameterNames: [...parameterNames],
    zero_vector: [...zero_vector],
    all_states,
    all_state_means,
    all_sojourn_times,
    state_means_matrix: [[...resting_vector]],
    states_present: ['resting'],
    substates: [[2]],
    sojourn_times: [resting_sojourn],
    transition_probs: [
      [0, 0.5, 0.5],
      [0.5, 0, 0.5],
      [0.5, 0.5, 0],
    ],
  };
}
