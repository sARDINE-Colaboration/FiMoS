export function buildTracksCsvFromRows(rows) {
  const lines = [];
  lines.push("Fish ID,X,Y,Z,Timestamp");
  rows.forEach(row => {
    lines.push(row.join(","));
  });
  return lines.join("\n");
}

export function buildStatesCsv(fishes) {
  const lines = [];
  lines.push("Fish ID,state,Timestamp,beta,v0,D_phi,D_theta,D_v,patch_strength,patch_dist,social_strength,social_align");
  fishes.forEach((fish, fishIndex) => {
    fish.states.forEach((state, idx) => {
      const timestamp = fish.timestamp_states[idx];
      const params = fish.parameters[idx] || [];
      lines.push([
        fishIndex + 1,
        state,
        timestamp,
        params[0],
        params[1],
        params[2],
        params[3],
        params[4],
        params[5],
        params[6],
        params[7],
        params[8],
      ].join(","));
    });
  });
  return lines.join("\n");
}

export function buildLakeTrianglesCsv(triangles, points) {
  const lines = [];
  lines.push("triangle_id,point_1_x,point_1_y,point_1_z,point_2_x,point_2_y,point_2_z,point_3_x,point_3_y,point_3_z");
  for (let i = 0; i < triangles.length; i += 3) {
    const triangle = triangles.slice(i, i + 3);
    const point1 = points[triangle[0]];
    const point2 = points[triangle[1]];
    const point3 = points[triangle[2]];
    lines.push(
      `${i / 3 + 1},${point1[0]},${point1[1]},${point1[2]},${point2[0]},${point2[1]},${point2[2]},${point3[0]},${point3[1]},${point3[2]}`
    );
  }
  return lines.join("\n");
}

export function buildDebugPointsCsvFromRows(rows) {
  const lines = [];
  lines.push("X,Y,Z,debugFlag");
  rows.forEach(row => {
    lines.push(row.join(","));
  });
  return lines.join("\n");
}
