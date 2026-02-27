export function depthIndicesWithConflicts(xyz, xyz_next, depth_map, depthResolution) {
  let xIndex = Math.round(xyz[0] / depthResolution);
  let yIndex = Math.round(xyz[1] / depthResolution);
  const N_x = depth_map[0].length;
  const N_y = depth_map.length;

  const xIndex_next = Math.round(xyz_next[0] / depthResolution);
  const yIndex_next = Math.round(xyz_next[1] / depthResolution);

  const xStart = Math.max(Math.min(xIndex, xIndex_next) - 1, 0);
  const xEnd = Math.min(Math.max(xIndex, xIndex_next) + 1, N_x - 1);
  const yStart = Math.max(Math.min(yIndex, yIndex_next) - 1, 0);
  const yEnd = Math.min(Math.max(yIndex, yIndex_next) + 1, N_y - 1);

  const subMatrix = depth_map.slice(yStart, yEnd + 1).map(row => row.slice(xStart, xEnd + 1));

  const conflictIndices = [];
  for (let i = 0; i < subMatrix.length; i++) {
    for (let j = 0; j < subMatrix[i].length; j++) {
      const idx_y = yStart + i;
      const idx_x = xStart + j;
      if (Number.isNaN(subMatrix[i][j])) continue;
      conflictIndices.push([idx_y, idx_x]);
    }
  }
  return conflictIndices;
}

export function getTrianglesFromPointIndices(indices, triangles, indices_to_triangle_start) {
  const triangle_start_indices = [];
  for (let i = 0; i < indices.length; i++) {
    try {
      triangle_start_indices.push(...indices_to_triangle_start[indices[i]]);
    } catch (e) {
      if (e instanceof TypeError) {
        console.log('ERROR, there is no entry in indices_to_triangle_start at: ', indices[i]);
        console.log('indices are: ', indices);
        console.log('indices_to_triangle_start.length: ', indices_to_triangle_start.length);
      } else {
        throw e;
      }
    }
  }
  const unique_triangle_start_indices = [...new Set(triangle_start_indices)];
  return unique_triangle_start_indices.map(index => ([
    triangles[index],
    triangles[index + 1],
    triangles[index + 2],
  ]));
}

export function getTriangleCoords(triangle, points) {
  const p1 = points[triangle[0]];
  const p2 = points[triangle[1]];
  const p3 = points[triangle[2]];
  return [p1, p2, p3];
}

export function crossProduct(a, b) {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ];
}

export function dotProduct(v1, v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

export function lineIntersectsTriangle_ppk(lineStart, lineEnd, triangle) {
  const [v0, v1, v2] = triangle;
  const lineDir_n = [lineStart[0] - lineEnd[0], lineStart[1] - lineEnd[1], lineStart[2] - lineEnd[2]];

  const edge1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
  const edge2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
  const normal = crossProduct(edge1, edge2);

  const lineDirN_dot_norm = dotProduct(lineDir_n, normal);
  if (Math.abs(lineDirN_dot_norm) < 1e-5) {
    return [1.5, normal];
  }

  const ls_minus_v0 = [lineStart[0] - v0[0], lineStart[1] - v0[1], lineStart[2] - v0[2]];
  const t = dotProduct(normal, ls_minus_v0) / lineDirN_dot_norm;
  if (t < 0 || t > 1) {
    return [1.5, normal];
  }

  const u = dotProduct(crossProduct(edge2, lineDir_n), ls_minus_v0) / lineDirN_dot_norm;
  if (u < 0 || u > 1) {
    return [1.5, normal];
  }

  const v = dotProduct(crossProduct(lineDir_n, edge1), ls_minus_v0) / lineDirN_dot_norm;
  if (u + v < 0 || u + v > 1) {
    return [1.5, normal];
  }
  return [t, normal];
}
