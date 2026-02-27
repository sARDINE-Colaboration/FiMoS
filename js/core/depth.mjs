function pointInsidePolygon(point, polygon) {
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

export function calculateDistanceToPolygon(geojson, options) {
  const { worldWidth, worldHeight, depthResolution } = options;
  const polygonCoordinates = geojson.features[0].geometry.coordinates;
  const matrix = [];

  for (let y = 0; y <= worldHeight; y += depthResolution) {
    const row = [];
    for (let x = 0; x <= worldWidth; x += depthResolution) {
      if (pointInsidePolygon([x, y], polygonCoordinates)) {
        const distance = pointToPolygonDistance([x, y], polygonCoordinates.flat());
        row.push(distance);
      } else {
        row.push(NaN);
      }
    }
    matrix.push(row);
  }
  return matrix;
}

export function getMatrixIndicesAndPoints(matrix, depthResolution) {
  const matrix_indices = [];
  const depth_points = [];
  let valid_value_counter = 0;

  for (let y = 0; y < matrix.length; y++) {
    const row_indices = [];
    for (let x = 0; x < matrix[y].length; x++) {
      const depth = matrix[y][x];
      if (!Number.isNaN(depth)) {
        row_indices.push(valid_value_counter);
        valid_value_counter++;
        const point = [x * depthResolution, y * depthResolution, depth];
        depth_points.push(point);
      } else {
        row_indices.push(NaN);
      }
    }
    matrix_indices.push(row_indices);
  }
  return [matrix_indices, depth_points];
}

export function sampleDepthMap(geojson, options) {
  const { worldWidth, worldHeight, depthResolution } = options;
  const polygonCoordinates = geojson.features[0].geometry.coordinates;
  const depthCoordinates = geojson.features[1].geometry.coordinates;
  const limitedDepthCoordinates = depthCoordinates.slice(0, 1000);
  const min_depth = Math.min(...limitedDepthCoordinates.map(coord => coord[2]));
  let sign_factor = 1;
  if (min_depth >= 0) {
    sign_factor = -1;
  }

  const matrix = [];
  for (let y = 0; y <= worldHeight; y += depthResolution) {
    const row = [];
    for (let x = 0; x <= worldWidth; x += depthResolution) {
      if (pointInsidePolygon([x, y], polygonCoordinates)) {
        const depth = sign_factor * getDepthOfClosestMultipoint([x, y], depthCoordinates);
        row.push(depth);
      } else {
        row.push(NaN);
      }
    }
    matrix.push(row);
  }
  return matrix;
}

function getDepthOfClosestMultipoint(point, multiPoint) {
  let minDistanceXY = Infinity;
  let depth = 0;
  multiPoint.forEach(mp => {
    const distance = Math.sqrt(
      Math.pow(point[0] - mp[0], 2) + Math.pow(point[1] - mp[1], 2)
    );
    if (distance < minDistanceXY) {
      minDistanceXY = distance;
      depth = mp[2];
    }
    minDistanceXY = Math.min(minDistanceXY, distance);
  });
  return depth;
}

function pointToPolygonDistance(point, polygon) {
  let minDistance = Infinity;
  for (let i = 0; i < polygon.length; i += 2) {
    const x1 = polygon[i];
    const y1 = polygon[i + 1];
    const x2 = polygon[(i + 2) % polygon.length];
    const y2 = polygon[(i + 3) % polygon.length];
    const distance = pointToSegmentDistance(point, [x1, y1], [x2, y2]);
    minDistance = Math.min(minDistance, distance);
  }
  return minDistance;
}

function pointToSegmentDistance(point, p1, p2) {
  const [x, y] = point;
  const [x1, y1] = p1;
  const [x2, y2] = p2;
  const A = x - x1;
  const B = y - y1;
  const C = x2 - x1;
  const D = y2 - y1;
  const dot = A * C + B * D;
  const lenSq = C * C + D * D;
  let param = -1;
  if (lenSq !== 0) param = dot / lenSq;
  let xx, yy;
  if (param < 0) {
    xx = x1;
    yy = y1;
  } else if (param > 1) {
    xx = x2;
    yy = y2;
  } else {
    xx = x1 + param * C;
    yy = y1 + param * D;
  }
  const dx = x - xx;
  const dy = y - yy;
  return Math.sqrt(dx * dx + dy * dy);
}

export function mapDepth(matrix, maxDepth) {
  const maxDepthValue = Math.max(...matrix.flat().filter(value => !Number.isNaN(value)));
  const scaledDepthMatrix = matrix.map(row =>
    row.map(value =>
      !Number.isNaN(value) ? (value / maxDepthValue) * maxDepth : NaN
    )
  );

  const noiseMatrix1 = generateNoise(matrix[0].length, matrix.length, 0.024, 10);
  const noiseMatrix2 = generateNoise(matrix[0].length, matrix.length, 0.03, 7);
  const noiseMatrix3 = generateNoise(matrix[0].length, matrix.length, 0.05, 5);

  const combinedMatrix = scaledDepthMatrix.map((row, y) =>
    row.map((value, x) => {
      const noise1 = noiseMatrix1[y][x];
      const noise2 = noiseMatrix2[y][x];
      const noise3 = noiseMatrix3[y][x];
      return !Number.isNaN(value) ? value - noise1 - noise2 - noise3 : NaN;
    })
  );

  const maxCombinedValue = Math.min(...combinedMatrix.flat().filter(value => !Number.isNaN(value)));
  const rescaledMatrix = combinedMatrix.map(row =>
    row.map(value =>
      !Number.isNaN(value) ? (value / maxCombinedValue) * maxDepth : NaN
    )
  );

  return rescaledMatrix;
}

function generateNoise(width, height, smoothness, amplitude) {
  const noise = [];
  for (let y = 0; y < height; y++) {
    noise[y] = [];
    for (let x = 0; x < width; x++) {
      noise[y][x] = PerlinNoise.noise(x * smoothness, y * smoothness, 0) * amplitude;
    }
  }
  return noise;
}

const PerlinNoise = new function () {
  this.noise = function (x, y, z) {
    const p = new Array(512);
    const permutation = [
      151, 160, 137, 91, 90, 15,
      131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23,
      190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33,
      88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166,
      77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244,
      102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196,
      135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123,
      5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42,
      223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9,
      129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228,
      251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107,
      49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254,
      138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180
    ];

    for (let i = 0; i < 256; i++) p[256 + i] = p[i] = permutation[i];

    let X = Math.floor(x) & 255;
    let Y = Math.floor(y) & 255;
    let Z = Math.floor(z) & 255;
    x -= Math.floor(x);
    y -= Math.floor(y);
    z -= Math.floor(z);
    const u = fade(x);
    const v = fade(y);
    const w = fade(z);
    const A = p[X] + Y, AA = p[A] + Z, AB = p[A + 1] + Z;
    const B = p[X + 1] + Y, BA = p[B] + Z, BB = p[B + 1] + Z;

    return scale(lerp(w, lerp(v, lerp(u, grad(p[AA], x, y, z),
      grad(p[BA], x - 1, y, z)),
      lerp(u, grad(p[AB], x, y - 1, z),
        grad(p[BB], x - 1, y - 1, z))),
      lerp(v, lerp(u, grad(p[AA + 1], x, y, z - 1),
        grad(p[BA + 1], x - 1, y, z - 1)),
        lerp(u, grad(p[AB + 1], x, y - 1, z - 1),
          grad(p[BB + 1], x - 1, y - 1, z - 1)))));
  };

  function fade(t) { return t * t * t * (t * (t * 6 - 15) + 10); }
  function lerp(t, a, b) { return a + t * (b - a); }
  function grad(hash, x, y, z) {
    const h = hash & 15;
    const u = h < 8 ? x : y;
    const v = h < 4 ? y : h === 12 || h === 14 ? x : z;
    return ((h & 1) === 0 ? u : -u) + ((h & 2) === 0 ? v : -v);
  }
  function scale(n) { return (1 + n) / 2; }
}();

function depthToColor(depth, maxDepth) {
  const temp = (depth / maxDepth);
  const blue = Math.max(0, Math.min(255, Math.floor(temp * 255)));
  const red = 255 - blue;
  const green = 0;
  return [red, green, blue, 255];
}

export function getDepthMapColors(depthInput, maxDepth) {
  const colors = depthInput.map(row => {
    return row.map(depth => {
      if (!Number.isNaN(depth)) {
        const color = depthToColor(depth, maxDepth);
        return color;
      }
      return [0, 0, 0, 0];
    });
  });
  return colors;
}

function pointsOfShorelineInterpolated(coordinates, depthResolution) {
  const shorelinePoints = [];
  for (let i = 0; i < coordinates.length - 1; i++) {
    const start = coordinates[i];
    const end = coordinates[i + 1];
    shorelinePoints.push([start[0], start[1], 0]);

    const dx = end[0] - start[0];
    const dy = end[1] - start[1];
    const distance = Math.sqrt(dx * dx + dy * dy);

    if (distance > depthResolution) {
      const numPoints = Math.ceil(distance / depthResolution);
      for (let j = 1; j < numPoints; j++) {
        const interpolatedX = start[0] + (dx * j) / numPoints;
        const interpolatedY = start[1] + (dy * j) / numPoints;
        shorelinePoints.push([interpolatedX, interpolatedY, 0]);
      }
    }
  }
  const lastPoint = coordinates[coordinates.length - 1];
  shorelinePoints.push([lastPoint[0], lastPoint[1], 0]);
  return shorelinePoints;
}

export function triangulateGround(depth_points, shore_points, depthResolution, delaunay) {
  if (!delaunay) {
    return [depth_points, []];
  }
  const shore_points_interpolated = pointsOfShorelineInterpolated(shore_points, depthResolution);
  const combinedPoints = [...depth_points, ...shore_points_interpolated];
  const coords2D = combinedPoints.map(p => [p[0], p[1]]);
  const delaunayTriangulation = delaunay.from(coords2D);
  const triangles = delaunayTriangulation.triangles;
  return [combinedPoints, triangles];
}

export function mapDepthPointIdxToTriangleStarts(depth_points_length, triangles) {
  const pointIdxToTriangleMap = Array.from({ length: depth_points_length }, () => []);
  for (let i = 0; i < triangles.length; i++) {
    const triangle_startIdx = Math.floor(i / 3) * 3;
    const pt_idx = triangles[i];
    if (pt_idx < depth_points_length) {
      pointIdxToTriangleMap[pt_idx].push(triangle_startIdx);
    }
  }
  for (let i = 0; i < pointIdxToTriangleMap.length; i++) {
    pointIdxToTriangleMap[i] = [...new Set(pointIdxToTriangleMap[i])];
  }
  return pointIdxToTriangleMap;
}
