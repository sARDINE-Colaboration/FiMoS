
function calculateDistanceToPolygon(geojson) {
    // Extract coordinates of the polygon boundary
    const polygonCoordinates = geojson.features[0].geometry.coordinates; // Assuming it's the first feature and a Polygon
    
    // Create an empty matrix to store distances
    const matrix = [];
    
    // Loop through each point inside the bounding box of the polygon
    for (let y = 0; y <= canvas.height; y+=depthResolution) {
        const row = [];
        const row_indices = [];
        for (let x = 0; x <= canvas.width; x+=depthResolution) {
            if(pointInsidePolygon([x, y], polygonCoordinates)){
                //polygonCoordinates must be flattened for the pointToPolygonDistance function **note!
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

// this function takes a matrix of depth values and returns a matrix of indices and a list of points
function getMatrixIndicesAndPoints(matrix) {
    // Create an empty matrix to store distances
    const matrix_indices = [];  
    const depth_points = [];
    var valid_value_counter = 0;
    
    // Loop through the whole matrix
    for (let y = 0; y < matrix.length; y++) {
        const row_indices = [];
        for (let x = 0; x < matrix[y].length; x++) {
            const depth = matrix[y][x];
            if(!isNaN(depth)){
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

function depthIndicesWithConflicts(xyz, xyz_next, depth_map){
    // get the indices of the points in the depth map with resolution depthResolution
    let xIndex = Math.round(xyz[0] / depthResolution);
    let yIndex = Math.round(xyz[1] / depthResolution);
    const N_x = depth_map[0].length;
    const N_y = depth_map.length;

    const xIndex_next = Math.round(xyz_next[0] / depthResolution);
    const yIndex_next = Math.round(xyz_next[1] / depthResolution);
    // now take the min and max of the indices to get the submatrix
    // ensure to take the next neighbors in as well, so we substract/add 1 to the start/end indices
    const xStart = Math.max(Math.min(xIndex, xIndex_next) - 1, 0);
    const xEnd = Math.min(Math.max(xIndex, xIndex_next) + 1, N_x - 1);
    const yStart = Math.max(Math.min(yIndex, yIndex_next) -1, 0);
    const yEnd = Math.min(Math.max(yIndex, yIndex_next) + 1, N_y - 1);
    // the indices of the next point + 1 (triangle edge) and subtract 1 from the current point
    // Extract the submatrix using slice (like numpy's depth_map[yStart:yEnd+1, xStart:xEnd+1])
    const subMatrix = depth_map.slice(yStart, yEnd + 1).map(row => row.slice(xStart, xEnd + 1));

    // now check if the current and next depth are 
    const conflictIndices = [];
    // return the indices of the submatrix whose depths are larger than the depth reference
    for (let i = 0; i < subMatrix.length; i++) {
        for (let j = 0; j < subMatrix[i].length; j++) {
            let idx_y = yStart + i;
            let idx_x = xStart + j;
            if (isNaN(subMatrix[i][j])) continue; // skip NaN values
            conflictIndices.push([idx_y, idx_x]);
        }
    }
    // now add those indices with a neighbor NaN
    return conflictIndices;
}

/**
* this function samples the depth for a grid with resolution depthResolution from the 
* second feature, which is of type == MultiPoint Z
* for each grid point it takes the depth from the closest point from the depth feature
*/
function sample_depth_map(geojson) {
    // Extract coordinates of the polygon boundary
    const polygonCoordinates = geojson.features[0].geometry.coordinates; // Assuming it's the first feature and a Polygon
    // Extract depth coordinates
    const depthCoordinates = geojson.features[1].geometry.coordinates; // Assuming it's the second feature and a MultiPopint Z 
    // ensure the depth is negative
    const limitedDepthCoordinates = depthCoordinates.slice(0, 1000); // Use only the first 1000 entries
    const min_depth = Math.min(...limitedDepthCoordinates.map(coord => coord[2])); // Get the minimum depth from the limited coordinates
    var sign_factor = 1;
    if (min_depth >= 0) {
        sign_factor = -1; // If the minimum depth is positive, we assume the depth values are negative
    }
    
    // Create an empty matrix to store distances
    const matrix = [];
    
    // Loop through each point inside the bounding box of the polygon
    for (let y = 0; y <= canvas.height; y+=depthResolution) {
        const row = [];
        for (let x = 0; x <= canvas.width; x+=depthResolution) {
            if(pointInsidePolygon([x, y], polygonCoordinates)){
                //polygonCoordinates must be flattened for the pointToPolygonDistance function **note!
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
    let depth = 0; // Initialize depth to 0, will be updated if a point is found
    multiPoint.forEach(mp => {
        const distance = Math.sqrt(
            Math.pow(point[0] - mp[0], 2) + Math.pow(point[1] - mp[1], 2)
        );
        if (distance < minDistanceXY) {
            minDistanceXY = distance;
            depth = mp[2]; // Update depth with the z-coordinate of the closest point
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
    if (lenSq !== 0) //in case of 0 length line
        param = dot / lenSq;
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

function mapDepth(matrix) {  
    // ATTENTION: THIS IS MESSY. IF maxDepth < 0, the maxCombinedValue are actually the min
    // TODO ..... write it better (assume only min values from the start)
    // Rescale depth matrix
    const maxDepthValue = Math.max(...matrix.flat().filter(value => !isNaN(value)));
    const scaledDepthMatrix = matrix.map(row => 
        row.map(value => 
            !isNaN(value) ? (value / maxDepthValue) * maxDepth : NaN
        )
    );

    // Generate noise with different smoothness parameters
    const noiseMatrix1 = generateNoise(matrix[0].length, matrix.length, 0.024, 10); // Adjust parameters as needed
    const noiseMatrix2 = generateNoise(matrix[0].length, matrix.length, 0.03, 7); // Adjust parameters as needed
    const noiseMatrix3 = generateNoise(matrix[0].length, matrix.length, 0.05, 5); // Adjust parameters as needed

    // Add noise to the scaled depth matrix
    const combinedMatrix = scaledDepthMatrix.map((row, y) => 
        row.map((value, x) => {
            const noise1 = noiseMatrix1[y][x];
            const noise2 = noiseMatrix2[y][x];
            const noise3 = noiseMatrix3[y][x];
            return !isNaN(value) ? value - noise1 - noise2 - noise3 : NaN;
        })
    );

    // Rescale combined matrix to keep depth between 0 and maxDepth
    const maxCombinedValue = Math.min(...combinedMatrix.flat().filter(value => !isNaN(value)));
    const rescaledMatrix = combinedMatrix.map(row => 
        row.map(value => 
            !isNaN(value) ? (value / maxCombinedValue) * maxDepth : NaN
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


PerlinNoise = new function() {

    this.noise = function(x, y, z) {

       var p = new Array(512)
       var permutation = [ 151,160,137,91,90,15,
       131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
       190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
       88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
       77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
       102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
       135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
       5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
       223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
       129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
       251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
       49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
       138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180
       ];
       for (var i=0; i < 256 ; i++) 
        p[256+i] = p[i] = permutation[i]; 

        var X = Math.floor(x) & 255,                  // FIND UNIT CUBE THAT
            Y = Math.floor(y) & 255,                  // CONTAINS POINT.
            Z = Math.floor(z) & 255;
        x -= Math.floor(x);                                // FIND RELATIVE X,Y,Z
        y -= Math.floor(y);                                // OF POINT IN CUBE.
        z -= Math.floor(z);
        var    u = fade(x),                                // COMPUTE FADE CURVES
               v = fade(y),                                // FOR EACH OF X,Y,Z.
               w = fade(z);
        var A = p[X  ]+Y, AA = p[A]+Z, AB = p[A+1]+Z,      // HASH COORDINATES OF
            B = p[X+1]+Y, BA = p[B]+Z, BB = p[B+1]+Z;      // THE 8 CUBE CORNERS,

        return scale(lerp(w, lerp(v, lerp(u, grad(p[AA  ], x  , y  , z   ),  // AND ADD
                                       grad(p[BA  ], x-1, y  , z   )), // BLENDED
                               lerp(u, grad(p[AB  ], x  , y-1, z   ),  // RESULTS
                                       grad(p[BB  ], x-1, y-1, z   ))),// FROM  8
                       lerp(v, lerp(u, grad(p[AA+1], x  , y  , z-1 ),  // CORNERS
                                       grad(p[BA+1], x-1, y  , z-1 )), // OF CUBE
                               lerp(u, grad(p[AB+1], x  , y-1, z-1 ),
                                       grad(p[BB+1], x-1, y-1, z-1 )))));
   }
   function fade(t) { return t * t * t * (t * (t * 6 - 15) + 10); }
   function lerp( t, a, b) { return a + t * (b - a); }
   function grad(hash, x, y, z) {
      var h = hash & 15;                      // CONVERT LO 4 BITS OF HASH CODE
      var u = h<8 ? x : y,                 // INTO 12 GRADIENT DIRECTIONS.
             v = h<4 ? y : h==12||h==14 ? x : z;
      return ((h&1) == 0 ? u : -u) + ((h&2) == 0 ? v : -v);
   } 
   function scale(n) { return (1 + n)/2; }
}

// Function to convert depth to color
function depthToColor(depth) {
    const temp = (depth / maxDepth)
    const blue = Math.max(0, Math.min(255, Math.floor(temp*255)));
    const red = 255 - blue;
    const green = 0;
    return [red, green, blue, 255]; // Return RGB values in an array
}

// Function to get depth map colors
function getDepthMapColors(depthInput) {
    const colors = depthInput.map(row => {
      return row.map(depth => {
        if (!isNaN(depth)) {
          const color = depthToColor(depth);
          return color;
        } else {
          return [0, 0, 0, 0]; // Transparent black color for NaN values
        }
      });
    });
    return colors;
}

function makeDepthMapImageData(colors) {
    // Rotate the colors matrix 90 degrees and flip on the vertical plane
    const rotatedColors = colors//.slice().reverse();

    // Get the 2D context of the canvas
    const ctx = canvas.getContext('2d');

    // Get canvas dimensions
    const width = canvas.width;
    const height = canvas.height;

    // Create ImageData object
    const imgData = ctx.createImageData(width, height);

    // Scale factor
    const scaleX = rotatedColors[0].length / width;
    const scaleY = rotatedColors.length / height;

    // Fill ImageData with the provided colors
    for (let y = 0; y < height; y++) {
      for (let x = 0; x < width; x++) {
        const color = rotatedColors[Math.floor(y * scaleY)][Math.floor(x * scaleX)];
        const index = (y * width + x) * 4;
        imgData.data[index] = color[0];     // Red channel
        imgData.data[index + 1] = color[1]; // Green channel
        imgData.data[index + 2] = color[2]; // Blue channel
        imgData.data[index + 3] = color[3]; // Alpha channel
      }
    }

    return imgData;
}

// function that returns the points of the shoreline but if the distance is larger than depthResolution, it interpolates the missing points
function pointsOfShorelineInterpolated(coordinates) {
    const shorelinePoints = [];
    for (let i = 0; i < coordinates.length - 1; i++) {
        const start = coordinates[i];
        const end = coordinates[i + 1];
        shorelinePoints.push([start[0], start[1], 0]); // Add the start point with z = 0
        
        // Calculate the distance between the two points
        const dx = end[0] - start[0];
        const dy = end[1] - start[1];
        const distance = Math.sqrt(dx * dx + dy * dy);
        
        // If the distance is greater than depthResolution, interpolate points
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
    shorelinePoints.push([lastPoint[0], lastPoint[1], 0]); // Add the last point
    return shorelinePoints;
}

// Performs Delaunay triangulation on a set of 3D ground points and returns a list of triangles.
function triangulate_ground(depth_points, shore_points){
    // interpolate shore points to ensure they are evenly spaced
    const shore_points_interpolated = pointsOfShorelineInterpolated(shore_points);
    // Combine depth points and shore points
    const combinedPoints = [...depth_points, ...shore_points_interpolated];

    // Extract 2D coordinates for triangulation
    const coords2D = combinedPoints.map(p => [p[0], p[1]]);

    // Create the Delaunay triangulation
    const delaunay = d3.Delaunay.from(coords2D);

    // Get the triangles (indices into coords2D/points)
    const triangles = delaunay.triangles; // flat array: [i0, i1, i2, i3, i4, i5, ...]
    // Each consecutive triple is a triangle: [i0, i1, i2] is the first triangle, etc.

    return [combinedPoints, triangles];
}

// Function to get a mapping of point index to triangles
function mapDepthPointIdxToTriangleStarts(depth_points_length, triangles){
    // create empty array of arrays, where each index corresponds to a point index
    const pointIdxToTriangleMap = Array.from({ length: depth_points_length }, () => []);
    // fill the map with triangle start indices for each point index
    for (let i = 0; i < triangles.length; i++) {
        let triangle_startIdx = Math.floor(i / 3) * 3; // Each triangle consists of 3 indices
        let pt_idx = triangles[i];
        if (pt_idx < depth_points_length) {
            pointIdxToTriangleMap[pt_idx].push(triangle_startIdx);
        }
    }
    // remove duplicates (should not be necessary)
    for (let i = 0; i < pointIdxToTriangleMap.length; i++) {
        pointIdxToTriangleMap[i] = [...new Set(pointIdxToTriangleMap[i])];
    }
    return pointIdxToTriangleMap;
}

function getTrianglesFromPointIndices(indices, triangles, indices_to_triangle_start) {
    // find location of indices in the triangles array
    const triangle_start_indices = [];
    for (let i = 0; i < indices.length; i++) {
        // Find all indices of idx[i] in triangles
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
    // remove duplicates
    const unique_triangle_start_indices = [...new Set(triangle_start_indices)];
    // Now extract the triangles from the unique indices
    const triangles_as_points_indices = unique_triangle_start_indices.map(index => {
        return [
            triangles[index],
            triangles[index + 1],
            triangles[index + 2]
        ];
    });
    // output is array of triangles, where each triangle is an array of point indices
    // e.g. [[0, 1, 2], [3, 4, 5], ...]
    return triangles_as_points_indices;
}

function getTriangleCoords(triangle, points) {
    // Extract the coordinates of the triangle vertices
    const p1 = points[triangle[0]];
    const p2 = points[triangle[1]];
    const p3 = points[triangle[2]];
    return [p1, p2, p3];
}

function crossProduct(a, b) {
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    ];
}

function dotProduct(v1, v2) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

function lineIntersectsTriangle(lineStart, lineEnd, triangle) {
    // Extract triangle vertices
    const [v0, v1, v2] = triangle;

    // Convert line segment to vector form (reversed direction)
    const lineDir_n = [lineStart[0] - lineEnd[0], lineStart[1] - lineEnd[1], lineStart[2] - lineEnd[2]];
    
    // Compute edges of the triangle
    const edge1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
    const edge2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
    // Compute the normal of the triangle
    const normal = crossProduct(edge1, edge2);

    // Check if the line is parallel to the triangle plane
    const lineDirN_dot_norm = dotProduct(lineDir_n, normal);
    if (Math.abs(lineDirN_dot_norm) < 1e-5) {
        return [1.5, normal]; // Line is parallel to the triangle plane
    }

    // Find intersection point with the triangle plane: following https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
    // computing t, u, v
    const ls_minus_v0 = [lineStart[0] - v0[0], lineStart[1] - v0[1], lineStart[2] - v0[2]];
    const t =  dotProduct(normal, ls_minus_v0) / lineDirN_dot_norm;
    if (t < 0 || t > 1) {
        return [1.5, normal]; // Intersection point is outside the segment
    }

    const u = dotProduct(crossProduct(edge2, lineDir_n), ls_minus_v0) / lineDirN_dot_norm;
    //if (u < 0 || u > 1) {
    if (u < -0.1 || u > 1.1) { // condition is relaxed, actually it is for u>1, but we allow a small margin
        return [1.5, normal]; // Intersection point is outside the segment
    }
    
    const v = dotProduct(crossProduct(lineDir_n, edge1), ls_minus_v0) / lineDirN_dot_norm;
    // if (u + v < 0 || u + v > 1) {
    if (u + v < -0.2 || u + v > 1.2) { // condition is relaxed, actually it is for u+v>1, but we allow a small margin
        return [1.5, normal]; // Intersection point is outside the segment
    }
    return [t, normal];
}