
function calculateDistanceToPolygon(geojson) {
    // Extract coordinates of the polygon boundary
    const polygonCoordinates = geojson.features[0].geometry.coordinates; // Assuming it's the first feature and a Polygon
    
    // Calculate bounding box
    const bbox = calculateBbox(geojson);
    
    // Create an empty matrix to store distances
    const matrix = [];
    
    // Explode polygon vertices
    const vertices = turf.explode(geojson);
    
    // Loop through each point inside the bounding box of the polygon
    for (let y = 0; y <= canvas.height; y+=depthResolution) {
        const row = [];
        for (let x = 0; x <= canvas.width; x+=depthResolution) {
            // Convert pixel coordinates to coordinates in the map's CRS
            const lng = x;
            const lat = y;

        if(pointInsidePolygon([lng, lat], polygonCoordinates)){
            //polygonCoordinates must be flattened for the pointToPolygonDistance function **note!
            const distance = pointToPolygonDistance([lng, lat], polygonCoordinates.flat());
            row.push(distance);
        }else{
            row.push(NaN);   
        }
    }
    
    matrix.push(row);
}

return matrix;
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
        return !isNaN(value) ? value + noise1 + noise2 + noise3 : NaN;
    })
);

// Rescale combined matrix to keep depth between 0 and maxDepth
const maxCombinedValue = Math.max(...combinedMatrix.flat().filter(value => !isNaN(value)));
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

/**
 * Performs Delaunay triangulation on a set of 3D ground points and returns a list of triangles.
 *
 * @param {Array<Array<number>>} points - An array of 3D points, where each point is an array of the form [x, y, z].
 * @returns {Array<Array<Array<number>>>} An array of triangles, each represented as an array of three [x, y, z] points.
 *
 * @note The input `points` array must contain at least three points, each formatted as [x, y, z].
 */
function triangulate_ground(points) {
    // Extract 2D coordinates for triangulation
    const coords2D = points.map(p => [p[0], p[1]]);

    // Create the Delaunay triangulation
    const delaunay = d3.Delaunay.from(coords2D);

    // Get the triangles (indices into coords2D/points)
    const triangles = delaunay.triangles; // flat array: [i0, i1, i2, i3, i4, i5, ...]
    // Each consecutive triple is a triangle: [i0, i1, i2] is the first triangle, etc.

    // To get the actual 3D triangles:
    const triangleList = [];
    for (let i = 0; i < triangles.length; i += 3) {
        triangleList.push([
            points[triangles[i]],
            points[triangles[i + 1]],
            points[triangles[i + 2]]
        ]);
    }
    // triangleList is now an array of triangles, each triangle is an array of 3 [x, y, z] points
    return triangleList;
}