export function buildDepthMapFromData(geojson, options) {
  const { worldWidth, worldHeight, depthResolution } = options;
  void worldWidth;
  void worldHeight;
  void depthResolution;
  void geojson;

  // TODO: port sample_depth_map + triangulation pipeline.
  return {
    depth_map: [],
    depth_points: [],
    triangles: [],
  };
}

export function buildDepthMapFromShoreDistance(geojson, options) {
  const { worldWidth, worldHeight, depthResolution } = options;
  void worldWidth;
  void worldHeight;
  void depthResolution;
  void geojson;

  // TODO: port calculateDistanceToPolygon + mapDepth pipeline.
  return {
    depth_map: [],
    depth_points: [],
    triangles: [],
  };
}
