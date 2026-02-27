export function calculateBbox(geojson) {
  let bbox = [Infinity, Infinity, -Infinity, -Infinity];
  if (geojson?.type === "FeatureCollection") {
    geojson.features.forEach(feature => {
      const coords = feature.geometry.coordinates;
      if (feature.geometry.type === "Point") {
        bbox[0] = Math.min(bbox[0], coords[0]);
        bbox[1] = Math.min(bbox[1], coords[1]);
        bbox[2] = Math.max(bbox[2], coords[0]);
        bbox[3] = Math.max(bbox[3], coords[1]);
      } else if (
        feature.geometry.type === "LineString" ||
        feature.geometry.type === "Polygon" ||
        feature.geometry.type === "MultiPoint"
      ) {
        coords.forEach(coord => {
          bbox[0] = Math.min(bbox[0], coord[0]);
          bbox[1] = Math.min(bbox[1], coord[1]);
          bbox[2] = Math.max(bbox[2], coord[0]);
          bbox[3] = Math.max(bbox[3], coord[1]);
        });
      }
    });
  }
  return bbox;
}

export function calculateMinMax2D(polygon) {
  let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
  polygon.forEach(coord => {
    minX = Math.min(minX, coord[0]);
    minY = Math.min(minY, coord[1]);
    maxX = Math.max(maxX, coord[0]);
    maxY = Math.max(maxY, coord[1]);
  });
  return [minX, minY, maxX, maxY];
}
