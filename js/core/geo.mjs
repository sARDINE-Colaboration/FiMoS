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

export function normalizeToLineString(feature) {
  if (feature && feature.geometry && feature.geometry.type !== "LineString") {
    const coordinates = feature.geometry.coordinates[0];
    feature.geometry = {
      type: "LineString",
      coordinates: coordinates,
    };
  }
  return feature;
}

export function rescaleCoordinates2D(coordinates, viewport, limitsOverride) {
  const limits = limitsOverride ?? (() => {
    const [minX, minY, maxX, maxY] = calculateMinMax2D(coordinates);
    return { minX, minY, maxX, maxY };
  })();
  const minX = limits.minX;
  const minY = limits.minY;
  const maxX = limits.maxX;
  const maxY = limits.maxY;
  let pixel_per_meter = viewport.width / (maxX - minX);
  let x_as_max_extent = true;
  let rescale_offset = (viewport.height - (maxY - minY) * pixel_per_meter) / 2;
  let rescaled_coords;

  if ((maxY - minY) * pixel_per_meter > viewport.height) {
    x_as_max_extent = false;
    pixel_per_meter = viewport.height / (maxY - minY);
    rescale_offset = (viewport.width - (maxX - minX) * pixel_per_meter) / 2;

    rescaled_coords = coordinates.map(coord => [
      (coord[0] - minX) * pixel_per_meter + rescale_offset,
      viewport.height - (coord[1] - minY) * pixel_per_meter,
    ]);
  } else {
    rescaled_coords = coordinates.map(coord => [
      (coord[0] - minX) * pixel_per_meter,
      viewport.height - ((coord[1] - minY) * pixel_per_meter + rescale_offset),
    ]);
  }

  const scalebar_length = Math.round(((maxX - minX) * 0.15) / 10) * 10;

  return {
    coords: rescaled_coords,
    pixel_per_meter,
    x_as_max_extent,
    rescale_offset,
    scalebar_length,
    limits: { minX, minY, maxX, maxY },
  };
}

export function rescaleCoordinates3D(coordinates, viewport, rescale, limits) {
  const minX = limits.minX;
  const minY = limits.minY;
  const maxX = limits.maxX;
  const maxY = limits.maxY;
  const pixel_per_meter = rescale.pixel_per_meter;
  const rescale_offset = rescale.rescale_offset;

  let rescaled_coords;
  if (rescale.x_as_max_extent) {
    rescaled_coords = coordinates.map(coord => [
      (coord[0] - minX) * pixel_per_meter,
      viewport.height - ((coord[1] - minY) * pixel_per_meter + rescale_offset),
      coord[2] * pixel_per_meter,
    ]);
  } else {
    rescaled_coords = coordinates.map(coord => [
      (coord[0] - minX) * pixel_per_meter + rescale_offset,
      viewport.height - (coord[1] - minY) * pixel_per_meter,
      coord[2] * pixel_per_meter,
    ]);
  }

  return rescaled_coords;
}

export function rescaleToGeoJSON(x, y, z, viewport, rescale, limits) {
  let scaledX;
  let scaledY;
  if (rescale.x_as_max_extent) {
    scaledX = x / rescale.pixel_per_meter + limits.minX;
    scaledY = (viewport.height - y + rescale.rescale_offset) / rescale.pixel_per_meter + limits.minY;
  } else {
    scaledX = (x - rescale.rescale_offset) / rescale.pixel_per_meter + limits.minX;
    scaledY = (viewport.height - y) / rescale.pixel_per_meter + limits.minY;
  }
  const scaledZ = z / rescale.pixel_per_meter;
  return [scaledX, scaledY, scaledZ];
}
