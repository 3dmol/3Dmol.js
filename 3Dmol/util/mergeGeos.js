// @ts-nocheck
// Miscellaneous functions and classes - to be incorporated into $3Dmol proper
/**
 *
 * @param {Geometry} geometry
 * @param {Mesh} mesh
 * @returns {undefined}
 */
export default function mergeGeos(geometry, mesh) {
  const meshGeo = mesh.geometry;

  if (meshGeo === undefined) return;

  geometry.geometryGroups.push(meshGeo.geometryGroups[0]);
}
