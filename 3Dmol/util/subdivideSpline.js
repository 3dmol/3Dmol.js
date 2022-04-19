// helper functions
// Catmull-Rom subdivision

import { Vector3 } from "../WebGL/math";

// eslint-disable-next-line camelcase
export default function subdivideSpline(_points, DIV) {
  // points as Vector3
  const ret = [];
  let points = _points;
  points = []; // Smoothing test
  points.push(_points[0]);

  let i;
  let lim;
  let size;
  let p0;
  let p1;
  let p2;
  let p3;
  let v0;
  let v1;

  for (i = 1, lim = _points.length - 1; i < lim; i++) {
    p1 = _points[i];
    p2 = _points[i + 1];
    if (p1.smoothen) {
      const np = new Vector3((p1.x + p2.x) / 2, (p1.y + p2.y) / 2, (p1.z + p2.z) / 2);
      np.atom = p1.atom;
      points.push(np);
    } else points.push(p1);
  }
  points.push(_points[_points.length - 1]);

  for (i = -1, size = points.length; i <= size - 3; i++) {
    p0 = points[i === -1 ? 0 : i];
    p1 = points[i + 1];
    p2 = points[i + 2];
    p3 = points[i === size - 3 ? size - 1 : i + 3];
    v0 = new Vector3().subVectors(p2, p0).multiplyScalar(0.5);
    v1 = new Vector3().subVectors(p3, p1).multiplyScalar(0.5);
    if (p2.skip) continue;

    for (let j = 0; j < DIV; j++) {
      const t = (1.0 / DIV) * j;
      const x =
        p1.x +
        t * v0.x +
        t * t * (-3 * p1.x + 3 * p2.x - 2 * v0.x - v1.x) +
        t * t * t * (2 * p1.x - 2 * p2.x + v0.x + v1.x);
      const y =
        p1.y +
        t * v0.y +
        t * t * (-3 * p1.y + 3 * p2.y - 2 * v0.y - v1.y) +
        t * t * t * (2 * p1.y - 2 * p2.y + v0.y + v1.y);
      const z =
        p1.z +
        t * v0.z +
        t * t * (-3 * p1.z + 3 * p2.z - 2 * v0.z - v1.z) +
        t * t * t * (2 * p1.z - 2 * p2.z + v0.z + v1.z);

      const pt = new Vector3(x, y, z);
      if (j < DIV / 2) {
        pt.atom = p1.atom;
      } else {
        pt.atom = p2.atom;
      }

      ret.push(pt);
    }
  }
  ret.push(points[points.length - 1]);
  return ret;
}
