import type { Camera } from '../Camera';
import { Ray, Matrix4, Vector3 } from "../math";
import { Sphere, Cylinder, Triangle } from "../shapes";

const descSort = (a: { distance: number; }, b: { distance: number; }) => {
  return a.distance - b.distance;
};



const viewProjectionMatrix = new Matrix4();

export class Raycaster {
  ray: Ray;
  near: number;
  far: number;
  precision = 0.0001;
  linePrecision = 0.2;
  constructor(origin: Vector3 | undefined, direction: Vector3 | undefined, far?: number, near?: number) {
    this.ray = new Ray(origin, direction);

    if (this.ray.direction.lengthSq() > 0) this.ray.direction.normalize();

    this.near = near || 0;
    this.far = far || Infinity;
  }

  set(origin: Vector3, direction: Vector3): void {
    this.ray.set(origin, direction);
  }

  setFromCamera(coords: { x: any; y: any; z: any; }, camera: Camera): void {
    if (!camera.ortho) {
      this.ray.origin.setFromMatrixPosition(camera.matrixWorld);
      this.ray.direction.set(coords.x, coords.y, coords.z);

      camera.projectionMatrixInverse.getInverse(camera.projectionMatrix);
      viewProjectionMatrix.multiplyMatrices(
        camera.matrixWorld,
        camera.projectionMatrixInverse
      );
      this.ray.direction.applyProjection(viewProjectionMatrix);
      this.ray.direction.sub(this.ray.origin).normalize();
    } else {
      this.ray.origin
        .set(
          coords.x,
          coords.y,
          (camera.near + camera.far) / (camera.near - camera.far)
        )
        .unproject(camera);
      this.ray.direction.set(0, 0, -1).transformDirection(camera.matrixWorld);
    }
  }

  intersectObjects(group: any, objects: string | any[]) {
    var intersects: any[] = [];
  
    for (var i = 0, l = objects.length; i < l; i++)
      intersectObject(group, objects[i], this, intersects);
  
    intersects.sort(descSort);
  
    return intersects;
  }
}

// [-1, 1]
const clamp = (x: number): number => {
  return Math.min(Math.max(x, -1), 1);
};

var sphere = new Sphere();
var cylinder = new Cylinder();
var triangle = new Triangle();
var w_0 = new Vector3(); // for cylinders, cylinder.c1 - ray.origin
var v1 = new Vector3(); // all purpose local vector
var v2 = new Vector3();
var v3 = new Vector3();
var matrixPosition = new Vector3();

//object is a Sphere or (Bounding) Box
export function intersectObject(group: { matrixWorld: Matrix4; }, clickable: { intersectionShape: any; boundingSphere: Sphere | undefined; }, raycaster: Raycaster, intersects: any[]) {
  matrixPosition.getPositionFromMatrix(group.matrixWorld);

  if (clickable.intersectionShape === undefined) return intersects;
  var intersectionShape = clickable.intersectionShape;
  var precision = raycaster.linePrecision;
  precision *= group.matrixWorld.getMaxScaleOnAxis();
  var precisionSq = precision * precision;

  //Check for intersection with clickable's bounding sphere, if it exists
  if (
    clickable.boundingSphere !== undefined &&
    clickable.boundingSphere instanceof Sphere
  ) {
    sphere.copy(clickable.boundingSphere);
    sphere.applyMatrix4(group.matrixWorld);
    if (!raycaster.ray.isIntersectionSphere(sphere)) {
      return intersects;
    }
  }

  //Iterate through intersection objects
  var i,
    il,
    norm,
    normProj,
    cylProj,
    rayProj,
    distance,
    closestDistSq,
    denom,
    discriminant,
    s,
    t,
    s_c,
    t_c;
  //triangle faces
  for (i = 0, il = intersectionShape.triangle.length; i < il; i++) {
    if (intersectionShape.triangle[i] instanceof Triangle) {
      triangle.copy(intersectionShape.triangle[i]);
      triangle.applyMatrix4(group.matrixWorld);

      norm = triangle.getNormal();

      normProj = raycaster.ray.direction.dot(norm);

      //face culling
      if (normProj >= 0) continue;

      w_0.subVectors(triangle.a, raycaster.ray.origin);

      distance = norm.dot(w_0) / normProj;

      if (distance < 0) continue;

      //intersects with plane, check if P inside triangle
      v1.copy(raycaster.ray.direction)
        .multiplyScalar(distance)
        .add(raycaster.ray.origin);
      v1.sub(triangle.a); // from pt a to intersection point P

      v2.copy(triangle.b).sub(triangle.a); // from pt a to b
      v3.copy(triangle.c).sub(triangle.a); // from pt a to c
      var b_dot_c = v2.dot(v3);
      var b_sq = v2.lengthSq();
      var c_sq = v3.lengthSq();

      // P = A + s(v2) + t(v3), inside trianle if 0 <= s, t <=1  and (s + t) <=0

      t =
        (b_sq * v1.dot(v3) - b_dot_c * v1.dot(v2)) /
        (b_sq * c_sq - b_dot_c * b_dot_c);

      if (t < 0 || t > 1) continue;

      s = (v1.dot(v2) - t * b_dot_c) / b_sq;

      if (s < 0 || s > 1 || s + t > 1) continue;
      else {
        intersects.push({ clickable: clickable, distance: distance });
      }
    }
  }
  //cylinders
  for (i = 0, il = intersectionShape.cylinder.length; i < il; i++) {
    if (intersectionShape.cylinder[i] instanceof Cylinder) {
      cylinder.copy(intersectionShape.cylinder[i]);
      cylinder.applyMatrix4(group.matrixWorld);

      w_0.subVectors(cylinder.c1, raycaster.ray.origin);

      cylProj = w_0.dot(cylinder.direction); // Dela
      rayProj = w_0.dot(raycaster.ray.direction); // Epsilon

      normProj = clamp(raycaster.ray.direction.dot(cylinder.direction)); // Beta

      denom = 1 - normProj * normProj;

      if (denom === 0.0) continue;

      s_c = (normProj * rayProj - cylProj) / denom;
      t_c = (rayProj - normProj * cylProj) / denom;

      v1.copy(cylinder.direction).multiplyScalar(s_c).add(cylinder.c1); // Q_c
      v2.copy(raycaster.ray.direction)
        .multiplyScalar(t_c)
        .add(raycaster.ray.origin); // P_c

      closestDistSq = v3.subVectors(v1, v2).lengthSq();
      var radiusSq = cylinder.radius * cylinder.radius;

      //Smoothing?
      //if (closestDistSq > radiusSq) radiusSq += precisionSq;

      // closest distance between ray and cylinder axis not greater than cylinder radius;
      // might intersect this cylinder between atom and bond midpoint
      if (closestDistSq <= radiusSq) {
        //Find points where ray intersects sides of cylinder
        discriminant =
          (normProj * cylProj - rayProj) * (normProj * cylProj - rayProj) -
          denom * (w_0.lengthSq() - cylProj * cylProj - radiusSq);

        // ray tangent to cylinder?
        if (discriminant <= 0) t = distance = Math.sqrt(closestDistSq);
        else
          t = distance =
            (rayProj - normProj * cylProj - Math.sqrt(discriminant)) / denom;

        //find closest intersection point; make sure it's between atom's position and cylinder midpoint

        s = normProj * t - cylProj;

        //does not intersect cylinder between atom and midpoint,
        // or intersects cylinder behind camera
        if (s < 0 || s * s > cylinder.lengthSq() || t < 0) continue;
        else intersects.push({ clickable: clickable, distance: distance });
      }
    }
  }
  //lines
  for (i = 0, il = intersectionShape.line.length; i < il; i += 2) {
    v1.copy(intersectionShape.line[i]);
    v1.applyMatrix4(group.matrixWorld);
    v2.copy(intersectionShape.line[i + 1]);
    v2.applyMatrix4(group.matrixWorld);

    v3.subVectors(v2, v1);
    var bondLengthSq = v3.lengthSq();
    v3.normalize();

    w_0.subVectors(v1, raycaster.ray.origin);

    var lineProj = w_0.dot(v3);
    rayProj = w_0.dot(raycaster.ray.direction);

    normProj = clamp(raycaster.ray.direction.dot(v3));

    denom = 1 - normProj * normProj;

    if (denom === 0.0) continue;

    s_c = (normProj * rayProj - lineProj) / denom;
    t_c = (rayProj - normProj * lineProj) / denom;

    v1.add(v3.multiplyScalar(s_c)); // Q_c
    v2.copy(raycaster.ray.direction)
      .multiplyScalar(t_c)
      .add(raycaster.ray.origin); // P_c

    closestDistSq = v3.subVectors(v2, v1).lengthSq();

    if (closestDistSq < precisionSq && s_c * s_c < bondLengthSq)
      intersects.push({ clickable: clickable, distance: t_c });
  }
  for (i = 0, il = intersectionShape.sphere.length; i < il; i++) {
    //sphere
    if (intersectionShape.sphere[i] instanceof Sphere) {
      sphere.copy(intersectionShape.sphere[i]);
      sphere.applyMatrix4(group.matrixWorld);

      if (raycaster.ray.isIntersectionSphere(sphere)) {
        v1.subVectors(sphere.center, raycaster.ray.origin);

        //distance from ray origin to point on the ray normal to sphere's center
        //must be less than sphere's radius (since ray intersects sphere)
        var distanceToCenter = v1.dot(raycaster.ray.direction);

        discriminant =
          distanceToCenter * distanceToCenter -
          (v1.lengthSq() - sphere.radius * sphere.radius);

        //Don't select if sphere center behind camera
        if (distanceToCenter < 0) return intersects;

        //ray tangent to sphere?
        if (discriminant <= 0) distance = distanceToCenter;
        //This is reversed if sphere is closer than ray origin.  Do we have
        //to worry about handling that case?
        else distance = distanceToCenter - Math.sqrt(discriminant);

        intersects.push({ clickable: clickable, distance: distance });
      }
    }
  }
  return intersects;
};
