// @ts-check

import { Vector3 } from "./Vector3";


// used in recast, distanceToPoint
let v1 = new Vector3();

// TODO: Remove methods we don't need (intersectPlane ??)
/** @constructor */
export class Ray {
  constructor(origin, direction) {

    this.origin = (origin !== undefined) ? origin : new Vector3();

    this.direction = (direction !== undefined) ? direction
      : new Vector3();

  }

  set(origin, direction) {

    this.origin.copy(origin);
    this.direction.copy(direction);

    return this;

  }

  copy(ray) {

    this.origin.copy(ray.origin);
    this.direction.copy(ray.direction);

    return this;

  }

  at(t, optionalTarget) {

    var result = optionalTarget || new Vector3();

    return result.copy(this.direction).multiplyScalar(t).add(this.origin);

  }

  //uses module level allocation of Vector3 v1
  recast(t) {
    this.origin.copy(this.at(t, v1));

    return this;
  }

  closestPointToPoint(point, optionalTarget) {

    var result = optionalTarget || new Vector3();
    result.subVectors(point, this.origin);
    var directionDistance = result.dot(this.direction);

    // returns a point on this ray
    return result.copy(this.direction).multiplyScalar(directionDistance)
      .add(this.origin);

  }

  //uses module level allocation of Vector3 v1
  distanceToPoint(point) {
    var directionDistance = v1.subVectors(point, this.origin).dot(
      this.direction);
    v1.copy(this.direction).multiplyScalar(directionDistance).add(
      this.origin);
    return v1.distanceTo(point);
  }

  isIntersectionCylinder() {

  }

  isIntersectionSphere(sphere) {
    return (this.distanceToPoint(sphere.center) <= sphere.radius);

  }

  isIntersectionPlane(plane) {

    var denominator = plane.normal.dot(this.direction);

    // plane and ray are not perpendicular
    if (denominator !== 0)
      return true;

    if (plane.distanceToPoint(this.origin) === 0)
      return true;

    return false;

  }

  distanceToPlane(plane) {

    var denominator = plane.normal.dot(this.direction);
    if (denominator === 0) {

      // line is coplanar
      if (plane.distanceToPoint(this.origin) === 0)
        return 0;

      // ray is parallel
      return undefined;
    }

    var t = -(this.origin.dot(plane.normal) + plane.constant) / denominator;

    return t;

  }

  intersectPlane(plane, optionalTarget) {

    var t = this.distanceToPlane(plane);

    if (t === undefined)
      return undefined;

    return this.at(t, optionalTarget);

  }

  applyMatrix4(matrix4) {

    this.direction.add(this.origin).applyMatrix4(matrix4);
    this.origin.applyMatrix4(matrix4);
    this.direction.sub(this.origin);

    return this;

  }

  equals(ray) {

    return ray.origin.equals(this.origin)
      && ray.direction.equals(this.direction);

  }

  clone() {

    return new Ray().copy(this);
  }
}