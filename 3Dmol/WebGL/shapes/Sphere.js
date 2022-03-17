// @ts-check

import { Vector3 } from "../math/Vector3";

//Intersection sphere for sphere, stick render
export class Sphere {
  constructor(center, radius) {
    this.center = center !== undefined ? center : new Vector3();
    this.radius = radius !== undefined ? radius : 0;
  }

  set(center, radius) {
    this.center.copy(center);
    this.radius = radius;

    return this;
  }

  copy(sphere) {
    this.center.copy(sphere.center);
    this.radius = sphere.radius;

    return this;
  }

  applyMatrix4(matrix) {
    this.center.applyMatrix4(matrix);
    this.radius = this.radius * matrix.getMaxScaleOnAxis();

    return this;
  }

  translate(offset) {
    this.center.add(offset);

    return this;
  }

  equals(sphere) {
    return sphere.center.equals(this.center) && sphere.radius === this.radius;
  }

  clone() {
    return new Sphere().copy(this);
  }
}
