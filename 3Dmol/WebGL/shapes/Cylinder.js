// @ts-check

import { Vector3 } from "../math/Vector3";

var vector = new Vector3();

//Bounding cylinder for stick render
export class Cylinder {
  constructor(c1, c2, radius) {
    this.c1 = c1 !== undefined ? c1 : new Vector3();
    this.c2 = c2 !== undefined ? c2 : new Vector3();
    this.direction = new Vector3().subVectors(this.c2, this.c1).normalize();
    this.radius = radius !== undefined ? radius : 0;
  }

  copy(cylinder) {
    this.c1.copy(cylinder.c1);
    this.c2.copy(cylinder.c2);
    this.direction.copy(cylinder.direction);
    this.radius = cylinder.radius;

    return this;
  }

  lengthSq() {
    return vector.subVectors(this.c2, this.c1).lengthSq();
  }

  applyMatrix4(matrix) {
    this.direction.add(this.c1).applyMatrix4(matrix);
    this.c1.applyMatrix4(matrix);
    this.c2.applyMatrix4(matrix);
    this.direction.sub(this.c1).normalize();
    this.radius = this.radius * matrix.getMaxScaleOnAxis();

    return this;
  }
}
