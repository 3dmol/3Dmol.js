import type { Matrix4 } from '../math';
import { Vector3 } from '../math';
let vector = new Vector3()

//Bounding cylinder for stick render
/** @constructor */
export class Cylinder {
  c1: Vector3
  c2: Vector3
  direction: Vector3
  radius: number

  constructor(c1?:Vector3, c2?:Vector3, radius?:number) {
    this.c1 = c1 !== undefined ? c1 : new Vector3();

    this.c2 = c2 !== undefined ? c2 : new Vector3();

    this.direction = new Vector3()
      .subVectors(this.c2, this.c1)
      .normalize();

    this.radius = radius !== undefined ? radius : 0;
  }

  copy(cylinder:Cylinder):Cylinder {
    this.c1.copy(cylinder.c1);
    this.c2.copy(cylinder.c2);
    this.direction.copy(cylinder.direction);
    this.radius = cylinder.radius;
    return this;
  }

  lengthSq():number { 
    return vector.subVectors(this.c2, this.c1).lengthSq();
  }

  applyMatrix4(matrix:Matrix4):Cylinder {
    this.direction.add(this.c1).applyMatrix4(matrix);
    this.c1.applyMatrix4(matrix);
    this.c2.applyMatrix4(matrix);
    this.direction.sub(this.c1).normalize();
    this.radius = this.radius * matrix.getMaxScaleOnAxis();

    return this;
  }
}
