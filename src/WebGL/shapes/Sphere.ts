import { Matrix4 } from '../math';
import { Vector3, XYZ } from '../math';
//Intersection sphere and box shapes.

//Intersection sphere for sphere, stick render
/** @class 
 *  @subcategory  Math
 * */ 
export class Sphere {
  center: Vector3
  radius: number
  box?: any;

  constructor(center: XYZ = {x:0,y:0,z:0}, radius = 0) {
    this.center = new Vector3(center.x,center.y,center.z);
    this.radius = radius;
  }

  set(center: Vector3, radius: number): Sphere {
    this.center.copy(center);
    this.radius = radius;
    return this;
  }

  copy(sphere: Sphere): Sphere {
    this.center.copy(sphere.center);
    this.radius = sphere.radius;
    return this;
  }

  applyMatrix4(matrix: Matrix4): Sphere {
    this.center.applyMatrix4(matrix);
    this.radius = this.radius * matrix.getMaxScaleOnAxis();
    return this;
  }

  translate(offset: Vector3): Sphere {
    this.center.add(offset);
    return this;
  }

  equals(sphere: Sphere): boolean {
    return sphere.center.equals(this.center) && sphere.radius === this.radius;
  }

  clone(): Sphere {
    return new Sphere().copy(this);
  }
}
