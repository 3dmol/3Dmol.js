import { Matrix4 } from '../math';
import { Vector3 } from '../math';
//Intersection sphere and box shapes.

//Intersection sphere for sphere, stick render
export class Sphere {
  center: Vector3
  radius: number

  constructor(center = new Vector3(), radius = 0) {
    this.center = center;
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
