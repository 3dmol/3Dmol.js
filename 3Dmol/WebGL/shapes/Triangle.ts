import type { Matrix4 } from "../math";
import { Vector3 } from "../math";

const v1 = new Vector3();

//plane specified by three points
export class Triangle {
  a: Vector3;
  b: Vector3;
  c: Vector3;

  constructor(a, b, c) {
    this.a = a !== undefined ? a : new Vector3();

    this.b = b !== undefined ? b : new Vector3();

    this.c = c !== undefined ? c : new Vector3();
  }

  copy(triangle: Triangle): Triangle {
    this.a.copy(triangle.a);
    this.b.copy(triangle.b);
    this.c.copy(triangle.c);

    return this;
  }

  applyMatrix4(matrix: Matrix4): Triangle {
    this.a.applyMatrix4(matrix);
    this.b.applyMatrix4(matrix);
    this.c.applyMatrix4(matrix);

    return this;
  }

  getNormal(): Vector3 {
    var norm = this.a.clone();
    norm.sub(this.b);
    v1.subVectors(this.c, this.b);

    norm.cross(v1);
    norm.normalize();

    return norm;
  }
}
