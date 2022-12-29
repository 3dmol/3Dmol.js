import type { Matrix4 } from "../math";
import { Vector3 } from "../math";

const v1 = new Vector3();


/**   plane specified by three points

 *  @class 
 *  @subcategory  Math
 *  
 * 
 */ 
export class Triangle {
  a: Vector3;
  b: Vector3;
  c: Vector3;

  constructor(a = new Vector3(), b = new Vector3(), c = new Vector3()) {
    this.a = a;
    this.b = b;
    this.c = c;
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
