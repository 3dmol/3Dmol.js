import type { Quaternion } from "./math";
import { Object3D } from "./core";
import { Matrix4, Vector3 } from "./math";
/*
 * Simplified Perspective Camera
 */

/* @constructor */
export class Camera extends Object3D {
  projectionMatrix = new Matrix4();
  projectionMatrixInverse = new Matrix4();
  matrixWorldInverse = new Matrix4();
  right: number;
  left: number;
  top: number;
  bottom: number;
  ortho: boolean;
  fov: number;
  aspect: number;
  near: number;
  far: number;
  z: number;
  constructor(fov = 50, aspect = 1, near = 0.1, far = 2000, ortho = false) {
    super();

    this.fov = fov;
    this.aspect = aspect;
    this.near = near;
    this.far = far;

    var center = this.position.z;
    this.right = center * Math.tan((Math.PI / 180) * fov);
    this.left = -this.right;
    this.top = this.right / this.aspect;
    this.bottom = -this.top;

    this.ortho = !!ortho;

    this.updateProjectionMatrix();
  }

  lookAt(vector: Vector3) {
    //Why is the parameter order switched (compared to Object3D)?
    this.matrix.lookAt(this.position, vector, this.up);

    if (this.rotationAutoUpdate) {
      if (this.useQuaternion === false && this.rotation instanceof Vector3) {
        this.rotation.setEulerFromRotationMatrix(this.matrix, this.eulerOrder);
      } else {
        console.error("Unimplemented math operation.");
      }
    }
  }

  updateProjectionMatrix() {
    if (this.ortho) {
      this.projectionMatrix.makeOrthographic(
        this.left,
        this.right,
        this.top,
        this.bottom,
        this.near,
        this.far
      );
    } else {
      this.projectionMatrix.makePerspective(
        this.fov,
        this.aspect,
        this.near,
        this.far
      );
    }

    this.projectionMatrixInverse.getInverse(this.projectionMatrix);
  }
}
