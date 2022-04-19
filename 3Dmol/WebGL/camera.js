import {Object3D} from './core';
import {Matrix4} from './math';

/*
 * Simplified Perspective Camera
 */
export default class Camera extends Object3D {
  projectionMatrix = new Matrix4();
  projectionMatrixInverse = new Matrix4();
  matrixWorldInverse = new Matrix4();

  constructor(fov, aspect, near, far, ortho) {
    super();

    this.fov = fov !== undefined ? fov : 50;
    this.aspect = aspect !== undefined ? aspect : 1;
    this.near = near !== undefined ? near : 0.1;
    this.far = far !== undefined ? far : 2000;

    const center = this.position.z;
    this.right = center * Math.tan((Math.PI / 180) * fov);
    this.left = -this.right;
    this.top = this.right / this.aspect;
    this.bottom = -this.top;

    this.ortho = !!ortho;

    this.updateProjectionMatrix();
  }

  lookAt(vector) {
    // Why is the parameter order switched (compared to Object3D)?
    this.matrix.lookAt(this.position, vector, this.up);

    if (this.rotationAutoUpdate) {
      if (this.useQuaternion === false)
        this.rotation.setEulerFromRotationMatrix(this.matrix, this.eulerOrder);
      else this.quaternion.copy(this.matrix.decompose()[1]);
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
      this.projectionMatrix.makePerspective(this.fov, this.aspect, this.near, this.far);
    }

    this.projectionMatrixInverse.getInverse(this.projectionMatrix);
  }
}
