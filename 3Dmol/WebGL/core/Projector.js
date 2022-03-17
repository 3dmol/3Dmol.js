// @ts-check

import { Matrix4 } from "../math/Matrix4";
import { Vector3 } from "../math/Vector3";
import { Camera } from "../objects";

//$3Dmol Projection
export class Projector {
  _viewProjectionMatrix = new Matrix4();

  /**
   * @param {Vector3} vector
   * @param {Camera} camera
   */
  projectVector(vector, camera) {
    camera.matrixWorldInverse.getInverse(camera.matrixWorld);

    this._viewProjectionMatrix.multiplyMatrices(
      camera.projectionMatrix,
      camera.matrixWorldInverse
    );

    return vector.applyProjection(this._viewProjectionMatrix);
  }

  /**
   * @param {Vector3} vector
   * @param {Camera} camera
   */
  unprojectVector(vector, camera) {
    camera.projectionMatrixInverse.getInverse(camera.projectionMatrix);

    this._viewProjectionMatrix.multiplyMatrices(
      camera.matrixWorld,
      camera.projectionMatrixInverse
    );

    return vector.applyProjection(this._viewProjectionMatrix);
  }
}
