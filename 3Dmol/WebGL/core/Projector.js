// @ts-check

import { Matrix4 } from "../math";

//$3Dmol Projection 
export class Projector {

  _viewProjectionMatrix = new Matrix4();

  projectVector = function (vector, camera) {

    camera.matrixWorldInverse.getInverse(camera.matrixWorld);

    this._viewProjectionMatrix.multiplyMatrices(camera.projectionMatrix, camera.matrixWorldInverse);

    return vector.applyProjection(this._viewProjectionMatrix);

  };

  unprojectVector = function (vector, camera) {

    camera.projectionMatrixInverse.getInverse(camera.projectionMatrix);

    this._viewProjectionMatrix.multiplyMatrices(camera.matrixWorld, camera.projectionMatrixInverse);

    return vector.applyProjection(this._viewProjectionMatrix);

  };

};