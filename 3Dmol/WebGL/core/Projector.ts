import type { Vector3 } from '../math';
import type { Camera } from '../Camera';
import { Matrix4 } from '../math';

const viewProjectionMatrix = new Matrix4();

//$3Dmol Projection 
export class Projector {
  static unprojectVector(vector: Vector3, camera: Camera) {

    camera.projectionMatrixInverse.getInverse(camera.projectionMatrix);

    viewProjectionMatrix.multiplyMatrices(camera.matrixWorld, camera.projectionMatrixInverse);

    return vector.applyProjection(viewProjectionMatrix);

  };

  static projectVector(vector: Vector3, camera: Camera) {

    camera.matrixWorldInverse.getInverse(camera.matrixWorld);

    viewProjectionMatrix.multiplyMatrices(camera.projectionMatrix, camera.matrixWorldInverse);

    return vector.applyProjection(viewProjectionMatrix);
  };

  projectVector(vector: Vector3, camera: Camera) {

    return Projector.projectVector(vector, camera);

  }

  unprojectVector(vector: Vector3, camera: Camera) {

    return Projector.unprojectVector(vector, camera);

  }
}