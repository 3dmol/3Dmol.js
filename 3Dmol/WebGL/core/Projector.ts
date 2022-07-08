import { Matrix4 } from './../math';

const viewProjectionMatrix = new Matrix4();

//$3Dmol Projection 
export class Projector {   
  static unprojectVector( vector, camera ) {

    camera.projectionMatrixInverse.getInverse(camera.projectionMatrix);

    viewProjectionMatrix.multiplyMatrices(camera.matrixWorld, camera.projectionMatrixInverse);

    return vector.applyProjection(viewProjectionMatrix );

};

  static projectVector( vector, camera ) {

    camera.matrixWorldInverse.getInverse( camera.matrixWorld );

    viewProjectionMatrix.multiplyMatrices( camera.projectionMatrix, camera.matrixWorldInverse );

    return vector.applyProjection( viewProjectionMatrix );

};
}