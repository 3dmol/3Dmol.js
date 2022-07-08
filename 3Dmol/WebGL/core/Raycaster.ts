import { Ray, Matrix4 } from "../math";
import { intersectObject } from "./intersectObject";

const descSort = (a, b) => {
  return a.distance - b.distance;
};



const viewProjectionMatrix = new Matrix4();

export class Raycaster {
  ray: Ray;
  near: number;
  far: number;
  precision = 0.0001;
  linePrecision = 0.2;
  constructor(origin, direction, far, near) {
    this.ray = new Ray(origin, direction);

    if (this.ray.direction.lengthSq() > 0) this.ray.direction.normalize();

    this.near = near || 0;
    this.far = far || Infinity;
  }

  set(origin, direction): void {
    this.ray.set(origin, direction);
  }

  setFromCamera(coords, camera): void {
    if (!camera.ortho) {
      this.ray.origin.setFromMatrixPosition(camera.matrixWorld);
      this.ray.direction.set(coords.x, coords.y, coords.z);

      camera.projectionMatrixInverse.getInverse(camera.projectionMatrix);
      viewProjectionMatrix.multiplyMatrices(
        camera.matrixWorld,
        camera.projectionMatrixInverse
      );
      this.ray.direction.applyProjection(viewProjectionMatrix);
      this.ray.direction.sub(this.ray.origin).normalize();
    } else {
      this.ray.origin
        .set(
          coords.x,
          coords.y,
          (camera.near + camera.far) / (camera.near - camera.far)
        )
        .unproject(camera);
      this.ray.direction.set(0, 0, -1).transformDirection(camera.matrixWorld);
    }
  }

  intersectObjects(group, objects) {
    var intersects = [];
  
    for (var i = 0, l = objects.length; i < l; i++)
      intersectObject(group, objects[i], this, intersects);
  
    intersects.sort(descSort);
  
    return intersects;
  }
}