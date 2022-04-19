import {Vector3} from '../WebGL/math';

/**
 *
 * @param {import("../specs").Vector3Like} [vector3Like]
 * @returns {import("../WebGL/math").Vector3}
 */
export default function vector3LikeToVector3(vector3Like) {
  if (!vector3Like) return new Vector3();
  return vector3Like instanceof Vector3
    ? vector3Like
    : new Vector3(vector3Like.x, vector3Like.y, vector3Like.z);
}
