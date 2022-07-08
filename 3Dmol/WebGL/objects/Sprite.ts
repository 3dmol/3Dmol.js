import { SpriteMaterial } from './../materials/SpriteMaterial';
import { Object3D } from "./../core/Object3D";
import { Vector3 } from '../math';
//Sprite object
/** @constructor */
export class Sprite extends Object3D {
  rotation3d: Vector3;
  constructor(material) {
    super();

    this.material = material !== undefined ? material : new SpriteMaterial();

    this.rotation3d = this.rotation as Vector3;
    this.rotation = 0;
  }

  updateMatrix() {
    this.matrix.setPosition(this.position);

    this.rotation3d.set(0, 0, this.rotation);
    this.matrix.setRotationFromEuler(this.rotation3d);

    if (this.scale.x !== 1 || this.scale.y !== 1) this.matrix.scale(this.scale);

    this.matrixWorldNeedsUpdate = true;
  }

  clone(object) {
    if (object === undefined) object = new Sprite(this.material);

    Object3D.prototype.clone.call(this, object);

    return object;
  }
}
