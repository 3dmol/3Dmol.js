import { Material } from './../materials/Material';
import { SpriteMaterial } from '../materials';
import { Object3D } from "../core";
import { Vector3 } from '../math';
//Sprite object
/* @constructor */
export class Sprite extends Object3D {
  rotation3d: Vector3;
  _modelViewMatrix: any;
  z?: number;
  material?: Material;

  constructor(material = new SpriteMaterial() as Material) {
    super();
    this.material = material as Material;
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

  clone<T extends this>(object = new Sprite(this.material)): Sprite  {
    Object3D.prototype.clone.call(this, object);
    return object;
  }
}
