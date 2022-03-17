// @ts-check

import { Object3D } from "../core/index";
import { SpriteMaterial } from "../materials/SpriteMaterial";

//Sprite object
/** @constructor */
// @ts-ignore
export class Sprite extends Object3D {
  /**
   * @param {SpriteMaterial} [material]
   */
  constructor(material) {
    super();

    this.material = material !== undefined ? material : new SpriteMaterial();

    this.rotation3d = this.rotation;
    this.rotation = 0;
  }

  updateMatrix() {
    this.matrix.setPosition(this.position);

    this.rotation3d.set(0, 0, this.rotation);
    this.matrix.setRotationFromEuler(this.rotation3d);

    if (this.scale.x !== 1 || this.scale.y !== 1) this.matrix.scale(this.scale);

    this.matrixWorldNeedsUpdate = true;
  }

  /**
   * @param {Sprite | undefined} [object]
   * ignored because type of rotation is different from object3d parent
   */
  // @ts-ignore
  clone(object) {
    if (object === undefined) object = new Sprite(this.material);

    // todo update the object3d type to reflect a sprites numeric rotation
    // @ts-ignore
    object = super.clone.call(this, object);

    return object;
  }
}
