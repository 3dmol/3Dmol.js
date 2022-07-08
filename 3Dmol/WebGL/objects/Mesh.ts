import { MeshBasicMaterial } from "../materials/MeshBasicMaterial";
import { Object3D } from "./../core/Object3D";
//Mesh Object
/** @constructor */
export class Mesh extends Object3D {
  constructor(geometry, material) {
    super();

    this.geometry = geometry;
    this.material =
      material !== undefined
        ? material
        : new MeshBasicMaterial({
            color: Math.random() * 0xffffff,
            wireframe: true,
          });
  }

  clone(object) {
    if (object === undefined)
      object = new Mesh(this.geometry, this.material);

    super.clone.call(this, object);

    return object;
  }
}
