import type { Material } from "../materials";
import type { Geometry } from "../core";
import { Object3D } from "../core";
//Mesh Object
/* @constructor */
export class Mesh extends Object3D {
  geometry: Geometry;
  material: Material;
  constructor(
    geometry: Geometry,
    material: Material
  ) {
    super();
    this.geometry = geometry;
    this.material = material;
  }

  clone(object: Mesh): Mesh {
    if (object === undefined) object = new Mesh(this.geometry, this.material);

    super.clone.call(this, object);

    return object;
  }
}
