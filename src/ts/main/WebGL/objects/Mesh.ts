import type { Material } from "../materials";
import type { Geometry } from "../core";
import { MeshBasicMaterial } from "../materials";
import { Object3D } from "../core";
//Mesh Object
/** @constructor */
export class Mesh extends Object3D {
  geometry: Geometry;
  material: Material;
  constructor(
    geometry: Geometry,
    material = new MeshBasicMaterial({
      color: Math.random() * 0xffffff,
      wireframe: true,
    }) as Material
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
