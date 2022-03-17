// @ts-check

import { Color } from "../core/Color";
import { Material } from "./Material";

//Outlined Mesh Lamert material
export class MeshOutlineMaterial extends Material {
  fog = true;
  shaderID = "outline";
  wireframe = false;
  /**
   * @param {{ color?: any; width?: any; pushback?: any; } | undefined} [parameters]
   */
  constructor(parameters) {
    super();
    parameters = parameters || {};
    this.outlineColor = parameters.color || new Color(0.0, 0.0, 0.0);
    this.outlineWidth = parameters.width || 0.1;
    this.outlinePushback = parameters.pushback || 1.0;
  }

  /**
   * @param {any} material
   */
  clone(material) {
    if (typeof material === "undefined") material = new MeshOutlineMaterial();
    super.clone.call(this, material);
    material.fog = this.fog;
    material.shaderID = this.shaderID;
    material.wireframe = this.wireframe;
    return material;
  }
}
