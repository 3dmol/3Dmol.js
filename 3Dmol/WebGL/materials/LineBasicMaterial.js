// @ts-check

import { Color } from "../core/Color";
import { Material } from "./Material";

//Line basic material
export class LineBasicMaterial extends Material {
  color = new Color(0xffffff);
  linewidth = 1;
  linecap = "round";
  linejoin = "round";
  vertexColors = false;
  fog = true;
  shaderID = "basic";
  /**
   * @param {Record<any,any>} [parameters]
   */
  constructor(parameters) {
    super();
    this.setValues(parameters || {});
  }

  clone() {
    var material = new LineBasicMaterial();

    super.clone.call(this, material);

    material.color.copy(this.color);
    return material;
  }
}
