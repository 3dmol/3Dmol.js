import { Color } from "../../colors";
import { Material } from "./Material";
//Line basic material
/* @constructor */
export class LineBasicMaterial extends Material {
  color = new Color(0xffffff);
  linewidth = 1;
  linecap = "round";
  linejoin = "round";
  vertexColors = false;
  fog = true;
  shaderID = "basic";
  constructor(parameters?: any) {
    super();
    this.setValues(parameters);
  }

  clone<T extends this>(material: T = new LineBasicMaterial() as T): T {

    super.clone.call(this, material);

    material.color.copy(this.color);
    return material as T;
  }
}
