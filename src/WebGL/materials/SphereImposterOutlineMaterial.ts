import { Color } from "../../colors";
import { ImposterMaterial } from "./ImposterMaterial";

export class SphereImposterOutlineMaterial extends ImposterMaterial {
  outlineColor: Color;
  outlineWidth: number;
  outlinePushback: number;
  
  constructor(parameters?: any) {
    super(parameters);
    parameters = parameters || {};

    this.shaderID = "sphereimposteroutline";
    this.outlineColor = parameters.color || new Color(0.0, 0.0, 0.0);
    this.outlineWidth = parameters.width || 0.1;
    this.outlinePushback = parameters.pushback || 1.0;

    this.setValues(parameters);
  }

  clone<T extends this>(material: T = new SphereImposterOutlineMaterial() as T): T {
    super.clone.call(this, material);
    material.outlineColor = this.outlineColor;
    material.outlineWidth = this.outlineWidth;
    material.outlinePushback = this.outlinePushback;
    return material;
  }
}
