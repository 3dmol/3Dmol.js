import { Color } from "../core/Color";
import { ImposterMaterial } from "./ImposterMaterial";

export class SphereImposterOutlineMaterial extends ImposterMaterial {
  outlineColor: any;
  outlineWidth: any;
  outlinePushback: any;
  
  constructor(parameters?: any) {
    super(parameters);
    parameters = parameters || {};

    this.shaderID = "sphereimposteroutline";
    this.outlineColor = parameters.color || new Color(0.0, 0.0, 0.0);
    this.outlineWidth = parameters.width || 0.1;
    this.outlinePushback = parameters.pushback || 1.0;

    this.setValues(parameters);
  }

  clone() {
    var material = new SphereImposterOutlineMaterial();
    super.clone.call(this, material);
    material.outlineColor = this.outlineColor;
    material.outlineWidth = this.outlineWidth;
    material.outlinePushback = this.outlinePushback;
    return material;
  }
}
