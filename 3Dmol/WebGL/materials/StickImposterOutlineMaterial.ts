import { Color } from "../core/Color";
import { ImposterMaterial } from "./ImposterMaterial";
export class StickImposterOutlineMaterial extends ImposterMaterial {
  outlineColor: any;
  outlineWidth: any;
  outlinePushback: any;
  constructor(parameters?: Record<string, any> & { color?: Color; width?: number; pushback?: number }) {
    super(parameters);
    parameters = parameters || {};

    this.shaderID = "stickimposteroutline";
    this.outlineColor = parameters.color || new Color(0.0, 0.0, 0.0);
    this.outlineWidth = parameters.width || 0.1;
    this.outlinePushback = parameters.pushback || 1.0;

    this.setValues(parameters);
  }

  clone(): StickImposterOutlineMaterial {
    var material = new StickImposterOutlineMaterial();
    super.clone.call(this, material);
    material.outlineColor = this.outlineColor;
    material.outlineWidth = this.outlineWidth;
    material.outlinePushback = this.outlinePushback;
    return material;
  }
}
