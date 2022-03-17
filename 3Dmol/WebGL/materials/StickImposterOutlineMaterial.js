// @ts-check

import { Color } from "../core/Color";
import { ImposterMaterial } from "./ImposterMaterial";

export class StickImposterOutlineMaterial extends ImposterMaterial {
    constructor(parameters) {
  
      super();
      parameters = parameters || {};
  
      this.shaderID = "stickimposteroutline";
      this.outlineColor = parameters.color || new Color(0.0, 0.0, 0.0);
      this.outlineWidth = parameters.width || 0.1;
      this.outlinePushback = parameters.pushback || 1.0;
  
      this.setValues(parameters);
  
    }
  
    clone() {
  
      var material = new StickImposterOutlineMaterial();
      super.clone.call(this, material);
      material.outlineColor = this.outlineColor;
      material.outlineWidth = this.outlineWidth;
      material.outlinePushback = this.outlinePushback;
      return material;
    }
  }