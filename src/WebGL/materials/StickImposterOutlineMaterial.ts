import { CC, Color, ColorSpec } from "../../colors";
import { ImposterMaterial } from "./ImposterMaterial";
export class StickImposterOutlineMaterial extends ImposterMaterial {
  shaderID = "stickimposteroutline";
  outlineColor = new Color(0.0, 0.0, 0.0);
  outlineWidth = 0.1;
  outlinePushback = 1.0;
  outlineMaxPixels = 0.0;

  constructor(
    parameters: Record<keyof StickImposterOutlineMaterial, unknown> & {
      width?: number;
      pushback?: number;
      maxpixels?: number;
    } = {} as any
  ) {
    super(parameters);
    if (parameters.color) this.outlineColor = CC.color(parameters.color as ColorSpec);
    if (parameters.width) this.outlineWidth = parameters.width as number;
    if (parameters.pushback)
      this.outlinePushback = parameters.pushback as number;
    if (parameters.maxpixels)
      this.outlineMaxPixels = parameters.maxpixels;

    this.setValues(parameters);
  }

  clone<T extends this>(material = new StickImposterOutlineMaterial() as T): T {
    super.clone.call(this, material);
    material.outlineColor = this.outlineColor;
    material.outlineWidth = this.outlineWidth;
    material.outlinePushback = this.outlinePushback;
    material.outlineMaxPixels = this.outlineMaxPixels;
    return material;
  }
}
