import { ImposterMaterial } from "./ImposterMaterial";
export class StickImposterMaterial extends ImposterMaterial {
  shaderID = "stickimposter";

  constructor(parameters?: Record<keyof ImposterMaterial, unknown>) {
    super(parameters);
    this.setValues(parameters);
  }

  clone<T extends this>(material = new StickImposterMaterial() as T): T {
    super.clone.call(this, material);
    return material;
  }
}
