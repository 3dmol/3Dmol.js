import { ImposterMaterial } from "./ImposterMaterial";
export class StickImposterMaterial extends ImposterMaterial {
  shaderID = "stickimposter";

  constructor(parameters?: any) {
    super(parameters);
    this.setValues(parameters);
  }

  clone<T extends StickImposterMaterial>(material: T): T;
  clone(material?: StickImposterMaterial): StickImposterMaterial;
  clone(material: StickImposterMaterial = new StickImposterMaterial()): StickImposterMaterial {
    super.clone.call(this, material);
    return material;
  }
}
