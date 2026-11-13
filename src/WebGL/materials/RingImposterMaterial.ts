import { ImposterMaterial } from "./ImposterMaterial";

export class RingImposterMaterial extends ImposterMaterial {
  shaderID = "ringimposter";
  constructor(parameters?: any) {
    super(parameters);
    this.setValues(parameters);
  }

  clone<T extends RingImposterMaterial>(material: T): T;
  clone(material?: RingImposterMaterial): RingImposterMaterial;
  clone(material: RingImposterMaterial = new RingImposterMaterial()): RingImposterMaterial {
    super.clone.call(this, material);
    return material;
  }
}
