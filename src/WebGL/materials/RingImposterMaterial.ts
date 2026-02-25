import { ImposterMaterial } from "./ImposterMaterial";

export class RingImposterMaterial extends ImposterMaterial {
  shaderID = "ringimposter";
  constructor(parameters?: any) {
    super(parameters);
    this.setValues(parameters);
  }

  clone<T extends this>(material = new RingImposterMaterial() as T): T {
    super.clone.call(this, material);
    return material;
  }
}
