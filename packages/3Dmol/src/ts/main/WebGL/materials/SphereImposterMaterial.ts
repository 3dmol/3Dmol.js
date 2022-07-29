import { ImposterMaterial } from "./ImposterMaterial";

export class SphereImposterMaterial extends ImposterMaterial {
  shaderID = "sphereimposter";
  constructor(parameters?: any) {
    super(parameters);
    this.setValues(parameters);
  }

  clone<T extends this>(material: T = new SphereImposterMaterial() as T): T {
    super.clone.call(this, material);
    return material;
  }
}
