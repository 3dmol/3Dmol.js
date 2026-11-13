import { ImposterMaterial } from "./ImposterMaterial";

export class SphereImposterMaterial extends ImposterMaterial {
  shaderID = "sphereimposter";
  constructor(parameters?: any) {
    super(parameters);
    this.setValues(parameters);
  }

  clone<T extends SphereImposterMaterial>(material: T): T;
  clone(material?: SphereImposterMaterial): SphereImposterMaterial;
  clone(material: SphereImposterMaterial = new SphereImposterMaterial()): SphereImposterMaterial {
    super.clone.call(this, material);
    return material;
  }
}
