import { ImposterMaterial } from "./ImposterMaterial";

export class SphereImposterMaterial extends ImposterMaterial {
  shaderID: string;
  constructor(parameters?: any) {
    super(parameters);

    this.shaderID = "sphereimposter";
    this.setValues(parameters);
  }

  clone() {
    var material = new SphereImposterMaterial();
    super.clone.call(this, material);
    return material;
  }
}
