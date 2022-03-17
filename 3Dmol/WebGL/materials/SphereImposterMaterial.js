// @ts-check
import { ImposterMaterial } from "./ImposterMaterial";

export class SphereImposterMaterial extends ImposterMaterial {
  shaderID = "sphereimposter";

  /**
   * @param {Record<any, any>} [parameters]
   */
  constructor(parameters) {
    super();
    this.setValues(parameters || {});
  }

  /**
   * @param {SphereImposterMaterial} material
   */
  clone(material) {
    var material = material || new SphereImposterMaterial();
    super.clone.call(this, material);
    return material;
  }
}
