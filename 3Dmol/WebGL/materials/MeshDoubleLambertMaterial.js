// @ts-check

import { MeshLambertMaterial } from "./MeshLambertMaterial";
import { DoubleSide } from "./sides";

//Double sided Mesh Lambert material
export class MeshDoubleLambertMaterial extends MeshLambertMaterial {
  shaderID = "lambertdouble";
  side = DoubleSide;

  clone() {
    var material = new MeshDoubleLambertMaterial();

    super.clone.call(this, material);

    return material;
  }
}
