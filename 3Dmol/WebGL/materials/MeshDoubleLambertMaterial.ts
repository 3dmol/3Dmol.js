import { Sides } from './../constants/Sides';
import { MeshLambertMaterial } from "./MeshLambertMaterial";
//Double sided Mesh Lambert material
/** @constructor */
export class MeshDoubleLambertMaterial extends MeshLambertMaterial {
  constructor(parameters?: any) {
    super(parameters);

    this.shaderID = "lambertdouble";
    this.side = Sides.DoubleSide;
  }

  clone() {
    var material = new MeshDoubleLambertMaterial();

    super.clone.call(this, material);

    return material;
  }
}
