import { DoubleSide } from "./../constants/Sides";
import { MeshLambertMaterial } from "./MeshLambertMaterial";
//Double sided Mesh Lambert material
/* @constructor */
export class MeshDoubleLambertMaterial extends MeshLambertMaterial {
  shaderID = "lambertdouble";
  side = DoubleSide;
  outline = false;
  
  constructor(parameters?: any) {
    super(parameters);
  }

  clone<T extends MeshDoubleLambertMaterial>(material: T): T;
  clone(material?: MeshDoubleLambertMaterial): MeshDoubleLambertMaterial;
  clone(material: MeshDoubleLambertMaterial = new MeshDoubleLambertMaterial()): MeshDoubleLambertMaterial {
    super.clone.call(this, material);
    return material;
  }
}
