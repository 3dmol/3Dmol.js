import { Color } from "../../colors";
import { Material } from "./Material";
//Outlined Mesh Lamert material
/* @constructor */
export class MeshOutlineMaterial extends Material {
  fog: boolean;
  shaderID: string;
  wireframe: boolean;
  outlineColor: any;
  outlineWidth: any;
  outlinePushback: any;
  constructor(parameters?: any) {
    super();
    parameters = parameters || {};
    this.fog = true;
    this.shaderID = "outline";
    this.wireframe = false;
    this.outlineColor = parameters.color || new Color(0.0, 0.0, 0.0);
    this.outlineWidth = parameters.width || 0.1;
    this.outlinePushback = parameters.pushback || 1.0;
  }
  clone<T extends this>(material: T = new MeshOutlineMaterial() as T): T {
    super.clone.call(this, material);
    material.fog = this.fog;
    material.shaderID = this.shaderID;
    material.wireframe = this.wireframe;
    return material;
  }
}
