import { Color } from './../core/Color';
import { Material } from "./Material";
//Line basic material
/** @constructor */
export class LineBasicMaterial extends Material {
  color: any;
  linewidth: number;
  linecap: string;
  linejoin: string;
  vertexColors: boolean;
  fog: boolean;
  shaderID: string;
  constructor(parameters?: any) {
    super();

    this.color = new Color(0xffffff);

    this.linewidth = 1;
    this.linecap = "round";
    this.linejoin = "round";

    this.vertexColors = false;

    this.fog = true;
    this.shaderID = "basic";
    this.setValues(parameters);
  }

  clone() {
    var material = new LineBasicMaterial();

    super.clone.call(this, material);

    material.color.copy(this.color);
    return material;
  }
}
