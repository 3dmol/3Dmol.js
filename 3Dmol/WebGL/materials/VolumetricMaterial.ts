//Volumetric material

import { Sides } from "../constants/Sides";
import { Color } from "../core/Color";
import { Material } from "./Material";

/** @constructor */
export class VolumetricMaterial extends Material {
  transparent: boolean;
  volumetric: boolean;
  color: any;
  transferfn: null;
  map: null;
  extent: never[];
  maxdepth: number;
  unit: number;
  texmatrix: null;
  transfermin: number;
  transfermax: number;
  subsamples: number;
  shaderID: string;
  side: number;
  constructor(parameters?: any) {
    super();

    this.transparent = false;
    this.volumetric = true;

    this.color = new Color(0xffffff);
    this.transferfn = null;
    this.map = null;
    this.volumetric = true;
    this.extent = [];
    this.maxdepth = 100.0;
    this.unit = 0;
    this.texmatrix = null;
    this.transfermin = -1.0;
    this.transfermax = 1.0;
    this.subsamples = 5.0;

    // this.fog = true; // TODO: to integrate the new shader with the fog stuff

    this.shaderID = "volumetric";
    this.side = Sides.FrontSide;

    this.setValues(parameters);
  }
;

clone() {

  var material = Object.assign(new VolumetricMaterial(), this);

  super.clone.call(this, material);
  return material;

};
}