//Volumetric material

import { FrontSide } from "../constants/Sides";
import { Color } from "../../colors";
import { Material } from "./Material";

/* @constructor */
export class VolumetricMaterial extends Material {
  transparent = false;
  volumetric = true;
  color = new Color(0xffffff);
  transferfn = null;
  map = undefined;
  extent = [];
  maxdepth = 100.0;
  unit = 0;
  texmatrix = null;
  transfermin = -1.0;
  transfermax = 1.0;
  subsamples = 5.0;
  shaderID = "volumetric";
  side = FrontSide;

  constructor(parameters?: any) {
    super();
    // this.fog = true; // TODO: to integrate the new shader with the fog stuff
    this.setValues(parameters);
  }
  clone<T extends this>(material = new VolumetricMaterial() as T): T {
    super.clone.call(this, material);
    material.transparent = this.transparent;
    material.volumetric = this.volumetric;
    material.color = this.color;
    material.transferfn = this.transferfn;
    material.map = this.map;
    material.extent = this.extent;
    material.maxdepth = this.maxdepth;
    material.unit = this.unit;
    material.texmatrix = this.texmatrix;
    material.transfermin = this.transfermin;
    material.transfermax = this.transfermax;
    material.subsamples = this.subsamples;
    material.shaderID = this.shaderID;
    material.side = this.side;
    return material;
  }
}
