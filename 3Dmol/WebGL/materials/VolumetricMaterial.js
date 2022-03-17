// @ts-check

import { Color } from "../core/Color";
import { Material } from "./Material";
import { FrontSide } from "./sides";

//Volumetric material
export class VolumetricMaterial extends Material {
  constructor(parameters) {

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
    this.side = FrontSide;

    this.setValues(parameters);
  }

  clone() {

    var material = Object.assign(new VolumetricMaterial(), this);

    super.clone.call(this, material);
    return material;

  }
}