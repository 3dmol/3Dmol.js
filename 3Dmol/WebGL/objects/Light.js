// @ts-check

import { Color, Object3D } from "../core/index";
import { Vector3 } from "../math/index";


export class Light extends Object3D {
  constructor(hex, intensity) {
    super();

    this.color = new Color(hex);
    this.position = new Vector3(0, 1, 0);
    this.target = new Object3D();

    this.intensity = (intensity !== undefined) ? intensity : 1;

    this.castShadow = false;
    this.onlyShadow = false;

  }
};
