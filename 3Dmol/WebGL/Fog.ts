/*
 * Fog Class
 */

import { Color } from "./core/Color";

/** @constructor */
export class Fog {
  name: string;
  color: any;
  near: any;
  far: any;
  constructor(hex, near, far) {
    this.name = "";

    this.color = new Color(hex);

    this.near = near !== undefined ? near : 1;
    this.far = far !== undefined ? far : 1000;
  }

  clone() {
    return new Fog(this.color.getHex(), this.near, this.far);
  }
}
