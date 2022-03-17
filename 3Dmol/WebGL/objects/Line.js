// @ts-check

import { Object3D } from "../core/index";
import { LineBasicMaterial } from "../materials/LineBasicMaterial";

export var LineStrip = 0;
export var LinePieces = 1;

//Line Object
export class Line extends Object3D {
  constructor(geometry, material, type) {
    super();

    this.geometry = geometry;
    //TODO: update material and type to webgl
    this.material =
      material !== undefined
        ? material
        : new LineBasicMaterial({ color: Math.random() * 0xffffff });
    this.type = type !== undefined ? type : LineStrip;
  }

  clone(object) {
    if (object === undefined)
      object = new Line(this.geometry, this.material, this.type);

    super.clone.call(this, object);

    return object;
  }
}
