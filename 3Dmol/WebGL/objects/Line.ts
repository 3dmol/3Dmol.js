import { LineBasicMaterial } from '../materials/LineBasicMaterial';
import { Object3D } from "../core/Object3D";
//Line Object
/** @constructor */
export class Line extends Object3D {
  type: any;
  constructor(geometry, material, type) {
    super();

    this.geometry = geometry;
    //TODO: update material and type to webgl
    this.material =
      material !== undefined
        ? material
        : new LineBasicMaterial({ color: Math.random() * 0xffffff });
    this.type = type !== undefined ? type : LineStyle.LineStrip;
  }

  clone(object) {
    if (object === undefined)
      object = new Line(this.geometry, this.material, this.type);

    Object3D.prototype.clone.call(this, object);

    return object;
  }
}

export enum LineStyle {
  LineStrip = 0,
  LinePieces = 1,
}
