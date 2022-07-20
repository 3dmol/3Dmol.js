import { Material, LineBasicMaterial } from "../materials";
import type { Geometry } from "../core";
import { Object3D } from "../core";

export enum LineStyle {
  LineStrip = 0,
  LinePieces = 1,
}

export class Line extends Object3D {
  type: any;
  geometry: Geometry;
  material: Material;
  constructor(
    geometry: Geometry,
    material: Material = new LineBasicMaterial({
      color: Math.random() * 0xffffff,
    }) as Material,
    type: LineStyle = LineStyle.LineStrip
  ) {
    super();

    this.geometry = geometry;
    //TODO: update material and type to webgl
    this.material = material;
    this.type = type;
  }

  clone<T extends this>(
    object = new Line(this.geometry, this.material, this.type) as T
  ): T {
    super.clone.call(this, object);
    return object;
  }
}


