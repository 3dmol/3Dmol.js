import { Object3D } from "./core/Object3D";
import { Color } from "./core/Color";
import { Vector3 } from "./math";

export class Light extends Object3D {
  color: Color;
  target: Object3D;
  intensity: any;
  castShadow: boolean;
  onlyShadow: boolean;
  constructor(hex, intensity) {
    super();

    this.color = new Color(hex);
    this.position = new Vector3(0, 1, 0);
    this.target = new Object3D();

    this.intensity = intensity !== undefined ? intensity : 1;

    this.castShadow = false;
    this.onlyShadow = false;
  }
}
