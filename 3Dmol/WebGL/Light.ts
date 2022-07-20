import { Color, ColorConstructorArg, Object3D } from "./core";
import { Vector3 } from "./math";

export class Light extends Object3D {
  color: Color;
  intensity: any;
  position = new Vector3(0, 1, 0);
  target = new Object3D();
  castShadow = false;
  onlyShadow = false;
  constructor(hex?: ColorConstructorArg, intensity: number = 1) {
    super();
    this.color = new Color(hex);
    this.intensity = intensity;
  }
}
