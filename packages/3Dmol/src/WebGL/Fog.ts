import { Color, ColorConstructorArg } from '../colors';

export class Fog {
  name = "";
  color: Color;
  near: number;
  far: number;
  constructor(hex?:ColorConstructorArg, near = 1, far = 1000) {
    this.color = new Color(hex);
    this.near = near;
    this.far = far;
  }

  clone() {
    return new Fog(this.color.getHex(), this.near, this.far);
  }
}
