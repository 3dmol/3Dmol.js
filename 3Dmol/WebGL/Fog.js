import Color from './core/Color';

/*
 * Fog Class
 */
export default class Fog {
  constructor(hex, near, far) {
    this.name = '';

    this.color = new Color(hex);

    this.near = near !== undefined ? near : 1;
    this.far = far !== undefined ? far : 1000;
  }

  clone() {
    return new Fog(this.color.getHex(), this.near, this.far);
  }
}
