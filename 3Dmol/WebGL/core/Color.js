// @ts-check
export class Color {
  /**
     * @param {number | undefined} [r]
     * @param {number | undefined} [g]
     * @param {number | undefined} [b]
     */
  constructor(r,g,b) {
    if (arguments.length > 1) {
      this.r = arguments[0] || 0.0;
      this.g = arguments[1] || 0.0;
      this.b = arguments[2] || 0.0;

      return this;
    }
    this.r = 0;
    this.g = 0;
    this.b = 0;
  }

  /**
   * @param {{ clone: () => any;}| number |{r: any; g: any; b: any; }} val
   */
  set(val) {
    if (val instanceof Color) return val.clone();
    else if (typeof val === "number") this.setHex(val);
    else if (
      typeof val === "object" &&
      "r" in val &&
      "g" in val &&
      "b" in val
    ) {
      this.r = val.r;
      this.g = val.g;
      this.b = val.b;
    }
  }

  /**
   * @param {number} hex
   */
  setHex(hex) {
    hex = Math.floor(hex);

    this.r = ((hex >> 16) & 255) / 255;
    this.g = ((hex >> 8) & 255) / 255;
    this.b = (hex & 255) / 255;

    return this;
  }

  getHex() {
    var R = Math.round(this.r * 255);
    var G = Math.round(this.g * 255);
    var B = Math.round(this.b * 255);
    return (R << 16) | (G << 8) | B;
  }

  clone() {
    return new Color(this.r, this.g, this.b);
  }

  /**
   * @param {Color} color
   */
  copy(color) {
    this.r = color.r;
    this.g = color.g;
    this.b = color.b;

    return this;
  }

  scaled() {
    var ret = {};
    ret.r = Math.round(this.r * 255);
    ret.g = Math.round(this.g * 255);
    ret.b = Math.round(this.b * 255);
    ret.a = 1.0;
    return ret;
  }
}
