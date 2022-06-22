export interface Colored {
  r: number;
  g: number;
  b: number;
  a?: number;
}

export class Color implements Colored {
  r: number = 0.0;
  g: number = 0.0;
  b: number = 0.0;

  constructor(r: number | Color | Colored, g?: number, b?: number) {
    if (arguments.length > 1 && typeof r === "number") {
      this.r = r || 0.0;
      this.g = g || 0.0;
      this.b = b || 0.0;

      return this;
    }

    return this.set(r);
  }

  set<T extends Colored>(val: Color | number | T): Color {
    if (val instanceof Color) return val.clone();
    else if (typeof val === "number") this.setHex(val);
    else if (typeof val === "object") {
      this.r = val?.r || 0.0;
      this.g = val?.g || 0.0;
      this.b = val?.b || 0.0;
    }
    return this;
  }

  setHex(hex: number): Color {
    hex = Math.floor(hex);

    this.r = ((hex >> 16) & 255) / 255;
    this.g = ((hex >> 8) & 255) / 255;
    this.b = (hex & 255) / 255;

    return this;
  }

  getHex(): number {
    var R = Math.round(this.r * 255);
    var G = Math.round(this.g * 255);
    var B = Math.round(this.b * 255);
    return (R << 16) | (G << 8) | B;
  }

  clone(): Color {
    return new Color(this.r, this.g, this.b);
  }

  copy(color: Color): Color {
    this.r = color.r;
    this.g = color.g;
    this.b = color.b;

    return this;
  }

  //return object that represents color components from 0 to 255
  scaled<T extends Colored>(): Colored {
    var ret: Partial<T> = {};
    ret.r = Math.round(this.r * 255);
    ret.g = Math.round(this.g * 255);
    ret.b = Math.round(this.b * 255);
    ret.a = 1.0;
    return ret as T;
  }
}