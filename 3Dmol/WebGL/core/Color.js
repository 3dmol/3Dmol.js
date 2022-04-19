export default class Color {
    r = 0.0;
  
    g = 0.0;
  
    b = 0.0;
  
    constructor(...args) {
      if (arguments.length > 1) {
        this.r = args[0] || 0.0;
        this.g = args[1] || 0.0;
        this.b = args[2] || 0.0;
      } else {
        this.set(args[0]);
      }
    }
  
    set(val) {
      if (val instanceof Color) return val.clone();
  
      if (typeof val === 'number') this.setHex(val);
      else if (typeof val === 'object' && 'r' in val && 'g' in val && 'b' in val) {
        this.r = val.r;
        this.g = val.g;
        this.b = val.b;
      }
      return this;
    }
  
    setHex(hex) {
      hex = Math.floor(hex);
  
      this.r = ((hex >> 16) & 255) / 255;
      this.g = ((hex >> 8) & 255) / 255;
      this.b = (hex & 255) / 255;
  
      return this;
    }
  
    getHex() {
      const R = Math.round(this.r * 255);
      const G = Math.round(this.g * 255);
      const B = Math.round(this.b * 255);
      return (R << 16) | (G << 8) | B;
    }
  
    clone() {
      return new Color(this.r, this.g, this.b);
    }
  
    copy(color) {
      this.r = color.r;
      this.g = color.g;
      this.b = color.b;
  
      return this;
    }
  
    // return object that represents color components from 0 to 255
    scaled() {
      const ret = {};
      ret.r = Math.round(this.r * 255);
      ret.g = Math.round(this.g * 255);
      ret.b = Math.round(this.b * 255);
      ret.a = 1.0;
      return ret;
    }
  }