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



// in an attempt to reduce memory overhead, cache all $3Dmol.Colors
// this makes things a little faster
export class CC {
  static rgbRegEx =  /rgb(a?)\(\s*([^ ,\)\t]+)\s*,\s*([^ ,\)\t]+)\s*,\s*([^ ,\)\t]+)/i;
  static cache: Record<number, Color> = {0:new Color(0)};
  static color(hex) {
      // Undefined values default to black
      if(!hex)
          return CC.cache[0];
      // cache hits
      if(typeof(CC.cache[hex]) !== "undefined") {
          return CC.cache[hex];
      }
      // arrays
      else if(hex && hex.constructor === Array) {
          // parse elements recursively
          return hex.map(CC.color);
      }
      // numbers and hex strings
      hex = CC.getHex(hex);
      if(typeof hex === 'number') {
          var c = new Color(hex);
          CC.cache[hex] = c;
          return c;
      } else {
          // pass through $3Dmol.Color & other objects
          return hex;
      }
  }

  static getHex(hex) {
      if (!isNaN(parseInt(hex)))
          return parseInt(hex);        
      else if (typeof(hex) === 'string') {
          hex = hex.trim();
          
          if(hex.length == 4 && hex[0] == '#') {
              hex = '#' + hex[1]+hex[1]+hex[2]+hex[2]+hex[3]+hex[3]; //expand to full hex number
          }
          
          if(hex.length == 7 && hex[0] == '#') {
              return parseInt(hex.substring(1),16);
          } 
          
          let m = CC.rgbRegEx.exec(hex);
          if(m) {
              if(m[1] != "") {
                  console.log("WARNING: Opacity value in rgba ignored.  Specify separately as opacity attribute.");
              }
              let ret = 0;
              for(let i = 2; i < 5; i++) {
                  ret *= 256;
                  let val = m[i].endsWith("%") ? 255*parseFloat(m[i])/100 : parseFloat(m[i]);
                  ret += Math.round(val);
              }
              return ret;
          }
          return (window as any)?.$3Dmol?.htmlColors[hex.toLowerCase()] || 0x000000;
          
      }
      return hex;
  }
  
};