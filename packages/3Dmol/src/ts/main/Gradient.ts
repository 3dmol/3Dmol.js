export abstract class GradientType {
  gradient?: string;
  abstract valueToHex(value: number, range: number): number;
  abstract range(): [number, number] | null;
}

export function normalizeValue(
  lo: number,
  hi: number,
  val: number
): { lo: number; hi: number; val: number } {
  if (hi >= lo) {
    if (val < lo) val = lo;
    if (val > hi) val = hi;
    return { lo: lo, hi: hi, val: val };
  } else {
    if (val > lo) val = lo;
    if (val < hi) val = hi;
    //flip the meaning of val, lo, hi
    val = lo - val + hi;
    return { lo: hi, hi: lo, val: val };
  }
}

//return a Gradient object, even if what is specified is descriptive
export function getGradient(grad) {
  if (grad instanceof GradientType) {
    return grad;
  } else if (
    grad.gradient !== undefined &&
    builtinGradients[grad.gradient]
  ) {
    let min = grad.min === undefined ? -1 : grad.min;
    let max = grad.max === undefined ? 1 : grad.max;
    if (grad.mid === undefined) {
      return new builtinGradients[grad.gradient](min, max);
    } else {
      return new builtinGradients[grad.gradient](min, max, grad.mid);
    }
  }
  return grad;
}

/**
 * Color scheme red to white to blue, for charges
 * Reverse gradients are supported when min>max so that the colors are displayed in reverse order.
 */
export class RWB extends GradientType {
  gradient = "RWB";
  min: number;
  max: number;
  mid?: number;
  mult: number;
  constructor(min: number | [number, number], max?: number, mid?: number) {
    super();
    this.mult = 1.0;
    this.mid = mid;
    this.min = min as number;
    this.max = max as number;
    if (typeof max == "undefined" && Array.isArray(min) && min.length >= 2) {
      //we were passed a single range
      this.max = min[1];
      this.min = min[0];
    } else if (!!min && !!max && !Array.isArray(min)) {
      this.min = min;
      this.max = max;
    }
  }

  //return range used for color mapping, null if none set
  range() {
    if (typeof this.min != "undefined" && typeof this.max != "undefined") {
      return [this.min, this.max] as [number, number];
    }
    return null;
  }

  //map value to hex color, range is provided
  valueToHex(val, range) {
    var lo, hi;
    val = this.mult * val; //reverse if necessary
    if (range) {
      lo = range[0];
      hi = range[1];
    } else {
      lo = this.min;
      hi = this.max;
    }

    if (val === undefined) return 0xffffff;

    var norm = normalizeValue(lo, hi, val);
    lo = norm.lo;
    hi = norm.hi;
    val = norm.val;

    var middle = (hi + lo) / 2;
    if (range && typeof range[2] != "undefined") middle = range[2];
    else if (typeof this.mid != "undefined")
      middle = this.mid; //allow user to specify midpoint
    else middle = (lo + hi) / 2;
    var scale, color;

    //scale bottom from red to white
    if (val < middle) {
      scale = Math.floor(255 * Math.sqrt((val - lo) / (middle - lo)));
      color = 0xff0000 + 0x100 * scale + scale;
      return color;
    } else if (val > middle) {
      //form white to blue
      scale = Math.floor(255 * Math.sqrt(1 - (val - middle) / (hi - middle)));
      color = 0x10000 * scale + 0x100 * scale + 0xff;
      return color;
    } else {
      //val == middle
      return 0xffffff;
    }
  }
}

/**
 * rainbow gradient, but without purple to match jmol
 * Reverse gradients are supported when min>max so that the colors are displayed in reverse order.
 */
export class ROYGB extends GradientType {
  gradient = "ROYGB";
  mult: number;
  max: number;
  min: number;
  constructor(min, max) {
    super();
    this.mult = 1.0;
    this.min = min;
    this.max = max;
    if (typeof max == "undefined" && Array.isArray(min) && min.length >= 2) {
      //we were passed a single range
      this.max = min[1];
      this.min = min[0];
    } else if (!!min && !!max && !Array.isArray(min)) {
      this.min = min;
      this.max = max;
    }
  };
  //map value to hex color, range is provided
  valueToHex(val, range) {
    var lo, hi;
    val = this.mult * val;
    if (range) {
      lo = range[0];
      hi = range[1];
    } else {
      lo = this.min;
      hi = this.max;
    }

    if (typeof val == "undefined") return 0xffffff;

    var norm = normalizeValue(lo, hi, val);
    lo = norm.lo;
    hi = norm.hi;
    val = norm.val;

    var mid = (lo + hi) / 2;
    var q1 = (lo + mid) / 2;
    var q3 = (mid + hi) / 2;

    var scale, color;
    if (val < q1) {
      //scale green up, red up, blue down
      scale = Math.floor(255 * Math.sqrt((val - lo) / (q1 - lo)));
      color = 0xff0000 + 0x100 * scale + 0;
      return color;
    } else if (val < mid) {
      //scale red down, green up, blue down
      scale = Math.floor(255 * Math.sqrt(1 - (val - q1) / (mid - q1)));
      color = 0x010000 * scale + 0xff00 + 0x0;
      return color;
    } else if (val < q3) {
      //scale blue up, red down, green up
      scale = Math.floor(255 * Math.sqrt((val - mid) / (q3 - mid)));
      color = 0x000000 + 0xff00 + 0x1 * scale;
      return color;
    } else {
      //scale green down, blue up, red down
      scale = Math.floor(255 * Math.sqrt(1 - (val - q3) / (hi - q3)));
      color = 0x000000 + 0x0100 * scale + 0xff;
      return color;
    }
  };

  //return range used for color mapping, null if none set
  range() {
    if (typeof this.min != "undefined" && typeof this.max != "undefined") {
      return [this.min, this.max] as [number, number];
    }
    return null;
  };
}
/**
 * @constructor Sinebow
 * rainbow gradient with constant saturation, all the way to purple!
 * Reverse gradients are supported when min>max so that the colors are displayed in reverse order.
 * 
 * @example $.get('data/1fas.pqr', function(data){
      viewer.addModel(data, "pqr");
      $.get("data/1fas.cube",function(volumedata){
          viewer.addSurface($3Dmol.SurfaceType.VDW, {
              opacity:0.85,
              voldata: new $3Dmol.VolumeData(volumedata, "cube"),
              volscheme: new $3Dmol.Gradient.Sinebow(2,0,1)
          },{});
          
      viewer.render();
      });
      viewer.zoomTo();
  });
 */
export class Sinebow extends GradientType {
  gradient = "Sinebow";
  mult: number;
  max: number;
  min: number;
  constructor(min, max) {
    super();
    this.mult = 1.0;
    this.min = min;
    this.max = max;
    if (typeof max == "undefined" && Array.isArray(min) && min.length >= 2) {
      //we were passed a single range
      this.max = min[1];
      this.min = min[0];
    }
    if (max < min) {
      //reverse the order
      this.mult = -1.0;
      this.min *= -1.0;
      this.max *= -1.0;
    }
  };

  //map value to hex color, range is provided
  valueToHex(val, range) {
    var lo, hi;
    val = this.mult * val;
    if (range) {
      lo = range[0];
      hi = range[1];
    } else {
      lo = this.min;
      hi = this.max;
    }

    if (typeof val == "undefined") return 0xffffff;
    var norm = Gradient.normalizeValue(lo, hi, val);
    lo = norm.lo;
    hi = norm.hi;
    val = norm.val;

    var scale = (val - lo) / (hi - lo);
    var h = (5 * scale) / 6.0 + 0.5;
    var r = Math.sin(Math.PI * h);
    r *= r * 255;
    var g = Math.sin(Math.PI * (h + 1 / 3.0));
    g *= g * 255;
    var b = Math.sin(Math.PI * (h + 2 / 3.0));
    b *= b * 255;

    return (
      0x10000 * Math.floor(r) + 0x100 * Math.floor(b) + 0x1 * Math.floor(g)
    );
  };

  //return range used for color mapping, null if none set
  range() {
    if (typeof this.min != "undefined" && typeof this.max != "undefined") {
      return [this.min, this.max] as [number, number];
    }
    return null;
  };
}

//map from names to gradient constructors
export const builtinGradients = {
  rwb: RWB,
  RWB: RWB,
  roygb: ROYGB,
  ROYGB: ROYGB,
  sinebow: Sinebow,
};

export class Gradient extends GradientType {
  static RWB = RWB;
  static ROYGB = ROYGB;
  static Sinebow = Sinebow;
  static builtinGradients = builtinGradients;
  static normalizeValue = normalizeValue;
  static getGradient = getGradient;
  // @ts-ignore
  valueToHex(_value: number, _range: number): number { }
  // @ts-ignore
  range(): [number, number] { }
}