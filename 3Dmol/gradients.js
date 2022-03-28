// color scheme mappings

/** Color mapping gradients
 * @interface
 * @param {number} min
 * @param {number} max
 */
$3Dmol.Gradient = function Gradient(/* min, max */) {};

/**
 * Map value to hex color
 * @param {number} val
 * @param {number} range
 * @returns {number}
 */
$3Dmol.Gradient.valueToHex = function valueToHex(/* val, range */) {};
// return range used for color mapping, null if none set
$3Dmol.Gradient.range = function range() {};

// if lo > hi, flip, also cap
$3Dmol.Gradient.normalizeValue = function normalizeValue(lo, hi, val) {
  if (hi >= lo) {
    if (val < lo) val = lo;
    if (val > hi) val = hi;
    return {lo, hi, val};
  }

  if (val > lo) val = lo;
  if (val < hi) val = hi;
  // flip the meaning of val, lo, hi
  val = lo - val + hi;
  return {lo: hi, hi: lo, val};
};

// return a Gradient object, even if what is specified is descriptive
$3Dmol.Gradient.getGradient = function getGradient(grad) {
  if (grad instanceof $3Dmol.Gradient) {
    return grad;
  }
  if (grad.gradient !== undefined && $3Dmol.Gradient.builtinGradients[grad.gradient]) {
    const min = grad.min === undefined ? -1 : grad.min;
    const max = grad.max === undefined ? 1 : grad.max;
    if (grad.mid === undefined) {
      return new $3Dmol.Gradient.builtinGradients[grad.gradient](min, max);
    }
    return new $3Dmol.Gradient.builtinGradients[grad.gradient](min, max, grad.mid);
  }
  return grad;
};

/**
 * Color scheme red to white to blue, for charges
 * Reverse gradients are supported when min>max so that the colors are displayed in reverse order.
 * @constructor
 * @implements {$3Dmol.Gradient}
 */
$3Dmol.Gradient.RWB = function RWB(min, max, mid) {
  const mult = 1.0;
  if (typeof max == 'undefined' && Array.isArray(min) && min.length >= 2) {
    // we were passed a single range
    max = min[1];
    min = min[0];
  }
  // map value to hex color, range is provided
  this.valueToHex = function (val, range) {
    let lo;
    let hi;
    val *= mult; // reverse if necessary
    if (range) {
      lo = range[0];
      hi = range[1];
    } else {
      lo = min;
      hi = max;
    }

    if (val === undefined) return 0xffffff;

    const norm = $3Dmol.Gradient.normalizeValue(lo, hi, val);
    lo = norm.lo;
    hi = norm.hi;
    val = norm.val;

    let middle = (hi + lo) / 2;
    if (range && typeof range[2] != 'undefined') middle = range[2];
    else if (typeof mid != 'undefined') middle = mid; // allow user to specify midpoint
    else middle = (lo + hi) / 2;
    let scale;
    let color;

    // scale bottom from red to white
    if (val < middle) {
      scale = Math.floor(255 * Math.sqrt((val - lo) / (middle - lo)));
      color = 0xff0000 + 0x100 * scale + scale;
      return color;
    }
    if (val > middle) {
      // form white to blue
      scale = Math.floor(255 * Math.sqrt(1 - (val - middle) / (hi - middle)));
      color = 0x10000 * scale + 0x100 * scale + 0xff;
      return color;
    } // val == middle
    return 0xffffff;
  };

  // return range used for color mapping, null if none set
  this.range = function () {
    if (typeof min != 'undefined' && typeof max != 'undefined') {
      return [min, max];
    }
    return null;
  };
};

/**
 * rainbow gradient, but without purple to match jmol
 * Reverse gradients are supported when min>max so that the colors are displayed in reverse order.
 * @constructor
 * @implements {$3Dmol.Gradient}
 */
$3Dmol.Gradient.ROYGB = function ROYGB(min, max) {
  const mult = 1.0;
  if (typeof max == 'undefined' && Array.isArray(min) && min.length >= 2) {
    // we were passed a single range
    max = min[1];
    min = min[0];
  }

  // map value to hex color, range is provided
  this.valueToHex = function (val, range) {
    let lo;
    let hi;
    val *= mult;
    if (range) {
      lo = range[0];
      hi = range[1];
    } else {
      lo = min;
      hi = max;
    }

    if (typeof val == 'undefined') return 0xffffff;

    const norm = $3Dmol.Gradient.normalizeValue(lo, hi, val);
    lo = norm.lo;
    hi = norm.hi;
    val = norm.val;

    const mid = (lo + hi) / 2;
    const q1 = (lo + mid) / 2;
    const q3 = (mid + hi) / 2;

    let scale;
    let color;
    if (val < q1) {
      // scale green up, red up, blue down
      scale = Math.floor(255 * Math.sqrt((val - lo) / (q1 - lo)));
      color = 0xff0000 + 0x100 * scale + 0;
      return color;
    }
    if (val < mid) {
      // scale red down, green up, blue down
      scale = Math.floor(255 * Math.sqrt(1 - (val - q1) / (mid - q1)));
      color = 0x010000 * scale + 0xff00 + 0x0;
      return color;
    }
    if (val < q3) {
      // scale blue up, red down, green up
      scale = Math.floor(255 * Math.sqrt((val - mid) / (q3 - mid)));
      color = 0x000000 + 0xff00 + 0x1 * scale;
      return color;
    }
    // scale green down, blue up, red down
    scale = Math.floor(255 * Math.sqrt(1 - (val - q3) / (hi - q3)));
    color = 0x000000 + 0x0100 * scale + 0xff;
    return color;
  };

  // return range used for color mapping, null if none set
  this.range = function () {
    if (typeof min != 'undefined' && typeof max != 'undefined') {
      return [min, max];
    }
    return null;
  };
};

/**
 * rainbow gradient with constant saturation, all the way to purple!
 * Reverse gradients are supported when min>max so that the colors are displayed in reverse order.
 * @constructor
 * @implements {$3Dmol.Gradient}
 */
$3Dmol.Gradient.Sinebow = function (min, max) {
  let mult = 1.0;
  if (typeof max == 'undefined' && Array.isArray(min) && min.length >= 2) {
    // we were passed a single range
    max = min[1];
    min = min[0];
  }
  if (max < min) {
    // reverse the order
    mult = -1.0;
    min *= -1.0;
    max *= -1.0;
  }
  // map value to hex color, range is provided
  this.valueToHex = function (val, range) {
    let lo;
    let hi;
    val *= mult;
    if (range) {
      lo = range[0];
      hi = range[1];
    } else {
      lo = min;
      hi = max;
    }

    if (typeof val == 'undefined') return 0xffffff;
    const norm = $3Dmol.Gradient.normalizeValue(lo, hi, val);
    lo = norm.lo;
    hi = norm.hi;
    val = norm.val;

    const scale = (val - lo) / (hi - lo);
    const h = (5 * scale) / 6.0 + 0.5;
    let r = Math.sin(Math.PI * h);
    r *= r * 255;
    let g = Math.sin(Math.PI * (h + 1 / 3.0));
    g *= g * 255;
    let b = Math.sin(Math.PI * (h + 2 / 3.0));
    b *= b * 255;

    return 0x10000 * Math.floor(r) + 0x100 * Math.floor(b) + 0x1 * Math.floor(g);
  };

  // return range used for color mapping, null if none set
  this.range = function () {
    if (typeof min != 'undefined' && typeof max != 'undefined') {
      return [min, max];
    }
    return null;
  };
};

// map from names to gradient constructors
$3Dmol.Gradient.builtinGradients = {
  rwb: $3Dmol.Gradient.RWB,
  RWB: $3Dmol.Gradient.RWB,
  roygb: $3Dmol.Gradient.ROYGB,
  ROYGB: $3Dmol.Gradient.ROYGB,
  sinebow: $3Dmol.Gradient.Sinebow,
};
