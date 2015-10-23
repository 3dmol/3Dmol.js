//color scheme mappings
var $3Dmol = $3Dmol || {};

/** Color mapping gradiens
 * @interface
 * @param {number} min
 * @param {number} max
 */
$3Dmol.Gradient = function(min, max) {};

/**
 * Map value to hex color
 * @param {number} val
 * @param {number} range
 * @returns {number}
 */
$3Dmol.Gradient.valueToHex = function(val, range) {};
//return range used for color mapping, null if none set
$3Dmol.Gradient.range = function() {};



/**
 * Color scheme red to white to blue, for charges
 * @constructor
 * @implements {$3Dmol.Gradient}
 */
$3Dmol.Gradient.RWB = function(min, max,mid) {
    var mult = 1.0;
    if(typeof(max) == 'undefined' && $.isArray(min) && min.length >= 2) {
        //we were passed a single range
        max = min[1];
        min = min[0];
    }
    if(max < min) { //reverse the order
        mult = -1.0;
        min *= -1.0;
        max *= -1.0;
    }
        
    //map value to hex color, range is provided
    this.valueToHex = function(val, range) {
        var lo, hi;
        val = mult*val; //reverse if necessary
        if(range) {
            lo = range[0];
            hi = range[1];
        }
        else {
            lo = min;
            hi = max;
        }
    
        if(val === undefined)
            return 0xffffff;
        
        if(val < lo) val = lo;
        if(val > hi) val = hi;
        
        var middle = (hi+lo)/2;
        if(typeof(mid) != 'undefined')
            middle = mid; //allow user to specify midpoint
        var scale, color;
        
        //scale bottom from red to white
        if(val <= middle) {
            scale = Math.floor(255*Math.sqrt((val-lo)/(middle-lo)));
            color = 0xff0000 + 0x100*scale + scale;
            return color;
        }
        else { //form white to blue
            scale = Math.floor(255*Math.sqrt((1-(val-middle)/(hi-middle))));
            color =  0x10000*scale+0x100*scale+0xff;
            return color;
        }
    };
    

    //return range used for color mapping, null if none set
    this.range = function() {
        if(typeof(min) != "undefined" && typeof(max) != "undefined") {
            return [min,max];
        }
        return null;
    };

};

/**
 * rainbow gradient, but without purple to match jmol
 * @constructor
 * @implements {$3Dmol.Gradient}
 */
$3Dmol.Gradient.ROYGB = function(min, max) {
    var mult = 1.0;
    if(typeof(max) == 'undefined' && $.isArray(min) && min.length >= 2) {
        //we were passed a single range
        max = min[1];
        min = min[0];
    }
    if(max < min) { //reverse the order
        mult = -1.0;
        min *= -1.0;
        max *= -1.0;
    }
    
    //map value to hex color, range is provided
    this.valueToHex = function(val, range) {
        var lo, hi;
        val = mult*val;
        if(range) {
            lo = range[0];
            hi = range[1];
        }
        else {
            lo = min;
            hi = max;
        }
    
        if(typeof(val) == "undefined")
            return 0xffffff;
        
        if(val < lo) val = lo;
        if(val > hi) val = hi;
        
        var mid = (lo+hi)/2;
        var q1 = (lo+mid)/2;
        var q3 = (mid+hi)/2;
        
        var scale, color;
        
        if(val < q1) { //scale green up, red up, blue down
            scale = Math.floor(255*Math.sqrt((val-lo)/(q1-lo)));
            color = 0xff0000 + 0x100*scale + 0;
            return color;
        }
        else if(val < mid) { //scale red down, green up, blue down
            scale = Math.floor(255*Math.sqrt((1-(val-q1)/(mid-q1))));
            color =  0x010000*scale+0xff00+0x0;
            return color;
        }
        else if(val < q3) { //scale blue up, red down, green up
            scale = Math.floor(255*Math.sqrt((val-mid)/(q3-mid)));
            color = 0x000000 + 0xff00 + 0x1*scale;
            return color;
        }
        else { //scale green down, blue up, red down
            scale = Math.floor(255*Math.sqrt((1-(val-q3)/(hi-q3))));
            color =  0x000000+0x0100*scale+0xff;
            return color;
        }        
    };
   

    //return range used for color mapping, null if none set
    this.range = function() {
        if(typeof(min) != "undefined" && typeof(max) != "undefined") {
            return [min,max];
        }
        return null;
    };

};

/**
 * rainbow gradient with constant saturation, all the way to purple!
 * @constructor
 * @implements {$3Dmol.Gradient}
 */
$3Dmol.Gradient.Sinebow = function(min, max) {
    var mult = 1.0;
    if(typeof(max) == 'undefined' && $.isArray(min) && min.length >= 2) {
        //we were passed a single range
        max = min[1];
        min = min[0];
    }
    if(max < min) { //reverse the order
        mult = -1.0;
        min *= -1.0;
        max *= -1.0;
    }
    //map value to hex color, range is provided
    this.valueToHex = function(val, range) {
        var lo, hi;
        val = mult*val;
        if(range) {
            lo = range[0];
            hi = range[1];
        }
        else {
            lo = min;
            hi = max;
        }
    
        if(typeof(val) == "undefined")
            return 0xffffff;
        
        if(val < lo) val = lo;
        if(val > hi) val = hi;
        
        var scale = (val-lo)/(hi-lo);
        var h = (5*scale/6.0+0.5);
        var r = Math.sin(Math.PI*h);
        r *= r*255;
        var g = Math.sin(Math.PI*(h+1/3.0));
        g *= g*255;
        var b = Math.sin(Math.PI*(h+2/3.0));
        b *= b*255;
        
        return 0x10000*Math.floor(r)+0x100*Math.floor(b)+0x1*Math.floor(g);
    };
    

    //return range used for color mapping, null if none set
    this.range = function() {
        if(typeof(min) != "undefined" && typeof(max) != "undefined") {
            return [min,max];
        }
        return null;
    };

};
