//this is only used for create the enum documentation in JSDoc
(function() {
/**
 * Color representation. 
 * @typedef ColorSpec
 * @prop {string} 0xAF10AB - any hex number
 * @prop {string} white   - 0xFFFFFF
 * @prop {string} silver  - 0xC0C0C0
 * @prop {string} gray    - 0x808080
 * @prop {string} grey    - 0x808080
 * @prop {string} black   - 0x000000
 * @prop {string} red     - 0xFF0000
 * @prop {string} maroon  - 0x800000
 * @prop {string} yellow  - 0xFFFF00
 * @prop {string} orange  - 0xFF6600
 * @prop {string} olive   - 0x808000
 * @prop {string} lime    - 0x00FF00
 * @prop {string} green   - 0x008000
 * @prop {string} aqua    - 0x00FFFF
 * @prop {string} cyan    - 0x00FFFF
 * @prop {string} teal    - 0x008080
 * @prop {string} blue    - 0x0000FF
 * @prop {string} navy    - 0x000080
 * @prop {string} fuchsia - 0xFF00FF
 * @prop {string} magenta - 0xFF00FF
 * @prop {string} purple  - 0x800080
 */

$3Dmol.elementColors.greenCarbon['C'] = 0x00ff00;


$3Dmol.elementColors.cyanCarbon['C'] = 0x00ffff;


$3Dmol.elementColors.magentaCarbon['C'] = 0xff00ff;


$3Dmol.elementColors.yellowCarbon['C'] = 0xffff00;


$3Dmol.elementColors.whiteCarbon['C'] = 0xffffff;


$3Dmol.elementColors.orangeCarbon['C'] = 0xff6600;


$3Dmol.elementColors.purpleCarbon['C'] = 0x800080;

$3Dmol.elementColors.blueCarbon['C'] = 0x0000ff;

 /**
 * Color scheme representation. 
 * @typedef ColorschemeSpec
 * @prop {string} greenCarbon   - 0x00FF00
 * @prop {string} cyanCarbon    - 0x00FFFF
 * @prop {string} magentaCarbon - 0xFF00FF
 * @prop {string} yellowCarbon  - 0xFFFF00
 * @prop {string} whiteCarbon   - 0xFFFFFF
 * @prop {string} orangeCarbon  - 0xFF6600
 * @prop {string} purpleCarbon  - 0x100080
 * @prop {string} blueCarbon    - 0x0000FF
 */
 
});

// in an attempt to reduce memory overhead, cache all $3Dmol.Colors
// this makes things a little faster
$3Dmol.CC = {
    cache : {0:new $3Dmol.Color(0)},
    color : function color_(hex) {
        // Undefined values default to black
        if(!hex)
            return this.cache[0];
        // cache hits
        if(typeof(this.cache[hex]) !== "undefined") {
            return this.cache[hex];
        }
        // arrays
        else if(hex && hex.constructor === Array) {
            // parse elements recursively
            return hex.map(color_,this);
        }
        // numbers and hex strings
        hex = this.getHex(hex);
        if(typeof hex === 'number') {
            var c = new $3Dmol.Color(hex);
            this.cache[hex] = c;
            return c;
        } else {
            // pass through $3Dmol.Color & other objects
            return hex;
        }
    },
 
    colorTab : {
        'white' : 0xFFFFFF,
        'silver' : 0xC0C0C0,
        'gray' : 0x808080,
        'grey' : 0x808080,
        'black' : 0x000000,
        'red' : 0xFF0000,
        'maroon' : 0x800000,
        'yellow' : 0xFFFF00,
        'orange' : 0xFF6600,
        'olive' : 0x808000,
        'lime' : 0x00FF00,
        'green' : 0x008000,
        'aqua' : 0x00FFFF,
        'cyan' : 0x00FFFF,
        'teal' : 0x008080,
        'blue' : 0x0000FF,
        'navy' : 0x000080,
        'fuchsia' : 0xFF00FF,
        'magenta' : 0xFF00FF,
        'purple' : 0x800080
    },    
    getHex : function(hex) {
        if (!isNaN(parseInt(hex)))
            return parseInt(hex);
        
        else if (typeof(hex) === 'string') {
            
            return this.colorTab[hex.trim().toLowerCase()] || 0x000000;
        }
        return hex;
    }
    
};


$3Dmol['CC'] = $3Dmol.CC;
$3Dmol['CC']['color'] = $3Dmol.CC.color;



/** Return proper color for atom given style
 * @param {AtomSpec} atom
 * @param {AtomStyle} style
 * @return {$3Dmol.Color}
 */
$3Dmol.getColorFromStyle = function(atom, style) {
    var color = atom.color;
    if (typeof (style.color) != "undefined" && style.color != "spectrum")
        color = style.color;
    if(typeof(style.colorscheme) != "undefined") {
        if(typeof($3Dmol.elementColors[style.colorscheme]) != "undefined") {
            //name of builtin colorscheme
            var scheme = $3Dmol.elementColors[style.colorscheme];
            if(typeof(scheme[atom.elem]) != "undefined") {
                color = scheme[atom.elem];
            }
        } else if(typeof(style.colorscheme[atom.elem]) != 'undefined') {
            //actual color scheme provided
            color = style.colorscheme[atom.elem];
        } else if(typeof(style.colorscheme.prop) != 'undefined' &&
                typeof(style.colorscheme.gradient) != 'undefined') {         
            //apply a property mapping
            var prop = style.colorscheme.prop;
            var scheme = style.colorscheme.gradient;
            var range = scheme.range() || [-1,1]; //sensible default
            var val = $3Dmol.getAtomProperty(atom, prop);
            if(val != null) {
                color = scheme.valueToHex(val, range);
            }
        }
    } 
    else if(typeof(style.colorfunc) != "undefined") {
    	//this is a user provided function for turning an atom into a color
    	color = style.colorfunc(atom);
    }
    
    var C = $3Dmol.CC.color(color);
    return C;
}

/** Preset element coloring - from individual element colors to entire mappings (e.g. '$3Dmol.elementColors.Jmol' colors atoms with Jmol stylings)
 * @struct
 */
$3Dmol.elementColors = $3Dmol.elementColors || {};

$3Dmol.elementColors.defaultColor = 0xff1493;

/** @property Jmol-like element colors*/
$3Dmol.elementColors.Jmol = {
        'H': 0xFFFFFF,
        'He': 0xD9FFFF,
        'HE': 0xD9FFFF,
        'Li': 0xCC80FF,
        'LI': 0xCC80FF,
        'B': 0xFFB5B5,
        'C': 0x909090,
        'N': 0x3050F8,
        'O': 0xFF0D0D,
        'F': 0x90E050,
        'Na': 0xAB5CF2,
        'NA': 0xAB5CF2,
        'Mg': 0x8AFF00,
        'MG': 0x8AFF00,
        'Al': 0xBFA6A6,
        'AL': 0xBFA6A6,
        'Si': 0xF0C8A0,
        'SI': 0xF0C8A0,
        'P': 0xFF8000,
        'S': 0xFFFF30,
        'Cl': 0x1FF01F,
        'CL': 0x1FF01F,
        'Ca': 0x3DFF00,
        'CA': 0x3DFF00,
        'Ti': 0xBFC2C7,
        'TI': 0xBFC2C7,
        'Cr': 0x8A99C7,
        'CR': 0x8A99C7,
        'Mn': 0x9C7AC7,
        'MN': 0x9C7AC7,
        'Fe': 0xE06633,
        'FE': 0xE06633,
        'Ni': 0x50D050,
        'NI': 0x50D050,
        'Cu': 0xC88033,
        'CU': 0xC88033,
        'Zn': 0x7D80B0,
        'ZN': 0x7D80B0,
        'Br': 0xA62929,
        'BR': 0xA62929,
        'Ag': 0xC0C0C0,
        'AG': 0xC0C0C0,
        'I': 0x940094,
        'Ba': 0x00C900,
        'BA': 0x00C900,
        'Au': 0xFFD123,
        'AU': 0xFFD123
};

/** @property rasmol-like element colors */
$3Dmol.elementColors.rasmol = {
        'H': 0xFFFFFF,
        'He': 0xFFC0CB,
        'HE': 0xFFC0CB,
        'Li': 0xB22222,
        'LI': 0xB22222,
        'B': 0x00FF00,
        'C': 0xC8C8C8,
        'N': 0x8F8FFF,
        'O': 0xF00000,
        'F': 0xDAA520,
        'Na': 0x0000FF,
        'NA': 0x0000FF,
        'Mg': 0x228B22,
        'MG': 0x228B22,
        'Al': 0x808090,
        'AL': 0x808090,
        'Si': 0xDAA520,
        'SI': 0xDAA520,
        'P': 0xFFA500,
        'S': 0xFFC832,
        'Cl': 0x00FF00,
        'CL': 0x00FF00,
        'Ca': 0x808090,
        'CA': 0x808090,
        'Ti': 0x808090,
        'TI': 0x808090,
        'Cr': 0x808090,
        'CR': 0x808090,
        'Mn': 0x808090,
        'MN': 0x808090,
        'Fe': 0xFFA500,
        'FE': 0xFFA500,
        'Ni': 0xA52A2A,
        'NI': 0xA52A2A,
        'Cu': 0xA52A2A,
        'CU': 0xA52A2A,
        'Zn': 0xA52A2A,
        'ZN': 0xA52A2A,
        'Br': 0xA52A2A,
        'BR': 0xA52A2A,
        'Ag': 0x808090,
        'AG': 0x808090,
        'I': 0xA020F0,
        'Ba': 0xFFA500,
        'BA': 0xFFA500,
        'Au': 0xDAA520,
        'AU': 0xDAA520    
};

$3Dmol.elementColors.defaultColors = $3Dmol.elementColors.rasmol;

$3Dmol.elementColors.greenCarbon = $.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.greenCarbon['C'] = 0x00ff00;

$3Dmol.elementColors.cyanCarbon =  $.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.cyanCarbon['C'] = 0x00ffff;

$3Dmol.elementColors.magentaCarbon =  $.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.magentaCarbon['C'] = 0xff00ff;

$3Dmol.elementColors.yellowCarbon =  $.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.yellowCarbon['C'] = 0xffff00;

$3Dmol.elementColors.whiteCarbon =  $.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.whiteCarbon['C'] = 0xffffff;

$3Dmol.elementColors.orangeCarbon =  $.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.orangeCarbon['C'] = 0xff6600;

$3Dmol.elementColors.purpleCarbon =  $.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.purpleCarbon['C'] = 0x800080;

$3Dmol.elementColors.blueCarbon =  $.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.blueCarbon['C'] = 0x0000ff;
