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
//duplicate ---------------------------------------------------------------------------------------------------------------------------------
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
 * @example viewer.setStyle({chain:'G'},{sphere:{colorscheme:'greenCarbon'}});
 * @typedef ColorschemeSpec
 * @prop {string} greenCarbon   - 0x00FF00
 * @prop {string} cyanCarbon    - 0x00FFFF
 * @prop {string} magentaCarbon - 0xFF00FF
 * @prop {string} yellowCarbon  - 0xFFFF00
 * @prop {string} whiteCarbon   - 0xFFFFFF
 * @prop {string} orangeCarbon  - 0xFF6600
 * @prop {string} purpleCarbon  - 0x100080
 * @prop {string} blueCarbon    - 0x0000FF
 * @prop {string} ssPyMOL - PyMol secondary colorscheme
 * @prop {string} ssJmol - Jmol secondary colorscheme
 * @prop {string} Jmol - Jmol primary colorscheme
 * @prop {string} default - default colorscheme
 * @prop {string} amino - amino acid colorscheme
 * @prop {string} shapely - shapely protien colorscheme
 * @prop {string} nucleic - nucleic acid colorscheme
 * @prop {string} chain - standard chain colorscheme
 * @prop {string} chainHetatm - chain Hetatm colorscheme
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




/** Preset secondary structure color scheme 
 * @struct
 */
$3Dmol.ssColors = $3Dmol.ssColors || {};
//names are in helix-sheet-coil order
$3Dmol.ssColors.pyMOL = {'h': 0xff0000, 's':  0xffff00, 'c': 0x00ff00};
$3Dmol.ssColors.Jmol = {'h': 0xff0080, 's': 0xffc800, 'c': 0xffffff};


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

$3Dmol.residues = {};

/** @property standard amino acid color scheme*/
$3Dmol.residues.amino ={
'Ala' : 0xC8C8C8,        
'Arg' : 0x145AFF,              
'Asn' : 0x00DCDC,              
'Asp' : 0xE60A0A,             
'Cys' : 0xE6E600,             
'Gln' : 0x00DCDC,            
'Glu' : 0xE60A0A,              
'Gly' : 0xEBEBEB,        
'His' :  0x8282D2,       
'Ile' :0x0F820F,             
'Leu':0x0F820F,              
'Lys' :0x145AFF,              
'Met' :0xE6E600,             
'Phe':0x3232AA,            
'Pro' :  0xDC9682,    
'Ser': 0xFA9600,           
'Thr': 0xFA9600,          
'Trp' :   0xB45AB4,    
'Tyr' :0x3232AA,            
'Val' :0x0F820F,            
'Asx'  :0xFF69B4,     
'Glx'  :0xFF69B4,       
'other' :0xBEA06E,

};

/** @property shapely amino acid color scheme*/
$3Dmol.residues.shapely ={
'Ala' : 0x8CFF8C,         
'Arg' : 0x00007C,              
'Asn' : 0xFF7C70,              
'Asp' : 0xA00042,             
'Cys' : 0xFFFF70,             
'Gln' : 0xFF4C4C,            
'Glu' : 0x660000,              
'Gly' : 0xFFFFFF,        
'His' :  0x7070FF,       
'Ile' :0x004C00,             
'Leu':0x455E45,              
'Lys' :0x4747B8,              
'Met' :0xB8A042,             
'Phe':0x534C52,            
'Pro' :  0x525252,    
'Ser': 0xFF7042,           
'Thr': 0xB84C00,          
'Trp' :   0x4F4600,    
'Tyr' :0x8C704C,            
'Val' :0xFF8CFF,            
'Asx'  :0xFF00FF,     
'Glx'  :0xFF00FF,       
'other' :0xFF00FF,

};

/** @property nucleic acid color scheme*/
$3Dmol.residues.nucleic = {
    'A':0xA0A0FF,    
    'G':  0xFF7070,    
    'I':  0x80FFFF,    
    'C':  0xFF8C4B,    
    'T':  0xA0FFA0,    
    'U':   0xFF8080
};
$3Dmol.chains ={};

/** @property chain based standard color scheme */
$3Dmol.chains.atom = {

'A' : 0xC0D0FF,        
'B' : 0xB0FFB0,              
'C' : 0xFFC0C8 ,              
'D' : 0xFFFF80,             
'E' : 0xFFC0FF,             
'F' : 0xB0F0F0,            
'G' : 0xFFD070,              
'H' : 0xF08080,        
'I' :  0xF5DEB3,       
'J' :0x00BFFF,             
'K':0xCD5C5C ,              
'L' :0x66CDAA,              
'M' :0x9ACD32,             
'N':0xEE82EE ,            
'O' :  0x00CED1,    
'P': 0x00FF7F,           
'Q': 0x3CB371 ,          
'R' :   0x00008B,    
'S' :0xBDB76B,            
'T' :0x006400,            
'U'  :0x800000,     
'V'  :0x808000,       
'W' :0x800080,
'X' :0x008080 ,
 'Y':0xB8860B,
 'Z':0xB22222,
};

/** @property hetatm color scheme */
$3Dmol.chains.hetatm = {

'A' : 0x90A0CF,        
'B' : 0x80CF98,              
'C' : 0xCF90B0 ,              
'D' : 0xCFCF70,             
'E' : 0xCF90CF,             
'F' : 0x80C0C0,            
'G' : 0xCFA060,              
'H' : 0xC05070,        
'I' : 0xC5AE83,       
'J' :0x00A7CF,             
'K':0xB54C4C ,              
'L' :0x56B592,              
'M' :0x8AB52A,             
'N':0xBE72BE ,            
'O' :  0x00B6A1,    
'P': 0x00CF6F,           
'Q': 0x349B61 ,          
'R' :   0x0000BB  ,    
'S' :0xA59F5B,            
'T' :0x009400,            
'U'  :0xB00000,     
'V'  :0xB0B000,       
'W' :0xB000B0,
'X' :0x00B0B0 ,
 'Y':0xE8B613,
 'Z':0xC23232,
};

/** @property built in color schemes 
* The user can pass all of these values directly as the colorscheme and they will use the respective colorscheme :
        'ssPyMOL' : {'prop':'ss', map:$3Dmol.ssColors.pyMOL},
        'ssJmol' :{'prop':'ss', map:$3Dmol.ssColors.Jmol},
        'Jmol' :{'prop':'elem', map:$3Dmol.elementColors.Jmol},
        'default' : {'prop': 'elem', map:$3Dmol.elementColors.defaultColors},
        'greenCarbon' : {'prop':'elem', map:$3Dmol.elementColors.greenCarbon},
        'cyanCarbon' : {'prop':'elem', map:$3Dmol.elementColors.cyanCarbon},
        'magentaCarbon' : {'prop':'elem', map:$3Dmol.elementColors.magentaCarbon},
        'yellowCarbon' : {'prop':'elem', map:$3Dmol.elementColors.yellowCarbon},
        'whiteCarbon' : {'prop':'elem', map:$3Dmol.elementColors.whiteCarbon},
        'orangeCarbon' : {'prop':'elem', map:$3Dmol.elementColors.orangeCarbon},
        'purpleCarbon' : {'prop':'elem', map:$3Dmol.elementColors.purpleCarbon},
        'blueCarbon' : {'prop':'elem', map:$3Dmol.elementColors.blueCarbon},
        'amino' : {'prop':'resAmino', map:$3Dmol.residues.amino},
        'shapely' :{'prop':'resShapely', map:$3Dmol.residues.shapely},
        'nucleic' :{'prop':'resNucleic', map:$3Dmol.residues.nucleic},
        'chain' :{'prop':'chain', map:$3Dmol.chains.atom},
        'chainHetatm' :{'prop':'chain', map:$3Dmol.chains.hetatm},

*/
$3Dmol.builtinColorSchemes = {
        'ssPyMOL' : {'prop':'ss', map:$3Dmol.ssColors.pyMOL},
        'ssJmol' :{'prop':'ss', map:$3Dmol.ssColors.Jmol},
        'Jmol' :{'prop':'elem', map:$3Dmol.elementColors.Jmol},
        'default' : {'prop': 'elem', map:$3Dmol.elementColors.defaultColors},
        'greenCarbon' : {'prop':'elem', map:$3Dmol.elementColors.greenCarbon},
        'cyanCarbon' : {'prop':'elem', map:$3Dmol.elementColors.cyanCarbon},
        'magentaCarbon' : {'prop':'elem', map:$3Dmol.elementColors.magentaCarbon},
        'yellowCarbon' : {'prop':'elem', map:$3Dmol.elementColors.yellowCarbon},
        'whiteCarbon' : {'prop':'elem', map:$3Dmol.elementColors.whiteCarbon},
        'orangeCarbon' : {'prop':'elem', map:$3Dmol.elementColors.orangeCarbon},
        'purpleCarbon' : {'prop':'elem', map:$3Dmol.elementColors.purpleCarbon},
        'blueCarbon' : {'prop':'elem', map:$3Dmol.elementColors.blueCarbon},
        'amino' : {'prop':'resAmino', map:$3Dmol.residues.amino},
        'shapely' :{'prop':'resShapely', map:$3Dmol.residues.shapely},
        'nucleic' :{'prop':'resNucleic', map:$3Dmol.residues.nucleic},
        'chain' :{'prop':'chain', map:$3Dmol.chains.atom},
        'chainHetatm' :{'prop':'chain', map:$3Dmol.chains.hetatm},
        
};

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
        if(typeof($3Dmol.builtinColorSchemes[style.colorscheme]) != "undefined") {
            //name of builtin colorscheme
            var scheme = $3Dmol.builtinColorSchemes[style.colorscheme].map;
            if(typeof(scheme[atom.elem]) != "undefined") {
                color = scheme[atom.elem];
            }
        } 
        else if(typeof($3Dmol.elementColors[style.colorscheme]) != "undefined") {
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
        } else if(typeof(style.colorscheme.prop) != 'undefined' &&
                typeof(style.colorscheme.map) != 'undefined') {         
            //apply a discrete property mapping
            var prop = style.colorscheme.prop;
            var val = $3Dmol.getAtomProperty(atom, prop);
            if( typeof style.colorscheme.map[val] != 'undefined' ) {
                color = style.colorscheme.map[val];
            }
        }
    } 
    else if(typeof(style.colorfunc) != "undefined") {
        //this is a user provided function for turning an atom into a color
        color = style.colorfunc(atom);
    }
    
    var C = $3Dmol.CC.color(color);
    return C;
};
