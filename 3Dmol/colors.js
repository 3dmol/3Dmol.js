(function() {
/**
 * Color representation. 
 * @typedef ColorSpec
 * @prop {string} 0xAF10AB - any hex number
 * @prop {string} <html color name>
 */

 /**
 
 * @typedef ColorschemeSpec
 * Built in colorschemes
 *
 * @example //Using a function in order to define the colors. 
  $3Dmol.download("pdb:4UAA",viewer,{},function(){
                  viewer.setBackgroundColor(0xffffffff);
                  var colorAsSnake = function(atom) {
                    return atom.resi % 2 ? 'white': 'green'
                  };

                  viewer.setStyle( {chain:'A'}, { cartoon: {colorfunc: colorAsSnake }});
                  viewer.setStyle( {chain:'B'}, { stick: {colorscheme: 'yellowCarbon'}});

                  viewer.render();
              });
 * @prop {string} <html color>Carbon   - use default element colors but with carbon set to specify html color string
 * @prop {string} ssPyMOL - PyMol secondary colorscheme
 * @prop {string} ssJmol - Jmol secondary colorscheme
 * @prop {string} Jmol - Jmol primary colorscheme
 * @prop {string} default - default colorscheme
 * @prop {string} amino - amino acid colorscheme
 * @prop {string} shapely - shapely protien colorscheme
 * @prop {string} nucleic - nucleic acid colorscheme
 * @prop {string} chain - standard chain colorscheme
 * @prop {string} chainHetatm - chain Hetatm colorscheme
 * @prop {string} prop - atomSpec property. Example 'b'. See AtomSpec.
 * @prop {Gradient} gradient - Allows the user to provide a gradient to the colorscheme.  Is either a $3Dmol.Gradient object or the name of a built-in gradient (rwb, roygb, sinebow)
 * @prop {min} - min value for gradient
 * @prop {max} - max value for gradient
 * @prop {mid} - mid point value for gradient (for rwb)
 * @prop {object} map - map of a certain AtomSpec property to a color of the form `{'prop': 'elem', map:$3Dmol.elementColors.greenCarbon}` Allows the user to provide a mapping of elements to colors to the colorscheme.  This can be done with any properties, and not just 'elem'.
 * @prop {function} colorfunc - Allows the user to provide a function for setting the colorschemes.
 */
 
})();

var htmlColors = $3Dmol.htmlColors = {
    "aliceblue" : 0xF0F8FF,
    "antiquewhite" : 0xFAEBD7,
    "aqua" : 0x00FFFF,
    "aquamarine" : 0x7FFFD4,
    "azure" : 0xF0FFFF,
    "beige" : 0xF5F5DC,
    "bisque" : 0xFFE4C4,
    "black" : 0x000000,
    "blanchedalmond" : 0xFFEBCD,
    "blue" : 0x0000FF,
    "blueviolet" : 0x8A2BE2,
    "brown" : 0xA52A2A,
    "burlywood" : 0xDEB887,
    "cadetblue" : 0x5F9EA0,
    "chartreuse" : 0x7FFF00,
    "chocolate" : 0xD2691E,
    "coral" : 0xFF7F50,
    "cornflowerblue" : 0x6495ED,
    "cornsilk" : 0xFFF8DC,
    "crimson" : 0xDC143C,
    "cyan" : 0x00FFFF,
    "darkblue" : 0x00008B,
    "darkcyan" : 0x008B8B,
    "darkgoldenrod" : 0xB8860B,
    "darkgray" : 0xA9A9A9,
    "darkgrey" : 0xA9A9A9,
    "darkgreen" : 0x006400,
    "darkkhaki" : 0xBDB76B,
    "darkmagenta" : 0x8B008B,
    "darkolivegreen" : 0x556B2F,
    "darkorange" : 0xFF8C00,
    "darkorchid" : 0x9932CC,
    "darkred" : 0x8B0000,
    "darksalmon" : 0xE9967A,
    "darkseagreen" : 0x8FBC8F,
    "darkslateblue" : 0x483D8B,
    "darkslategray" : 0x2F4F4F,
    "darkslategrey" : 0x2F4F4F,
    "darkturquoise" : 0x00CED1,
    "darkviolet" : 0x9400D3,
    "deeppink" : 0xFF1493,
    "deepskyblue" : 0x00BFFF,
    "dimgray" : 0x696969,
    "dimgrey" : 0x696969,
    "dodgerblue" : 0x1E90FF,
    "firebrick" : 0xB22222,
    "floralwhite" : 0xFFFAF0,
    "forestgreen" : 0x228B22,
    "fuchsia" : 0xFF00FF,
    "gainsboro" : 0xDCDCDC,
    "ghostwhite" : 0xF8F8FF,
    "gold" : 0xFFD700,
    "goldenrod" : 0xDAA520,
    "gray" : 0x808080,
    "grey" : 0x808080,
    "green" : 0x008000,
    "greenyellow" : 0xADFF2F,
    "honeydew" : 0xF0FFF0,
    "hotpink" : 0xFF69B4,
    "indianred" : 0xCD5C5C,
    "indigo" : 0x4B0082,
    "ivory" : 0xFFFFF0,
    "khaki" : 0xF0E68C,
    "lavender" : 0xE6E6FA,
    "lavenderblush" : 0xFFF0F5,
    "lawngreen" : 0x7CFC00,
    "lemonchiffon" : 0xFFFACD,
    "lightblue" : 0xADD8E6,
    "lightcoral" : 0xF08080,
    "lightcyan" : 0xE0FFFF,
    "lightgoldenrodyellow" : 0xFAFAD2,
    "lightgray" : 0xD3D3D3,
    "lightgrey" : 0xD3D3D3,
    "lightgreen" : 0x90EE90,
    "lightpink" : 0xFFB6C1,
    "lightsalmon" : 0xFFA07A,
    "lightseagreen" : 0x20B2AA,
    "lightskyblue" : 0x87CEFA,
    "lightslategray" : 0x778899,
    "lightslategrey" : 0x778899,
    "lightsteelblue" : 0xB0C4DE,
    "lightyellow" : 0xFFFFE0,
    "lime" : 0x00FF00,
    "limegreen" : 0x32CD32,
    "linen" : 0xFAF0E6,
    "magenta" : 0xFF00FF,
    "maroon" : 0x800000,
    "mediumaquamarine" : 0x66CDAA,
    "mediumblue" : 0x0000CD,
    "mediumorchid" : 0xBA55D3,
    "mediumpurple" : 0x9370DB,
    "mediumseagreen" : 0x3CB371,
    "mediumslateblue" : 0x7B68EE,
    "mediumspringgreen" : 0x00FA9A,
    "mediumturquoise" : 0x48D1CC,
    "mediumvioletred" : 0xC71585,
    "midnightblue" : 0x191970,
    "mintcream" : 0xF5FFFA,
    "mistyrose" : 0xFFE4E1,
    "moccasin" : 0xFFE4B5,
    "navajowhite" : 0xFFDEAD,
    "navy" : 0x000080,
    "oldlace" : 0xFDF5E6,
    "olive" : 0x808000,
    "olivedrab" : 0x6B8E23,
    "orange" : 0xFFA500,
    "orangered" : 0xFF4500,
    "orchid" : 0xDA70D6,
    "palegoldenrod" : 0xEEE8AA,
    "palegreen" : 0x98FB98,
    "paleturquoise" : 0xAFEEEE,
    "palevioletred" : 0xDB7093,
    "papayawhip" : 0xFFEFD5,
    "peachpuff" : 0xFFDAB9,
    "peru" : 0xCD853F,
    "pink" : 0xFFC0CB,
    "plum" : 0xDDA0DD,
    "powderblue" : 0xB0E0E6,
    "purple" : 0x800080,
    "rebeccapurple" : 0x663399,
    "red" : 0xFF0000,
    "rosybrown" : 0xBC8F8F,
    "royalblue" : 0x4169E1,
    "saddlebrown" : 0x8B4513,
    "salmon" : 0xFA8072,
    "sandybrown" : 0xF4A460,
    "seagreen" : 0x2E8B57,
    "seashell" : 0xFFF5EE,
    "sienna" : 0xA0522D,
    "silver" : 0xC0C0C0,
    "skyblue" : 0x87CEEB,
    "slateblue" : 0x6A5ACD,
    "slategray" : 0x708090,
    "slategrey" : 0x708090,
    "snow" : 0xFFFAFA,
    "springgreen" : 0x00FF7F,
    "steelblue" : 0x4682B4,
    "tan" : 0xD2B48C,
    "teal" : 0x008080,
    "thistle" : 0xD8BFD8,
    "tomato" : 0xFF6347,
    "turquoise" : 0x40E0D0,
    "violet" : 0xEE82EE,
    "wheat" : 0xF5DEB3,
    "white" : 0xFFFFFF,
    "whitesmoke" : 0xF5F5F5,
    "yellow" : 0xFFFF00,
    "yellowgreen" : 0x9ACD32
};
// in an attempt to reduce memory overhead, cache all $3Dmol.Colors
// this makes things a little faster
$3Dmol.CC = {
    rgbRegEx : /rgb(a?)\(\s*([^ ,\)\t]+)\s*,\s*([^ ,\)\t]+)\s*,\s*([^ ,\)\t]+)/i,
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
 
    getHex : function(hex) {
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
            
            let m = this.rgbRegEx.exec(hex);
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
            return htmlColors[hex.toLowerCase()] || 0x000000;
            
        }
        return hex;
    }
    
};


$3Dmol.CC = $3Dmol.CC;
$3Dmol.CC.color = $3Dmol.CC.color;




/** Preset secondary structure color scheme 
 * @struct
 */
$3Dmol.ssColors = $3Dmol.ssColors || {};
//names are in helix-sheet-coil order
$3Dmol.ssColors.pyMol = {'h': 0xff0000, 's':  0xffff00, 'c': 0x00ff00};
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
        'Be': 0xC2FF00,
        'BE': 0xC2FF00,
        'B': 0xFFB5B5,
        'C': 0x909090,
        'N': 0x3050F8,
        'O': 0xFF0D0D,
        'F': 0x90E050,
        'Ne': 0xB3E3F5,
        'NE': 0xB3E3F5,
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
        'Ar': 0x80D1E3,
        'AR': 0x80D1E3,
        'K': 0x8F40D4,
        'Ca': 0x3DFF00,
        'CA': 0x3DFF00,
        'Sc': 0xE6E6E6,
        'SC': 0xE6E6E6,
        'Ti': 0xBFC2C7,
        'TI': 0xBFC2C7,
        'V': 0xA6A6AB,
        'Cr': 0x8A99C7,
        'CR': 0x8A99C7,
        'Mn': 0x9C7AC7,
        'MN': 0x9C7AC7,
        'Fe': 0xE06633,
        'FE': 0xE06633,
        'Co': 0xF090A0,
        'CO': 0xF090A0,
        'Ni': 0x50D050,
        'NI': 0x50D050,
        'Cu': 0xC88033,
        'CU': 0xC88033,
        'Zn': 0x7D80B0,
        'ZN': 0x7D80B0,
        'Ga': 0xC28F8F,
        'GA': 0xC28F8F,
        'Ge': 0x668F8F,
        'GE': 0x668F8F,
        'As': 0xBD80E3,
        'AS': 0xBD80E3,
        'Se': 0xFFA100,
        'SE': 0xFFA100,
        'Br': 0xA62929,
        'BR': 0xA62929,
        'Kr': 0x5CB8D1,
        'KR': 0x5CB8D1,
        'Rb': 0x702EB0,
        'RB': 0x702EB0,
        'Sr': 0x00FF00,
        'SR': 0x00FF00,
        'Y': 0x94FFFF,
        'Zr': 0x94E0E0,
        'ZR': 0x94E0E0,
        'Nb': 0x73C2C9,
        'NB': 0x73C2C9,
        'Mo': 0x54B5B5,
        'MO': 0x54B5B5,
        'Tc': 0x3B9E9E,
        'TC': 0x3B9E9E,
        'Ru': 0x248F8F,
        'RU': 0x248F8F,
        'Rh': 0x0A7D8C,
        'RH': 0x0A7D8C,
        'Pd': 0x006985,
        'PD': 0x006985,
        'Ag': 0xC0C0C0,
        'AG': 0xC0C0C0,
        'Cd': 0xFFD98F,
        'CD': 0xFFD98F,
        'In': 0xA67573,
        'IN': 0xA67573,
        'Sn': 0x668080,
        'SN': 0x668080,
        'Sb': 0x9E63B5,
        'SB': 0x9E63B5,
        'Te': 0xD47A00,
        'TE': 0xD47A00,
        'I': 0x940094,
        'Xe': 0x429EB0,
        'XE': 0x429EB0,
        'Cs': 0x57178F,
        'CS': 0x57178F,
        'Ba': 0x00C900,
        'BA': 0x00C900,
        'La': 0x70D4FF,
        'LA': 0x70D4FF,
        'Ce': 0xFFFFC7,
        'CE': 0xFFFFC7,
        'Pr': 0xD9FFC7,
        'PR': 0xD9FFC7,
        'Nd': 0xC7FFC7,
        'ND': 0xC7FFC7,
        'Pm': 0xA3FFC7,
        'PM': 0xA3FFC7,
        'Sm': 0x8FFFC7,
        'SM': 0x8FFFC7,
        'Eu': 0x61FFC7,
        'EU': 0x61FFC7,
        'Gd': 0x45FFC7,
        'GD': 0x45FFC7,
        'Tb': 0x30FFC7,
        'TB': 0x30FFC7,
        'Dy': 0x1FFFC7,
        'DY': 0x1FFFC7,
        'Ho': 0x00FF9C,
        'HO': 0x00FF9C,
        'Er': 0x00E675,
        'ER': 0x00E675,
        'Tm': 0x00D452,
        'TM': 0x00D452,
        'Yb': 0x00BF38,
        'YB': 0x00BF38,
        'Lu': 0x00AB24,
        'LU': 0x00AB24,
        'Hf': 0x4DC2FF,
        'HF': 0x4DC2FF,
        'Ta': 0x4DA6FF,
        'TA': 0x4DA6FF,
        'W': 0x2194D6,
        'Re': 0x267DAB,
        'RE': 0x267DAB,
        'Os': 0x266696,
        'OS': 0x266696,
        'Ir': 0x175487,
        'IR': 0x175487,
        'Pt': 0xD0D0E0,
        'PT': 0xD0D0E0,
        'Au': 0xFFD123,
        'AU': 0xFFD123,
        'Hg': 0xB8B8D0,
        'HG': 0xB8B8D0,
        'Tl': 0xA6544D,
        'TL': 0xA6544D,
        'Pb': 0x575961,
        'PB': 0x575961,
        'Bi': 0x9E4FB5,
        'BI': 0x9E4FB5,
        'Po': 0xAB5C00,
        'PO': 0xAB5C00,
        'At': 0x754F45,
        'AT': 0x754F45,
        'Rn': 0x428296,
        'RN': 0x428296,
        'Fr': 0x420066,
        'FR': 0x420066,
        'Ra': 0x007D00,
        'RA': 0x007D00,
        'Ac': 0x70ABFA,
        'AC': 0x70ABFA,
        'Th': 0x00BAFF,
        'TH': 0x00BAFF,
        'Pa': 0x00A1FF,
        'PA': 0x00A1FF,
        'U': 0x008FFF,
        'Np': 0x0080FF,
        'NP': 0x0080FF,
        'Pu': 0x006BFF,
        'PU': 0x006BFF,
        'Am': 0x545CF2,
        'AM': 0x545CF2,
        'Cm': 0x785CE3,
        'CM': 0x785CE3,
        'Bk': 0x8A4FE3,
        'BK': 0x8A4FE3,
        'Cf': 0xA136D4,
        'CF': 0xA136D4,
        'Es': 0xB31FD4,
        'ES': 0xB31FD4,
        'Fm': 0xB31FBA,
        'FM': 0xB31FBA,
        'Md': 0xB30DA6,
        'MD': 0xB30DA6,
        'No': 0xBD0D87,
        'NO': 0xBD0D87,
        'Lr': 0xC70066,
        'LR': 0xC70066,
        'Rf': 0xCC0059,
        'RF': 0xCC0059,
        'Db': 0xD1004F,
        'DB': 0xD1004F,
        'Sg': 0xD90045,
        'SG': 0xD90045,
        'Bh': 0xE00038,
        'BH': 0xE00038,
        'Hs': 0xE6002E,
        'HS': 0xE6002E,
        'Mt': 0xEB0026,
        'MT': 0xEB0026
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

$3Dmol.elementColors.greenCarbon = $3Dmol.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.greenCarbon.C = 0x00ff00; //bright green

$3Dmol.elementColors.cyanCarbon =  $3Dmol.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.cyanCarbon.C = 0x00ffff;

$3Dmol.elementColors.magentaCarbon =  $3Dmol.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.magentaCarbon.C = 0xff00ff;

$3Dmol.elementColors.yellowCarbon =  $3Dmol.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.yellowCarbon.C = 0xffff00;

$3Dmol.elementColors.whiteCarbon =  $3Dmol.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.whiteCarbon.C = 0xffffff;

$3Dmol.elementColors.orangeCarbon =  $3Dmol.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.orangeCarbon.C = 0xffa500;

$3Dmol.elementColors.purpleCarbon =  $3Dmol.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.purpleCarbon.C = 0x800080;

$3Dmol.elementColors.blueCarbon =  $3Dmol.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.blueCarbon.C = 0x0000ff;


$3Dmol.residues = {};

/** @property standard amino acid color scheme*/
$3Dmol.residues.amino ={
'ALA' : 0xC8C8C8,        
'ARG' : 0x145AFF,              
'ASN' : 0x00DCDC,              
'ASP' : 0xE60A0A,             
'CYS' : 0xE6E600,             
'GLN' : 0x00DCDC,            
'GLU' : 0xE60A0A,              
'GLY' : 0xEBEBEB,        
'HIS' :  0x8282D2,       
'ILE' :0x0F820F,             
'LEU':0x0F820F,              
'LYS' :0x145AFF,              
'MET' :0xE6E600,             
'PHE':0x3232AA,            
'PRO' :  0xDC9682,    
'SER': 0xFA9600,           
'THR': 0xFA9600,          
'TRP' :   0xB45AB4,    
'TYR' :0x3232AA,            
'VAL' :0x0F820F,            
'ASX'  :0xFF69B4,     
'GLX'  :0xFF69B4,   

};

/** @property shapely amino acid color scheme*/
$3Dmol.residues.shapely ={
'ALA' : 0x8CFF8C,         
'ARG' : 0x00007C,              
'ASN' : 0xFF7C70,              
'ASP' : 0xA00042,             
'CYS' : 0xFFFF70,             
'GLN' : 0xFF4C4C,            
'GLU' : 0x660000,              
'GLY' : 0xFFFFFF,        
'HIS' :  0x7070FF,       
'ILE' :0x004C00,             
'LEU':0x455E45,              
'LYS' :0x4747B8,              
'MET' :0xB8A042,             
'PHE':0x534C52,            
'PRO' :  0x525252,    
'SER': 0xFF7042,           
'THR': 0xB84C00,          
'TRP' :   0x4F4600,    
'TYR' :0x8C704C,            
'VAL' :0xFF8CFF,            
'ASX'  :0xFF00FF,     
'GLX'  :0xFF00FF,    

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
* The user can pass all of these values directly as the colorscheme and they will use the respective colorscheme */
$3Dmol.builtinColorSchemes = {
        'ssPyMol' : {'prop':'ss', map:$3Dmol.ssColors.pyMol},
        'ssJmol' :{'prop':'ss', map:$3Dmol.ssColors.Jmol},
        'Jmol' :{'prop':'elem', map:$3Dmol.elementColors.Jmol},
        'greenCarbon': {'prop': 'elem', map:$3Dmol.elementColors.greenCarbon},
        'default' : {'prop': 'elem', map:$3Dmol.elementColors.defaultColors},
        'amino' : {'prop':'resn', map:$3Dmol.residues.amino},
        'shapely' :{'prop':'resn', map:$3Dmol.residues.shapely},
        'nucleic' :{'prop':'resn', map:$3Dmol.residues.nucleic},
        'chain' :{'prop':'chain', map:$3Dmol.chains.atom},
        'chainHetatm' :{'prop':'chain', map:$3Dmol.chains.hetatm},
        
};

/** Return proper color for atom given style
 * @param {AtomSpec} atom
 * @param {AtomStyle} style
 * @return {$3Dmol.Color}
 */
 
$3Dmol.getColorFromStyle = function(atom, style) {
    var scheme = style.colorscheme;  
    if(typeof($3Dmol.builtinColorSchemes[scheme]) != "undefined") {
        scheme = $3Dmol.builtinColorSchemes[scheme];
    } else if(typeof(scheme) == 'string' && scheme.endsWith('Carbon')) {
        //any color you want of carbon
        var ccolor = scheme.substring(0,scheme.lastIndexOf("Carbon")).toLowerCase();
        if(typeof(htmlColors[ccolor]) != "undefined") {
            var newscheme = $3Dmol.extend({},$3Dmol.elementColors.defaultColors);
            newscheme.C = htmlColors[ccolor];
            $3Dmol.builtinColorSchemes[scheme] = {'prop': 'elem', map:newscheme};
            scheme = $3Dmol.builtinColorSchemes[scheme];
        }        
    }
    
    var color = atom.color;
    if (typeof (style.color) != "undefined" && style.color != "spectrum")
        color = style.color;
    if(typeof(scheme) != "undefined") {
        var prop, val;
        if(typeof($3Dmol.elementColors[scheme]) != "undefined") {
            //name of builtin colorscheme
            scheme = $3Dmol.elementColors[scheme];
            if(typeof(scheme[atom[scheme.prop]]) != "undefined") {
                color = scheme.map[atom[scheme.prop]];
            }
        } else if(typeof(scheme[atom[scheme.prop]]) != 'undefined') {
            //actual color scheme provided
            color = scheme.map[atom[scheme.prop]];
        } else if(typeof(scheme.prop) != 'undefined' &&
                typeof(scheme.gradient) != 'undefined') {         
            //apply a property mapping
            prop = scheme.prop;
            var grad = scheme.gradient; //redefining scheme
            if(typeof($3Dmol.Gradient.builtinGradients[grad]) != "undefined") {
                grad = new $3Dmol.Gradient.builtinGradients[grad](scheme.min, scheme.max, scheme.mid);
            }
            
            var range = grad.range() || [-1,1]; //sensible default
            val = $3Dmol.getAtomProperty(atom, prop);
            if(val != null) {
                color = grad.valueToHex(val, range);
            }
        } else if(typeof(scheme.prop) != 'undefined' &&
                typeof(scheme.map) != 'undefined') {         
            //apply a discrete property mapping
            prop = scheme.prop;
            val = $3Dmol.getAtomProperty(atom, prop);
            if( typeof scheme.map[val] != 'undefined' ) {
                color = scheme.map[val];
            }
        } else if(typeof(style.colorscheme[atom.elem]) != 'undefined') {
            //actual color scheme provided
            color = style.colorscheme[atom.elem];
        } else {
            console.log("Could not interpret colorscheme "+scheme);
        } 
    } 
    else if(typeof(style.colorfunc) != "undefined") {
        //this is a user provided function for turning an atom into a color
        color = style.colorfunc(atom);
    }
    
    var C = $3Dmol.CC.color(color);
    return C;
};
