import { builtinGradients } from "./Gradient";
import { getAtomProperty } from "./utils/getAtomProperty";
import { CC } from "./WebGL";

export const htmlColors = {
  aliceblue: 0xf0f8ff,
  antiquewhite: 0xfaebd7,
  aqua: 0x00ffff,
  aquamarine: 0x7fffd4,
  azure: 0xf0ffff,
  beige: 0xf5f5dc,
  bisque: 0xffe4c4,
  black: 0x000000,
  blanchedalmond: 0xffebcd,
  blue: 0x0000ff,
  blueviolet: 0x8a2be2,
  brown: 0xa52a2a,
  burlywood: 0xdeb887,
  cadetblue: 0x5f9ea0,
  chartreuse: 0x7fff00,
  chocolate: 0xd2691e,
  coral: 0xff7f50,
  cornflowerblue: 0x6495ed,
  cornsilk: 0xfff8dc,
  crimson: 0xdc143c,
  cyan: 0x00ffff,
  darkblue: 0x00008b,
  darkcyan: 0x008b8b,
  darkgoldenrod: 0xb8860b,
  darkgray: 0xa9a9a9,
  darkgrey: 0xa9a9a9,
  darkgreen: 0x006400,
  darkkhaki: 0xbdb76b,
  darkmagenta: 0x8b008b,
  darkolivegreen: 0x556b2f,
  darkorange: 0xff8c00,
  darkorchid: 0x9932cc,
  darkred: 0x8b0000,
  darksalmon: 0xe9967a,
  darkseagreen: 0x8fbc8f,
  darkslateblue: 0x483d8b,
  darkslategray: 0x2f4f4f,
  darkslategrey: 0x2f4f4f,
  darkturquoise: 0x00ced1,
  darkviolet: 0x9400d3,
  deeppink: 0xff1493,
  deepskyblue: 0x00bfff,
  dimgray: 0x696969,
  dimgrey: 0x696969,
  dodgerblue: 0x1e90ff,
  firebrick: 0xb22222,
  floralwhite: 0xfffaf0,
  forestgreen: 0x228b22,
  fuchsia: 0xff00ff,
  gainsboro: 0xdcdcdc,
  ghostwhite: 0xf8f8ff,
  gold: 0xffd700,
  goldenrod: 0xdaa520,
  gray: 0x808080,
  grey: 0x808080,
  green: 0x008000,
  greenyellow: 0xadff2f,
  honeydew: 0xf0fff0,
  hotpink: 0xff69b4,
  indianred: 0xcd5c5c,
  indigo: 0x4b0082,
  ivory: 0xfffff0,
  khaki: 0xf0e68c,
  lavender: 0xe6e6fa,
  lavenderblush: 0xfff0f5,
  lawngreen: 0x7cfc00,
  lemonchiffon: 0xfffacd,
  lightblue: 0xadd8e6,
  lightcoral: 0xf08080,
  lightcyan: 0xe0ffff,
  lightgoldenrodyellow: 0xfafad2,
  lightgray: 0xd3d3d3,
  lightgrey: 0xd3d3d3,
  lightgreen: 0x90ee90,
  lightpink: 0xffb6c1,
  lightsalmon: 0xffa07a,
  lightseagreen: 0x20b2aa,
  lightskyblue: 0x87cefa,
  lightslategray: 0x778899,
  lightslategrey: 0x778899,
  lightsteelblue: 0xb0c4de,
  lightyellow: 0xffffe0,
  lime: 0x00ff00,
  limegreen: 0x32cd32,
  linen: 0xfaf0e6,
  magenta: 0xff00ff,
  maroon: 0x800000,
  mediumaquamarine: 0x66cdaa,
  mediumblue: 0x0000cd,
  mediumorchid: 0xba55d3,
  mediumpurple: 0x9370db,
  mediumseagreen: 0x3cb371,
  mediumslateblue: 0x7b68ee,
  mediumspringgreen: 0x00fa9a,
  mediumturquoise: 0x48d1cc,
  mediumvioletred: 0xc71585,
  midnightblue: 0x191970,
  mintcream: 0xf5fffa,
  mistyrose: 0xffe4e1,
  moccasin: 0xffe4b5,
  navajowhite: 0xffdead,
  navy: 0x000080,
  oldlace: 0xfdf5e6,
  olive: 0x808000,
  olivedrab: 0x6b8e23,
  orange: 0xffa500,
  orangered: 0xff4500,
  orchid: 0xda70d6,
  palegoldenrod: 0xeee8aa,
  palegreen: 0x98fb98,
  paleturquoise: 0xafeeee,
  palevioletred: 0xdb7093,
  papayawhip: 0xffefd5,
  peachpuff: 0xffdab9,
  peru: 0xcd853f,
  pink: 0xffc0cb,
  plum: 0xdda0dd,
  powderblue: 0xb0e0e6,
  purple: 0x800080,
  rebeccapurple: 0x663399,
  red: 0xff0000,
  rosybrown: 0xbc8f8f,
  royalblue: 0x4169e1,
  saddlebrown: 0x8b4513,
  salmon: 0xfa8072,
  sandybrown: 0xf4a460,
  seagreen: 0x2e8b57,
  seashell: 0xfff5ee,
  sienna: 0xa0522d,
  silver: 0xc0c0c0,
  skyblue: 0x87ceeb,
  slateblue: 0x6a5acd,
  slategray: 0x708090,
  slategrey: 0x708090,
  snow: 0xfffafa,
  springgreen: 0x00ff7f,
  steelblue: 0x4682b4,
  tan: 0xd2b48c,
  teal: 0x008080,
  thistle: 0xd8bfd8,
  tomato: 0xff6347,
  turquoise: 0x40e0d0,
  violet: 0xee82ee,
  wheat: 0xf5deb3,
  white: 0xffffff,
  whitesmoke: 0xf5f5f5,
  yellow: 0xffff00,
  yellowgreen: 0x9acd32,
};

export type ColorSpec = keyof typeof htmlColors | number;

/**
 * Color representation.
 * @typedef ColorSpec
 * @prop {string} 0xAF10AB - any hex number
 * @prop {string} <html color name>
 */
/**
 
* @typedef ColorschemeSpec
* Built in colorschemes
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
* @prop {Gradient} gradient - Allows the user to provide a gradient to the colorscheme.  Is either a Gradient object or the name of a built-in gradient (rwb, roygb, sinebow)
* @prop {min} - min value for gradient
* @prop {max} - max value for gradient
* @prop {mid} - mid point value for gradient (for rwb)
* @prop {object} map - map of a certain AtomSpec property to a color of the form `{'prop': 'elem', map:elementColors.greenCarbon}` Allows the user to provide a mapping of elements to colors to the colorscheme.  This can be done with any properties, and not just 'elem'.
* @prop {function} colorfunc - Allows the user to provide a function for setting the colorschemes.
*/

/** Preset secondary structure color scheme
 * @struct
 */
export const ssColors = {
  //names are in helix-sheet-coil order
  pyMol: { h: 0xff0000, s: 0xffff00, c: 0x00ff00 },
  Jmol: { h: 0xff0080, s: 0xffc800, c: 0xffffff },
};

const rasmol = {
  H: 0xffffff,
  He: 0xffc0cb,
  HE: 0xffc0cb,
  Li: 0xb22222,
  LI: 0xb22222,
  B: 0x00ff00,
  C: 0xc8c8c8,
  N: 0x8f8fff,
  O: 0xf00000,
  F: 0xdaa520,
  Na: 0x0000ff,
  NA: 0x0000ff,
  Mg: 0x228b22,
  MG: 0x228b22,
  Al: 0x808090,
  AL: 0x808090,
  Si: 0xdaa520,
  SI: 0xdaa520,
  P: 0xffa500,
  S: 0xffc832,
  Cl: 0x00ff00,
  CL: 0x00ff00,
  Ca: 0x808090,
  CA: 0x808090,
  Ti: 0x808090,
  TI: 0x808090,
  Cr: 0x808090,
  CR: 0x808090,
  Mn: 0x808090,
  MN: 0x808090,
  Fe: 0xffa500,
  FE: 0xffa500,
  Ni: 0xa52a2a,
  NI: 0xa52a2a,
  Cu: 0xa52a2a,
  CU: 0xa52a2a,
  Zn: 0xa52a2a,
  ZN: 0xa52a2a,
  Br: 0xa52a2a,
  BR: 0xa52a2a,
  Ag: 0x808090,
  AG: 0x808090,
  I: 0xa020f0,
  Ba: 0xffa500,
  BA: 0xffa500,
  Au: 0xdaa520,
  AU: 0xdaa520,
};

/** Preset element coloring - from individual element colors to entire mappings (e.g. 'elementColors.Jmol' colors atoms with Jmol stylings)
 * @struct
 */
export const elementColors = {
  defaultColor: 0xff1493,
  /** Jmol-like element colors*/
  Jmol: {
    H: 0xffffff,
    He: 0xd9ffff,
    HE: 0xd9ffff,
    Li: 0xcc80ff,
    LI: 0xcc80ff,
    Be: 0xc2ff00,
    BE: 0xc2ff00,
    B: 0xffb5b5,
    C: 0x909090,
    N: 0x3050f8,
    O: 0xff0d0d,
    F: 0x90e050,
    Ne: 0xb3e3f5,
    NE: 0xb3e3f5,
    Na: 0xab5cf2,
    NA: 0xab5cf2,
    Mg: 0x8aff00,
    MG: 0x8aff00,
    Al: 0xbfa6a6,
    AL: 0xbfa6a6,
    Si: 0xf0c8a0,
    SI: 0xf0c8a0,
    P: 0xff8000,
    S: 0xffff30,
    Cl: 0x1ff01f,
    CL: 0x1ff01f,
    Ar: 0x80d1e3,
    AR: 0x80d1e3,
    K: 0x8f40d4,
    Ca: 0x3dff00,
    CA: 0x3dff00,
    Sc: 0xe6e6e6,
    SC: 0xe6e6e6,
    Ti: 0xbfc2c7,
    TI: 0xbfc2c7,
    V: 0xa6a6ab,
    Cr: 0x8a99c7,
    CR: 0x8a99c7,
    Mn: 0x9c7ac7,
    MN: 0x9c7ac7,
    Fe: 0xe06633,
    FE: 0xe06633,
    Co: 0xf090a0,
    CO: 0xf090a0,
    Ni: 0x50d050,
    NI: 0x50d050,
    Cu: 0xc88033,
    CU: 0xc88033,
    Zn: 0x7d80b0,
    ZN: 0x7d80b0,
    Ga: 0xc28f8f,
    GA: 0xc28f8f,
    Ge: 0x668f8f,
    GE: 0x668f8f,
    As: 0xbd80e3,
    AS: 0xbd80e3,
    Se: 0xffa100,
    SE: 0xffa100,
    Br: 0xa62929,
    BR: 0xa62929,
    Kr: 0x5cb8d1,
    KR: 0x5cb8d1,
    Rb: 0x702eb0,
    RB: 0x702eb0,
    Sr: 0x00ff00,
    SR: 0x00ff00,
    Y: 0x94ffff,
    Zr: 0x94e0e0,
    ZR: 0x94e0e0,
    Nb: 0x73c2c9,
    NB: 0x73c2c9,
    Mo: 0x54b5b5,
    MO: 0x54b5b5,
    Tc: 0x3b9e9e,
    TC: 0x3b9e9e,
    Ru: 0x248f8f,
    RU: 0x248f8f,
    Rh: 0x0a7d8c,
    RH: 0x0a7d8c,
    Pd: 0x006985,
    PD: 0x006985,
    Ag: 0xc0c0c0,
    AG: 0xc0c0c0,
    Cd: 0xffd98f,
    CD: 0xffd98f,
    In: 0xa67573,
    IN: 0xa67573,
    Sn: 0x668080,
    SN: 0x668080,
    Sb: 0x9e63b5,
    SB: 0x9e63b5,
    Te: 0xd47a00,
    TE: 0xd47a00,
    I: 0x940094,
    Xe: 0x429eb0,
    XE: 0x429eb0,
    Cs: 0x57178f,
    CS: 0x57178f,
    Ba: 0x00c900,
    BA: 0x00c900,
    La: 0x70d4ff,
    LA: 0x70d4ff,
    Ce: 0xffffc7,
    CE: 0xffffc7,
    Pr: 0xd9ffc7,
    PR: 0xd9ffc7,
    Nd: 0xc7ffc7,
    ND: 0xc7ffc7,
    Pm: 0xa3ffc7,
    PM: 0xa3ffc7,
    Sm: 0x8fffc7,
    SM: 0x8fffc7,
    Eu: 0x61ffc7,
    EU: 0x61ffc7,
    Gd: 0x45ffc7,
    GD: 0x45ffc7,
    Tb: 0x30ffc7,
    TB: 0x30ffc7,
    Dy: 0x1fffc7,
    DY: 0x1fffc7,
    Ho: 0x00ff9c,
    HO: 0x00ff9c,
    Er: 0x00e675,
    ER: 0x00e675,
    Tm: 0x00d452,
    TM: 0x00d452,
    Yb: 0x00bf38,
    YB: 0x00bf38,
    Lu: 0x00ab24,
    LU: 0x00ab24,
    Hf: 0x4dc2ff,
    HF: 0x4dc2ff,
    Ta: 0x4da6ff,
    TA: 0x4da6ff,
    W: 0x2194d6,
    Re: 0x267dab,
    RE: 0x267dab,
    Os: 0x266696,
    OS: 0x266696,
    Ir: 0x175487,
    IR: 0x175487,
    Pt: 0xd0d0e0,
    PT: 0xd0d0e0,
    Au: 0xffd123,
    AU: 0xffd123,
    Hg: 0xb8b8d0,
    HG: 0xb8b8d0,
    Tl: 0xa6544d,
    TL: 0xa6544d,
    Pb: 0x575961,
    PB: 0x575961,
    Bi: 0x9e4fb5,
    BI: 0x9e4fb5,
    Po: 0xab5c00,
    PO: 0xab5c00,
    At: 0x754f45,
    AT: 0x754f45,
    Rn: 0x428296,
    RN: 0x428296,
    Fr: 0x420066,
    FR: 0x420066,
    Ra: 0x007d00,
    RA: 0x007d00,
    Ac: 0x70abfa,
    AC: 0x70abfa,
    Th: 0x00baff,
    TH: 0x00baff,
    Pa: 0x00a1ff,
    PA: 0x00a1ff,
    U: 0x008fff,
    Np: 0x0080ff,
    NP: 0x0080ff,
    Pu: 0x006bff,
    PU: 0x006bff,
    Am: 0x545cf2,
    AM: 0x545cf2,
    Cm: 0x785ce3,
    CM: 0x785ce3,
    Bk: 0x8a4fe3,
    BK: 0x8a4fe3,
    Cf: 0xa136d4,
    CF: 0xa136d4,
    Es: 0xb31fd4,
    ES: 0xb31fd4,
    Fm: 0xb31fba,
    FM: 0xb31fba,
    Md: 0xb30da6,
    MD: 0xb30da6,
    No: 0xbd0d87,
    NO: 0xbd0d87,
    Lr: 0xc70066,
    LR: 0xc70066,
    Rf: 0xcc0059,
    RF: 0xcc0059,
    Db: 0xd1004f,
    DB: 0xd1004f,
    Sg: 0xd90045,
    SG: 0xd90045,
    Bh: 0xe00038,
    BH: 0xe00038,
    Hs: 0xe6002e,
    HS: 0xe6002e,
    Mt: 0xeb0026,
    MT: 0xeb0026,
  },
  /** rasmol-like element colors */
  rasmol,
  defaultColors: {
    ...rasmol,
  },
  greenCarbon: {
    ...rasmol,
    C: 0x00ff00,
  },
  cyanCarbon: {
    ...rasmol,
    C: 0x00ffff,
  },
  magentaCarbon: {
    ...rasmol,
    C: 0xff00ff,
  },
  yellowCarbon: {
    ...rasmol,
    C: 0xffff00,
  },
  whiteCarbon: {
    ...rasmol,
    C: 0xffffff,
  },
  orangeCarbon: {
    ...rasmol,
    C: 0xffa500,
  },
  purpleCarbon: {
    ...rasmol,
    C: 0x800080,
  },
  blueCarbon: {
    ...rasmol,
    C: 0x0000ff,
  },
};

export const residues = {
  /** @property standard amino acid color scheme*/

  amino: {
    ALA: 0xc8c8c8,
    ARG: 0x145aff,
    ASN: 0x00dcdc,
    ASP: 0xe60a0a,
    CYS: 0xe6e600,
    GLN: 0x00dcdc,
    GLU: 0xe60a0a,
    GLY: 0xebebeb,
    HIS: 0x8282d2,
    ILE: 0x0f820f,
    LEU: 0x0f820f,
    LYS: 0x145aff,
    MET: 0xe6e600,
    PHE: 0x3232aa,
    PRO: 0xdc9682,
    SER: 0xfa9600,
    THR: 0xfa9600,
    TRP: 0xb45ab4,
    TYR: 0x3232aa,
    VAL: 0x0f820f,
    ASX: 0xff69b4,
    GLX: 0xff69b4,
  },

  /** @property shapely amino acid color scheme*/
  shapely: {
    ALA: 0x8cff8c,
    ARG: 0x00007c,
    ASN: 0xff7c70,
    ASP: 0xa00042,
    CYS: 0xffff70,
    GLN: 0xff4c4c,
    GLU: 0x660000,
    GLY: 0xffffff,
    HIS: 0x7070ff,
    ILE: 0x004c00,
    LEU: 0x455e45,
    LYS: 0x4747b8,
    MET: 0xb8a042,
    PHE: 0x534c52,
    PRO: 0x525252,
    SER: 0xff7042,
    THR: 0xb84c00,
    TRP: 0x4f4600,
    TYR: 0x8c704c,
    VAL: 0xff8cff,
    ASX: 0xff00ff,
    GLX: 0xff00ff,
  },

  /** @property nucleic acid color scheme*/
  nucleic: {
    A: 0xa0a0ff,
    G: 0xff7070,
    I: 0x80ffff,
    C: 0xff8c4b,
    T: 0xa0ffa0,
    U: 0xff8080,
  },
};

export const chains = {
  /** @property chain based standard color scheme */
  atom: {
    A: 0xc0d0ff,
    B: 0xb0ffb0,
    C: 0xffc0c8,
    D: 0xffff80,
    E: 0xffc0ff,
    F: 0xb0f0f0,
    G: 0xffd070,
    H: 0xf08080,
    I: 0xf5deb3,
    J: 0x00bfff,
    K: 0xcd5c5c,
    L: 0x66cdaa,
    M: 0x9acd32,
    N: 0xee82ee,
    O: 0x00ced1,
    P: 0x00ff7f,
    Q: 0x3cb371,
    R: 0x00008b,
    S: 0xbdb76b,
    T: 0x006400,
    U: 0x800000,
    V: 0x808000,
    W: 0x800080,
    X: 0x008080,
    Y: 0xb8860b,
    Z: 0xb22222,
  },

  /** @property hetatm color scheme */
  hetatm: {
    A: 0x90a0cf,
    B: 0x80cf98,
    C: 0xcf90b0,
    D: 0xcfcf70,
    E: 0xcf90cf,
    F: 0x80c0c0,
    G: 0xcfa060,
    H: 0xc05070,
    I: 0xc5ae83,
    J: 0x00a7cf,
    K: 0xb54c4c,
    L: 0x56b592,
    M: 0x8ab52a,
    N: 0xbe72be,
    O: 0x00b6a1,
    P: 0x00cf6f,
    Q: 0x349b61,
    R: 0x0000bb,
    S: 0xa59f5b,
    T: 0x009400,
    U: 0xb00000,
    V: 0xb0b000,
    W: 0xb000b0,
    X: 0x00b0b0,
    Y: 0xe8b613,
    Z: 0xc23232,
  },
};

/**
 * built in color schemes
 * The user can pass all of these values directly as the colorscheme and they will use the respective colorscheme 
 * @example window.$3Dmol.download("pdb:4UAA",viewer,{},function(){
 *    viewer.setBackgroundColor(0xffffffff);
 *    var colorAsSnake = function(atom) {
 *      return atom.resi % 2 ? 'white': 'green'
 *    };
 *    viewer.setStyle( {chain:'A'}, { cartoon: {colorfunc: colorAsSnake }});
 *    viewer.setStyle( {chain:'B'}, { stick: {colorscheme: 'yellowCarbon'}});
 *    viewer.render();
 *  });
  */
export const builtinColorSchemes = {
  ssPyMol: { prop: "ss", map: ssColors.pyMol },
  ssJmol: { prop: "ss", map: ssColors.Jmol },
  Jmol: { prop: "elem", map: elementColors.Jmol },
  amino: { prop: "resn", map: residues.amino },
  shapely: { prop: "resn", map: residues.shapely },
  nucleic: { prop: "resn", map: residues.nucleic },
  chain: { prop: "chain", map: chains.atom },
  rasmol: { prop: "elem", map: elementColors.rasmol },
  default: { prop: "elem", map: elementColors.defaultColors },
  greenCarbon: { prop: "elem", map: elementColors.greenCarbon },
  chainHetatm: { prop: "chain", map: chains.hetatm },
  cyanCarbon: { prop: "elem", map: elementColors.cyanCarbon },
  magentaCarbon: { prop: "elem", map: elementColors.magentaCarbon },
  purpleCarbon: { prop: "elem", map: elementColors.purpleCarbon },
  whiteCarbon: { prop: "elem", map: elementColors.whiteCarbon },
  orangeCarbon: { prop: "elem", map: elementColors.orangeCarbon },
  yellowCarbon: { prop: "elem", map: elementColors.yellowCarbon },
  blueCarbon: { prop: "elem", map: elementColors.blueCarbon },
};

/** Return proper color for atom given style
 * @param {AtomSpec} atom
 * @param {AtomStyle} style
 * @return {Color}
 */
export function getColorFromStyle(atom, style) {
  var scheme = style.colorscheme;
  if (typeof builtinColorSchemes[scheme] != "undefined") {
    scheme = builtinColorSchemes[scheme];
  } else if (typeof scheme == "string" && scheme.endsWith("Carbon")) {
    //any color you want of carbon
    var ccolor = scheme
      .substring(0, scheme.lastIndexOf("Carbon"))
      .toLowerCase();
    if (typeof htmlColors[ccolor] != "undefined") {
      var newscheme = { ...elementColors.defaultColors };
      newscheme.C = htmlColors[ccolor];
      builtinColorSchemes[scheme] = { prop: "elem", map: newscheme };
      scheme = builtinColorSchemes[scheme];
    }
  }

  var color = atom.color;
  if (typeof style.color != "undefined" && style.color != "spectrum")
    color = style.color;
  if (typeof scheme != "undefined") {
    var prop, val;
    if (typeof elementColors[scheme] != "undefined") {
      //name of builtin colorscheme
      scheme = elementColors[scheme];
      if (typeof scheme[atom[scheme.prop]] != "undefined") {
        color = scheme.map[atom[scheme.prop]];
      }
    } else if (typeof scheme[atom[scheme.prop]] != "undefined") {
      //actual color scheme provided
      color = scheme.map[atom[scheme.prop]];
    } else if (
      typeof scheme.prop != "undefined" &&
      typeof scheme.gradient != "undefined"
    ) {
      //apply a property mapping
      prop = scheme.prop;
      var grad = scheme.gradient; //redefining scheme
      if (typeof builtinGradients[grad] != "undefined") {
        grad = new builtinGradients[grad](
          scheme.min,
          scheme.max,
          scheme.mid
        );
      }

      var range = grad.range() || [-1, 1]; //sensible default
      val = getAtomProperty(atom, prop);
      if (val != null) {
        color = grad.valueToHex(val, range);
      }
    } else if (
      typeof scheme.prop != "undefined" &&
      typeof scheme.map != "undefined"
    ) {
      //apply a discrete property mapping
      prop = scheme.prop;
      val = getAtomProperty(atom, prop);
      if (typeof scheme.map[val] != "undefined") {
        color = scheme.map[val];
      }
    } else if (typeof style.colorscheme[atom.elem] != "undefined") {
      //actual color scheme provided
      color = style.colorscheme[atom.elem];
    } else {
      console.log("Could not interpret colorscheme " + scheme);
    }
  } else if (typeof style.colorfunc != "undefined") {
    //this is a user provided function for turning an atom into a color
    color = style.colorfunc(atom);
  }

  var C = CC.color(color);
  return C;
};
