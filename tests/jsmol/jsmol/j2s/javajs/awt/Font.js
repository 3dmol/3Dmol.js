Clazz.declarePackage ("javajs.awt");
Clazz.load (null, "javajs.awt.Font", ["JU.AU"], function () {
c$ = Clazz.decorateAsClass (function () {
this.fid = 0;
this.fontFace = null;
this.fontStyle = null;
this.fontSizeNominal = 0;
this.idFontFace = 0;
this.idFontStyle = 0;
this.fontSize = 0;
this.font = null;
this.fontMetrics = null;
this.manager = null;
this.ascent = 0;
this.descent = 0;
this.isBold = false;
this.isItalic = false;
Clazz.instantialize (this, arguments);
}, javajs.awt, "Font");
Clazz.makeConstructor (c$, 
 function (manager, fid, idFontFace, idFontStyle, fontSize, fontSizeNominal, graphics) {
this.manager = manager;
this.fid = fid;
this.fontFace = javajs.awt.Font.fontFaces[idFontFace];
this.fontStyle = javajs.awt.Font.fontStyles[idFontStyle];
this.idFontFace = idFontFace;
this.idFontStyle = idFontStyle;
this.fontSize = fontSize;
this.isBold = (idFontStyle & 1) == 1;
this.isItalic = (idFontStyle & 2) == 2;
this.fontSizeNominal = fontSizeNominal;
this.font = manager.newFont (javajs.awt.Font.fontFaces[idFontFace], this.isBold, this.isItalic, fontSize);
this.fontMetrics = manager.getFontMetrics (this, graphics);
this.descent = manager.getFontDescent (this.fontMetrics);
this.ascent = manager.getFontAscent (this.fontMetrics);
}, "javajs.api.FontManager,~N,~N,~N,~N,~N,~O");
c$.getFont3D = Clazz.defineMethod (c$, "getFont3D", 
function (fontID) {
return javajs.awt.Font.font3ds[fontID & 0xFF];
}, "~N");
c$.createFont3D = Clazz.defineMethod (c$, "createFont3D", 
function (fontface, fontstyle, fontsize, fontsizeNominal, manager, graphicsForMetrics) {
if (fontsize > 0xFF) fontsize = 0xFF;
var fontsizeX16 = (Clazz.floatToInt (fontsize)) << 4;
var fontkey = ((fontface & 3) | ((fontstyle & 3) << 2) | (fontsizeX16 << 4));
for (var i = javajs.awt.Font.fontkeyCount; --i > 0; ) if (fontkey == javajs.awt.Font.fontkeys[i] && javajs.awt.Font.font3ds[i].fontSizeNominal == fontsizeNominal) return javajs.awt.Font.font3ds[i];

var fontIndexNext = javajs.awt.Font.fontkeyCount++;
if (fontIndexNext == javajs.awt.Font.fontkeys.length) javajs.awt.Font.fontkeys = JU.AU.arrayCopyI (javajs.awt.Font.fontkeys, fontIndexNext + 8);
javajs.awt.Font.font3ds = JU.AU.arrayCopyObject (javajs.awt.Font.font3ds, fontIndexNext + 8);
var font3d =  new javajs.awt.Font (manager, fontIndexNext, fontface, fontstyle, fontsize, fontsizeNominal, graphicsForMetrics);
javajs.awt.Font.font3ds[fontIndexNext] = font3d;
javajs.awt.Font.fontkeys[fontIndexNext] = fontkey;
return font3d;
}, "~N,~N,~N,~N,javajs.api.FontManager,~O");
c$.getFontFaceID = Clazz.defineMethod (c$, "getFontFaceID", 
function (fontface) {
return ("Monospaced".equalsIgnoreCase (fontface) ? 2 : "Serif".equalsIgnoreCase (fontface) ? 1 : 0);
}, "~S");
c$.getFontStyleID = Clazz.defineMethod (c$, "getFontStyleID", 
function (fontstyle) {
for (var i = 4; --i >= 0; ) if (javajs.awt.Font.fontStyles[i].equalsIgnoreCase (fontstyle)) return i;

return -1;
}, "~S");
Clazz.defineMethod (c$, "getAscent", 
function () {
return this.ascent;
});
Clazz.defineMethod (c$, "getDescent", 
function () {
return this.descent;
});
Clazz.defineMethod (c$, "getHeight", 
function () {
return this.getAscent () + this.getDescent ();
});
Clazz.defineMethod (c$, "getFontMetrics", 
function () {
return this.fontMetrics;
});
Clazz.defineMethod (c$, "stringWidth", 
function (text) {
return this.manager.fontStringWidth (this, text);
}, "~S");
Clazz.defineMethod (c$, "getInfo", 
function () {
return this.fontSizeNominal + " " + this.fontFace + " " + this.fontStyle;
});
Clazz.defineStatics (c$,
"FONT_ALLOCATION_UNIT", 8,
"fontkeyCount", 1,
"fontkeys",  Clazz.newIntArray (8, 0));
c$.font3ds = c$.prototype.font3ds =  new Array (8);
Clazz.defineStatics (c$,
"FONT_FACE_SANS", 0,
"FONT_FACE_SERIF", 1,
"FONT_FACE_MONO", 2,
"fontFaces", ["SansSerif", "Serif", "Monospaced", ""],
"FONT_STYLE_PLAIN", 0,
"FONT_STYLE_BOLD", 1,
"FONT_STYLE_ITALIC", 2,
"FONT_STYLE_BOLDITALIC", 3,
"fontStyles", ["Plain", "Bold", "Italic", "BoldItalic"]);
});
