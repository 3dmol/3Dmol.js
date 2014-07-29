Clazz.declarePackage ("J.g3d");
Clazz.load (["java.util.Hashtable"], "J.g3d.TextRenderer", ["JU.CU"], function () {
c$ = Clazz.decorateAsClass (function () {
this.height = 0;
this.ascent = 0;
this.width = 0;
this.mapWidth = 0;
this.size = 0;
this.tmap = null;
this.isInvalid = false;
Clazz.instantialize (this, arguments);
}, J.g3d, "TextRenderer");
c$.clearFontCache = Clazz.defineMethod (c$, "clearFontCache", 
function () {
if (J.g3d.TextRenderer.working) return;
J.g3d.TextRenderer.htFont3d.clear ();
J.g3d.TextRenderer.htFont3dAntialias.clear ();
});
c$.plot = Clazz.defineMethod (c$, "plot", 
function (x, y, z, argb, bgargb, text, font3d, g3d, jmolRenderer, antialias) {
if (text.length == 0) return 0;
if (text.indexOf ("<su") >= 0 || text.indexOf ("<color") >= 0) return J.g3d.TextRenderer.plotByCharacter (x, y, z, argb, bgargb, text, font3d, g3d, jmolRenderer, antialias);
var offset = font3d.getAscent ();
y -= offset;
var text3d = J.g3d.TextRenderer.getPlotText3D (x, y, g3d, text, font3d, antialias);
if (text3d.isInvalid) return text3d.width;
if (antialias && (argb & 0xC0C0C0) == 0) {
argb = argb | 0x040404;
}if (jmolRenderer != null || (x < 0 || x + text3d.width > g3d.width || y < 0 || y + text3d.height > g3d.height)) J.g3d.TextRenderer.plotClipped (x, y, z, argb, bgargb, g3d, jmolRenderer, text3d.mapWidth, text3d.height, text3d.tmap);
 else J.g3d.TextRenderer.plotUnclipped (x, y, z, argb, bgargb, g3d, text3d.mapWidth, text3d.height, text3d.tmap);
return text3d.width;
}, "~N,~N,~N,~N,~N,~S,javajs.awt.Font,J.g3d.Graphics3D,J.api.JmolRendererInterface,~B");
c$.plotByCharacter = Clazz.defineMethod (c$, "plotByCharacter", 
 function (x, y, z, argb, bgargb, text, font3d, g3d, jmolRenderer, antialias) {
var w = 0;
var len = text.length;
var suboffset = Math.round (font3d.getHeight () * 0.25);
var supoffset = -Math.round (font3d.getHeight () * 0.3);
var argb0 = 0;
for (var i = 0; i < len; i++) {
if (text.charAt (i) == '<') {
if (i + 5 < len && text.substring (i, i + 6).equals ("<color")) {
argb0 = argb;
var pt = text.indexOf (">", i);
if (pt < 0) continue;
argb = JU.CU.getArgbFromString (text.substring (i + 7, pt).trim ());
i = pt;
continue;
}if (i + 7 < len && text.substring (i, i + 8).equals ("</color>")) {
i += 7;
argb = argb0;
continue;
}if (i + 4 < len && text.substring (i, i + 5).equals ("<sub>")) {
i += 4;
y += suboffset;
continue;
}if (i + 4 < len && text.substring (i, i + 5).equals ("<sup>")) {
i += 4;
y += supoffset;
continue;
}if (i + 5 < len && text.substring (i, i + 6).equals ("</sub>")) {
i += 5;
y -= suboffset;
continue;
}if (i + 5 < len && text.substring (i, i + 6).equals ("</sup>")) {
i += 5;
y -= supoffset;
continue;
}}var width = J.g3d.TextRenderer.plot (x + w, y, z, argb, bgargb, text.substring (i, i + 1), font3d, g3d, jmolRenderer, antialias);
w += width;
}
return w;
}, "~N,~N,~N,~N,~N,~S,javajs.awt.Font,J.g3d.Graphics3D,J.api.JmolRendererInterface,~B");
c$.plotUnclipped = Clazz.defineMethod (c$, "plotUnclipped", 
 function (x, y, z, argb, bgargb, g3d, textWidth, textHeight, tmap) {
var offset = 0;
var zbuf = g3d.zbuf;
var renderWidth = g3d.width;
var pbufOffset = y * renderWidth + x;
for (var i = 0; i < textHeight; i++) {
for (var j = 0; j < textWidth; j++) {
var shade = tmap[offset++];
if (shade != 0 && z < zbuf[pbufOffset]) g3d.shadeTextPixel (pbufOffset, z, argb, bgargb, shade);
pbufOffset++;
}
pbufOffset += (renderWidth - textWidth);
}
}, "~N,~N,~N,~N,~N,J.g3d.Graphics3D,~N,~N,~A");
c$.plotClipped = Clazz.defineMethod (c$, "plotClipped", 
 function (x, y, z, argb, bgargb, g3d, jmolRenderer, textWidth, textHeight, tmap) {
if (jmolRenderer == null) jmolRenderer = g3d;
var offset = 0;
for (var i = 0; i < textHeight; i++) {
for (var j = 0; j < textWidth; j++) {
var shade = tmap[offset++];
if (shade != 0) jmolRenderer.plotImagePixel (argb, x + j, y + i, z, shade, bgargb);
}
}
}, "~N,~N,~N,~N,~N,J.g3d.Graphics3D,J.api.JmolRendererInterface,~N,~N,~A");
Clazz.makeConstructor (c$, 
 function (text, font3d) {
this.ascent = font3d.getAscent ();
this.height = font3d.getHeight ();
this.width = font3d.stringWidth (text);
if (this.width == 0) return;
this.mapWidth = this.width;
this.size = this.mapWidth * this.height;
}, "~S,javajs.awt.Font");
c$.getPlotText3D = Clazz.defineMethod (c$, "getPlotText3D", 
 function (x, y, g3d, text, font3d, antialias) {
J.g3d.TextRenderer.working = true;
var ht = (antialias ? J.g3d.TextRenderer.htFont3dAntialias : J.g3d.TextRenderer.htFont3d);
var htForThisFont = ht.get (font3d);
var text3d = null;
var newFont = false;
var newText = false;
if (htForThisFont != null) {
text3d = htForThisFont.get (text);
} else {
htForThisFont =  new java.util.Hashtable ();
newFont = true;
}if (text3d == null) {
text3d =  new J.g3d.TextRenderer (text, font3d);
newText = true;
}text3d.isInvalid = (text3d.width == 0 || x + text3d.width <= 0 || x >= g3d.width || y + text3d.height <= 0 || y >= g3d.height);
if (text3d.isInvalid) return text3d;
if (newFont) ht.put (font3d, htForThisFont);
if (newText) {
text3d.setTranslucency (text, font3d, g3d);
htForThisFont.put (text, text3d);
}J.g3d.TextRenderer.working = false;
return text3d;
}, "~N,~N,J.g3d.Graphics3D,~S,javajs.awt.Font,~B");
Clazz.defineMethod (c$, "setTranslucency", 
 function (text, font3d, g3d) {
var pixels = g3d.apiPlatform.getTextPixels (text, font3d, g3d.platform.getGraphicsForTextOrImage (this.mapWidth, this.height), g3d.platform.offscreenImage, this.mapWidth, this.height, this.ascent);
if (pixels == null) return;
this.tmap =  Clazz.newByteArray (this.size, 0);
for (var i = pixels.length; --i >= 0; ) {
var p = pixels[i] & 0xFF;
if (p != 0) {
this.tmap[i] = J.g3d.TextRenderer.translucency[p >> 5];
}}
}, "~S,javajs.awt.Font,J.g3d.Graphics3D");
Clazz.defineStatics (c$,
"translucency", [7, 6, 5, 4, 3, 2, 1, 8],
"working", false);
c$.htFont3d = c$.prototype.htFont3d =  new java.util.Hashtable ();
c$.htFont3dAntialias = c$.prototype.htFont3dAntialias =  new java.util.Hashtable ();
});
