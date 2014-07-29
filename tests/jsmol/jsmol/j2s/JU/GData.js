Clazz.declarePackage ("JU");
Clazz.load (["J.api.JmolGraphicsInterface"], "JU.GData", ["javajs.awt.Font", "JU.AU", "$.P3", "$.V3", "JU.C", "$.Shader"], function () {
c$ = Clazz.decorateAsClass (function () {
this.apiPlatform = null;
this.translucentCoverOnly = false;
this.windowWidth = 0;
this.windowHeight = 0;
this.displayMinX = 0;
this.displayMaxX = 0;
this.displayMinY = 0;
this.displayMaxY = 0;
this.antialiasThisFrame = false;
this.antialiasEnabled = false;
this.inGreyscaleMode = false;
this.changeableColixMap = null;
this.backgroundImage = null;
this.newWindowWidth = 0;
this.newWindowHeight = 0;
this.newAntialiasing = false;
this.bgcolor = 0;
this.xLast = 0;
this.yLast = 0;
this.slab = 0;
this.depth = 0;
this.width = 0;
this.height = 0;
this.zSlab = 0;
this.zDepth = 0;
this.zShadePower = 3;
this.colixCurrent = 0;
this.argbCurrent = 0;
this.bufferSize = 0;
this.shader = null;
this.zShadeR = 0;
this.zShadeG = 0;
this.zShadeB = 0;
this.graphicsForMetrics = null;
this.$isPass2 = false;
this.ambientOcclusion = 0;
Clazz.instantialize (this, arguments);
}, JU, "GData", null, J.api.JmolGraphicsInterface);
Clazz.prepareFields (c$, function () {
this.changeableColixMap =  Clazz.newShortArray (16, 0);
});
Clazz.defineMethod (c$, "setTranslucentCoverOnly", 
function (TF) {
this.translucentCoverOnly = TF;
}, "~B");
Clazz.defineMethod (c$, "getTranslucentCoverOnly", 
function () {
return this.translucentCoverOnly;
});
Clazz.makeConstructor (c$, 
function () {
this.shader =  new JU.Shader ();
});
Clazz.defineMethod (c$, "initialize", 
function (apiPlatform) {
this.apiPlatform = apiPlatform;
}, "javajs.api.GenericPlatform");
Clazz.overrideMethod (c$, "setDepth", 
function (depthValue) {
this.depth = depthValue < 0 ? 0 : depthValue;
}, "~N");
Clazz.overrideMethod (c$, "setSlab", 
function (slabValue) {
this.slab = slabValue < 0 ? 0 : slabValue;
}, "~N");
Clazz.defineMethod (c$, "setZShade", 
function (zShade, zSlab, zDepth, zPower) {
if (zShade) this.setZShade2 (zSlab, zDepth, zPower);
}, "~B,~N,~N,~N");
Clazz.defineMethod (c$, "setZShade2", 
function (zSlab, zDepth, zPower) {
this.zShadeR = this.bgcolor & 0xFF;
this.zShadeG = (this.bgcolor & 0xFF00) >> 8;
this.zShadeB = (this.bgcolor & 0xFF0000) >> 16;
this.zSlab = zSlab < 0 ? 0 : zSlab;
this.zDepth = zDepth < 0 ? 0 : zDepth;
this.zShadePower = zPower;
}, "~N,~N,~N");
Clazz.overrideMethod (c$, "setAmbientOcclusion", 
function (value) {
this.ambientOcclusion = value;
}, "~N");
Clazz.overrideMethod (c$, "getRenderWidth", 
function () {
return this.width;
});
Clazz.overrideMethod (c$, "getRenderHeight", 
function () {
return this.height;
});
Clazz.overrideMethod (c$, "getSlab", 
function () {
return this.slab;
});
Clazz.overrideMethod (c$, "getDepth", 
function () {
return this.depth;
});
Clazz.defineMethod (c$, "isDisplayAntialiased", 
function () {
return this.antialiasEnabled;
});
Clazz.overrideMethod (c$, "isAntialiased", 
function () {
return this.antialiasThisFrame;
});
Clazz.defineMethod (c$, "getChangeableColix", 
function (id, argb) {
if (id >= this.changeableColixMap.length) this.changeableColixMap = JU.AU.arrayCopyShort (this.changeableColixMap, id + 16);
if (this.changeableColixMap[id] == 0) this.changeableColixMap[id] = JU.C.getColix (argb);
return (id | -32768);
}, "~N,~N");
Clazz.defineMethod (c$, "changeColixArgb", 
function (id, argb) {
if (id < this.changeableColixMap.length && this.changeableColixMap[id] != 0) this.changeableColixMap[id] = JU.C.getColix (argb);
}, "~N,~N");
Clazz.defineMethod (c$, "getBgColixes", 
function (bgcolixes) {
return bgcolixes;
}, "~A");
Clazz.overrideMethod (c$, "getColorArgbOrGray", 
function (colix) {
if (colix < 0) colix = this.changeableColixMap[colix & 2047];
return (this.inGreyscaleMode ? JU.C.getArgbGreyscale (colix) : JU.C.getArgb (colix));
}, "~N");
Clazz.defineMethod (c$, "getShades", 
function (colix) {
if (colix < 0) colix = this.changeableColixMap[colix & 2047];
return (this.inGreyscaleMode ? this.shader.getShadesG (colix) : this.shader.getShades (colix));
}, "~N");
Clazz.defineMethod (c$, "setGreyscaleMode", 
function (greyscaleMode) {
this.inGreyscaleMode = greyscaleMode;
}, "~B");
Clazz.defineMethod (c$, "getSpecularPower", 
function () {
return this.shader.specularPower;
});
Clazz.defineMethod (c$, "setSpecularPower", 
function (val) {
if (val < 0) {
this.setSpecularExponent (-val);
return;
}if (this.shader.specularPower == val) return;
this.shader.specularPower = val;
this.shader.intenseFraction = val / 100;
this.shader.flushCaches ();
}, "~N");
Clazz.defineMethod (c$, "getSpecularPercent", 
function () {
return this.shader.specularPercent;
});
Clazz.defineMethod (c$, "setSpecularPercent", 
function (val) {
if (this.shader.specularPercent == val) return;
this.shader.specularPercent = val;
this.shader.specularFactor = val / 100;
this.shader.flushCaches ();
}, "~N");
Clazz.defineMethod (c$, "getSpecularExponent", 
function () {
return this.shader.specularExponent;
});
Clazz.defineMethod (c$, "setSpecularExponent", 
function (val) {
if (this.shader.specularExponent == val) return;
this.shader.specularExponent = val;
this.shader.phongExponent = Clazz.doubleToInt (Math.pow (2, val));
this.shader.usePhongExponent = false;
this.shader.flushCaches ();
}, "~N");
Clazz.defineMethod (c$, "getPhongExponent", 
function () {
return this.shader.phongExponent;
});
Clazz.defineMethod (c$, "setPhongExponent", 
function (val) {
if (this.shader.phongExponent == val && this.shader.usePhongExponent) return;
this.shader.phongExponent = val;
var x = (Math.log (val) / Math.log (2));
this.shader.usePhongExponent = (x != Clazz.floatToInt (x));
if (!this.shader.usePhongExponent) this.shader.specularExponent = Clazz.floatToInt (x);
this.shader.flushCaches ();
}, "~N");
Clazz.defineMethod (c$, "getDiffusePercent", 
function () {
return this.shader.diffusePercent;
});
Clazz.defineMethod (c$, "setDiffusePercent", 
function (val) {
if (this.shader.diffusePercent == val) return;
this.shader.diffusePercent = val;
this.shader.diffuseFactor = val / 100;
this.shader.flushCaches ();
}, "~N");
Clazz.defineMethod (c$, "getAmbientPercent", 
function () {
return this.shader.ambientPercent;
});
Clazz.defineMethod (c$, "setAmbientPercent", 
function (val) {
if (this.shader.ambientPercent == val) return;
this.shader.ambientPercent = val;
this.shader.ambientFraction = val / 100;
this.shader.flushCaches ();
}, "~N");
Clazz.defineMethod (c$, "getSpecular", 
function () {
return this.shader.specularOn;
});
Clazz.defineMethod (c$, "setSpecular", 
function (val) {
if (this.shader.specularOn == val) return;
this.shader.specularOn = val;
this.shader.flushCaches ();
}, "~B");
Clazz.defineMethod (c$, "setCel", 
function (val) {
this.shader.setCel (val, this.shader.celPower, this.bgcolor);
}, "~B");
Clazz.defineMethod (c$, "getCel", 
function () {
return this.shader.celOn;
});
Clazz.defineMethod (c$, "getCelPower", 
function () {
return this.shader.celPower;
});
Clazz.defineMethod (c$, "setCelPower", 
function (celPower) {
this.shader.setCel (this.shader.celOn || this.shader.celPower == 0, celPower, this.bgcolor);
}, "~N");
Clazz.defineMethod (c$, "getLightSource", 
function () {
return this.shader.lightSource;
});
Clazz.defineMethod (c$, "isClipped3", 
function (x, y, z) {
return (x < 0 || x >= this.width || y < 0 || y >= this.height || z < this.slab || z > this.depth);
}, "~N,~N,~N");
Clazz.defineMethod (c$, "isClipped", 
function (x, y) {
return (x < 0 || x >= this.width || y < 0 || y >= this.height);
}, "~N,~N");
Clazz.overrideMethod (c$, "isInDisplayRange", 
function (x, y) {
return (x >= this.displayMinX && x < this.displayMaxX && y >= this.displayMinY && y < this.displayMaxY);
}, "~N,~N");
Clazz.overrideMethod (c$, "isClippedXY", 
function (diameter, x, y) {
var r = (diameter + 1) >> 1;
return (x < -r || x >= this.width + r || y < -r || y >= this.height + r);
}, "~N,~N,~N");
Clazz.overrideMethod (c$, "isClippedZ", 
function (z) {
return (z != -2147483648 && (z < this.slab || z > this.depth));
}, "~N");
Clazz.defineMethod (c$, "clipCode3", 
function (x, y, z) {
var code = 0;
if (x < 0) code |= 8;
 else if (x >= this.width) code |= 4;
if (y < 0) code |= 2;
 else if (y >= this.height) code |= 1;
if (z < this.slab) code |= 32;
 else if (z > this.depth) code |= 16;
return code;
}, "~N,~N,~N");
Clazz.defineMethod (c$, "clipCode", 
function (z) {
var code = 0;
if (z < this.slab) code |= 32;
 else if (z > this.depth) code |= 16;
return code;
}, "~N");
Clazz.defineMethod (c$, "getFont3D", 
function (fontSize) {
return javajs.awt.Font.createFont3D (0, 0, fontSize, fontSize, this.apiPlatform, this.graphicsForMetrics);
}, "~N");
Clazz.defineMethod (c$, "getFont3DFS", 
function (fontFace, fontSize) {
return javajs.awt.Font.createFont3D (javajs.awt.Font.getFontFaceID (fontFace), 0, fontSize, fontSize, this.apiPlatform, this.graphicsForMetrics);
}, "~S,~N");
Clazz.defineMethod (c$, "getFontFidFS", 
function (fontFace, fontSize) {
return this.getFont3DFSS (fontFace, "Bold", fontSize).fid;
}, "~S,~N");
Clazz.defineMethod (c$, "getFont3DFSS", 
function (fontFace, fontStyle, fontSize) {
var iStyle = javajs.awt.Font.getFontStyleID (fontStyle);
if (iStyle < 0) iStyle = 0;
return javajs.awt.Font.createFont3D (javajs.awt.Font.getFontFaceID (fontFace), iStyle, fontSize, fontSize, this.apiPlatform, this.graphicsForMetrics);
}, "~S,~S,~N");
Clazz.overrideMethod (c$, "getFont3DScaled", 
function (font, scale) {
var newScale = font.fontSizeNominal * scale;
return (newScale == font.fontSize ? font : javajs.awt.Font.createFont3D (font.idFontFace, font.idFontStyle, newScale, font.fontSizeNominal, this.apiPlatform, this.graphicsForMetrics));
}, "javajs.awt.Font,~N");
Clazz.overrideMethod (c$, "getFontFid", 
function (fontSize) {
return this.getFont3D (fontSize).fid;
}, "~N");
c$.getFontStyleID = Clazz.defineMethod (c$, "getFontStyleID", 
function (fontStyle) {
return javajs.awt.Font.getFontStyleID (fontStyle);
}, "~S");
Clazz.defineMethod (c$, "setBackgroundTransparent", 
function (TF) {
}, "~B");
Clazz.defineMethod (c$, "setBackgroundArgb", 
function (argb) {
this.bgcolor = argb;
this.setCel (this.shader.celOn);
}, "~N");
Clazz.defineMethod (c$, "setBackgroundImage", 
function (image) {
this.backgroundImage = image;
}, "~O");
Clazz.defineMethod (c$, "setWindowParameters", 
function (width, height, antialias) {
this.setWinParams (width, height, antialias);
}, "~N,~N,~B");
Clazz.defineMethod (c$, "setWinParams", 
function (width, height, antialias) {
this.newWindowWidth = width;
this.newWindowHeight = height;
this.newAntialiasing = antialias;
}, "~N,~N,~B");
Clazz.defineMethod (c$, "setNewWindowParametersForExport", 
function () {
this.windowWidth = this.newWindowWidth;
this.windowHeight = this.newWindowHeight;
this.setWidthHeight (false);
});
Clazz.defineMethod (c$, "setWidthHeight", 
function (isAntialiased) {
this.width = this.windowWidth;
this.height = this.windowHeight;
if (isAntialiased) {
this.width <<= 1;
this.height <<= 1;
}this.xLast = this.width - 1;
this.yLast = this.height - 1;
this.displayMinX = -(this.width >> 1);
this.displayMaxX = this.width - this.displayMinX;
this.displayMinY = -(this.height >> 1);
this.displayMaxY = this.height - this.displayMinY;
this.bufferSize = this.width * this.height;
}, "~B");
Clazz.defineMethod (c$, "beginRendering", 
function (stereoRotationMatrix, translucentMode, isImageWrite, renderLow) {
}, "JU.M3,~B,~B,~B");
Clazz.defineMethod (c$, "endRendering", 
function () {
});
Clazz.defineMethod (c$, "snapshotAnaglyphChannelBytes", 
function () {
});
Clazz.defineMethod (c$, "getScreenImage", 
function (isImageWrite) {
return null;
}, "~B");
Clazz.defineMethod (c$, "releaseScreenImage", 
function () {
});
Clazz.defineMethod (c$, "applyAnaglygh", 
function (stereoMode, stereoColors) {
}, "J.c.STER,~A");
Clazz.defineMethod (c$, "setPass2", 
function (antialias) {
return false;
}, "~B");
Clazz.defineMethod (c$, "destroy", 
function () {
});
Clazz.defineMethod (c$, "clearFontCache", 
function () {
});
Clazz.defineMethod (c$, "plotImage", 
function (x, y, z, image, jmolRenderer, bgcolix, width, height) {
}, "~N,~N,~N,~O,J.api.JmolRendererInterface,~N,~N,~N");
Clazz.defineMethod (c$, "plotText", 
function (x, y, z, colorArgbOrGray, bgColor, text, font3d, jmolRenderer) {
}, "~N,~N,~N,~N,~N,~S,javajs.awt.Font,J.api.JmolRendererInterface");
Clazz.defineMethod (c$, "renderBackground", 
function (jmolRenderer) {
}, "J.api.JmolRendererInterface");
Clazz.defineMethod (c$, "getFont3DCurrent", 
function () {
return null;
});
Clazz.defineMethod (c$, "setFont", 
function (font3d) {
}, "javajs.awt.Font");
Clazz.defineMethod (c$, "setFontFid", 
function (fid) {
}, "~N");
Clazz.defineMethod (c$, "setColor", 
function (color) {
this.argbCurrent = color;
}, "~N");
Clazz.defineMethod (c$, "isPass2", 
function () {
return this.$isPass2;
});
Clazz.defineMethod (c$, "setC", 
function (colix) {
return true;
}, "~N");
Clazz.defineMethod (c$, "isDirectedTowardsCamera", 
function (normix) {
return true;
}, "~N");
Clazz.defineMethod (c$, "getTransformedVertexVectors", 
function () {
return null;
});
Clazz.defineMethod (c$, "setNoisySurfaceShade", 
function (pointA, pointB, pointC) {
}, "JU.P3i,JU.P3i,JU.P3i");
c$.roundInt = Clazz.defineMethod (c$, "roundInt", 
function (a) {
{
return a;
}}, "~N");
Clazz.defineMethod (c$, "clear", 
function () {
});
Clazz.overrideMethod (c$, "renderAllStrings", 
function (jmolRenderer) {
}, "~O");
c$.getScreenOctant = Clazz.defineMethod (c$, "getScreenOctant", 
function (pt) {
var i = 0;
if (pt.x < 0) i |= 1;
if (pt.y < 0) i |= 2;
if (pt.z < 0) i |= 4;
return i;
}, "JU.P3");
Clazz.defineMethod (c$, "addRenderer", 
function (tok) {
}, "~N");
c$.getHermiteList = Clazz.defineMethod (c$, "getHermiteList", 
function (tension, p0, p1, p2, p3, p4, list, index0, n, isPt) {
var nPoints = n + 1;
var fnPoints = n - 1;
var x1 = p1.x;
var y1 = p1.y;
var z1 = p1.z;
var x2 = p2.x;
var y2 = p2.y;
var z2 = p2.z;
var xT1 = ((x2 - p0.x) * tension) / 8;
var yT1 = ((y2 - p0.y) * tension) / 8;
var zT1 = ((z2 - p0.z) * tension) / 8;
var xT2 = ((p3.x - x1) * tension) / 8;
var yT2 = ((p3.y - y1) * tension) / 8;
var zT2 = ((p3.z - z1) * tension) / 8;
var xT3 = ((p4.x - x2) * tension) / 8;
var yT3 = ((p4.y - y2) * tension) / 8;
var zT3 = ((p4.z - z2) * tension) / 8;
list[index0] = p1;
for (var i = 0; i < nPoints; i++) {
var s = i / fnPoints;
if (i == nPoints - 1) {
x1 = x2;
y1 = y2;
z1 = z2;
x2 = p3.x;
y2 = p3.y;
z2 = p3.z;
xT1 = xT2;
yT1 = yT2;
zT1 = zT2;
xT2 = xT3;
yT2 = yT3;
zT2 = zT3;
s -= 1;
}var s2 = s * s;
var s3 = s2 * s;
var h1 = 2 * s3 - 3 * s2 + 1;
var h2 = -2 * s3 + 3 * s2;
var h3 = s3 - 2 * s2 + s;
var h4 = s3 - s2;
var x = (h1 * x1 + h2 * x2 + h3 * xT1 + h4 * xT2);
var y = (h1 * y1 + h2 * y2 + h3 * yT1 + h4 * yT2);
var z = (h1 * z1 + h2 * z2 + h3 * zT1 + h4 * zT2);
list[index0 + i] = (isPt ? JU.P3.new3 (x, y, z) : JU.V3.new3 (x, y, z));
}
}, "~N,JU.T3,JU.T3,JU.T3,JU.T3,JU.T3,~A,~N,~N,~B");
Clazz.defineStatics (c$,
"ENDCAPS_NONE", 0,
"ENDCAPS_OPEN", 1,
"ENDCAPS_FLAT", 2,
"ENDCAPS_SPHERICAL", 3,
"ENDCAPS_OPENEND", 4,
"EXPORT_RAYTRACER", 2,
"EXPORT_CARTESIAN", 1,
"EXPORT_NOT", 0,
"yGT", 1,
"yLT", 2,
"xGT", 4,
"xLT", 8,
"zGT", 16,
"zLT", 32);
});
