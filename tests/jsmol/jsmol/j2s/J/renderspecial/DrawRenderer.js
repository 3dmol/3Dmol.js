Clazz.declarePackage ("J.renderspecial");
Clazz.load (["J.render.MeshRenderer", "JU.BS", "$.P3", "$.P3i", "$.V3"], "J.renderspecial.DrawRenderer", ["JU.A4", "$.M3", "J.shapespecial.Draw", "JU.C", "$.GData", "$.Measure"], function () {
c$ = Clazz.decorateAsClass (function () {
this.drawType = null;
this.dmesh = null;
this.controlHermites = null;
this.pt0 = null;
this.pt1 = null;
this.pt2 = null;
this.vTemp = null;
this.vTemp2 = null;
this.pt0f = null;
this.pt0i = null;
this.bsHandles = null;
Clazz.instantialize (this, arguments);
}, J.renderspecial, "DrawRenderer", J.render.MeshRenderer);
Clazz.prepareFields (c$, function () {
this.pt0 =  new JU.P3 ();
this.pt1 =  new JU.P3 ();
this.pt2 =  new JU.P3 ();
this.vTemp =  new JU.V3 ();
this.vTemp2 =  new JU.V3 ();
this.pt0f =  new JU.P3 ();
this.pt0i =  new JU.P3i ();
this.bsHandles =  new JU.BS ();
});
Clazz.overrideMethod (c$, "render", 
function () {
this.needTranslucent = false;
this.imageFontScaling = this.vwr.getImageFontScaling ();
var draw = this.shape;
for (var i = draw.meshCount; --i >= 0; ) if (this.renderMesh (this.dmesh = draw.meshes[i])) this.renderInfo ();

return this.needTranslucent;
});
Clazz.overrideMethod (c$, "isPolygonDisplayable", 
function (i) {
return J.shapespecial.Draw.isPolygonDisplayable (this.dmesh, i) && (this.dmesh.modelFlags == null || this.dmesh.bsMeshesVisible.get (i));
}, "~N");
Clazz.overrideMethod (c$, "renderMesh", 
function (mesh) {
if (mesh.connections != null) {
if (mesh.connections[0] < 0) return false;
mesh.vs =  new Array (4);
mesh.vc = 4;
var c = mesh.connections;
for (var i = 0; i < 4; i++) {
mesh.vs[i] = (c[i] < 0 ? mesh.vs[i - 1] : this.vwr.getAtomPoint3f (c[i]));
}
mesh.recalcAltVertices = true;
}return this.renderMesh2 (mesh);
}, "J.shape.Mesh");
Clazz.overrideMethod (c$, "render2", 
function (isExport) {
this.drawType = this.dmesh.drawType;
this.diameter = this.dmesh.diameter;
this.width = this.dmesh.width;
if (this.mesh.connections != null) this.getConnectionPoints ();
if (this.mesh.lineData != null) {
this.drawLineData (this.mesh.lineData);
return;
}var isDrawPickMode = (this.vwr.getPickingMode () == 4);
var nPoints = this.vertexCount;
var isCurved = ((this.drawType === J.shapespecial.Draw.EnumDrawType.CURVE || this.drawType === J.shapespecial.Draw.EnumDrawType.ARROW || this.drawType === J.shapespecial.Draw.EnumDrawType.ARC) && this.vertexCount >= 2);
var isSegments = (this.drawType === J.shapespecial.Draw.EnumDrawType.LINE_SEGMENT);
if (this.width > 0 && isCurved) {
this.pt1f.set (0, 0, 0);
var n = (this.drawType === J.shapespecial.Draw.EnumDrawType.ARC ? 2 : this.vertexCount);
for (var i = 0; i < n; i++) this.pt1f.add (this.vertices[i]);

this.pt1f.scale (1 / n);
this.tm.transformPtScr (this.pt1f, this.pt1i);
this.diameter = Clazz.floatToInt (this.vwr.tm.scaleToScreen (this.pt1i.z, Clazz.doubleToInt (Math.floor (this.width * 1000))));
if (this.diameter == 0) this.diameter = 1;
}if ((this.dmesh.isVector) && this.dmesh.haveXyPoints) {
var ptXY = 0;
for (var i = 0; i < 2; i++) if (this.vertices[i].z == 3.4028235E38 || this.vertices[i].z == -3.4028235E38) ptXY += i + 1;

if (--ptXY < 2) {
this.renderXyArrow (ptXY);
return;
}}var tension = 5;
switch (this.drawType) {
default:
this.render2b (false);
break;
case J.shapespecial.Draw.EnumDrawType.CIRCULARPLANE:
if (this.dmesh.scale > 0) this.width *= this.dmesh.scale;
this.render2b (false);
break;
case J.shapespecial.Draw.EnumDrawType.CIRCLE:
this.tm.transformPtScr (this.vertices[0], this.pt1i);
if (this.diameter == 0 && this.width == 0) this.width = 1.0;
if (this.dmesh.scale > 0) this.width *= this.dmesh.scale;
if (this.width > 0) this.diameter = Clazz.floatToInt (this.vwr.tm.scaleToScreen (this.pt1i.z, Clazz.doubleToInt (Math.floor (this.width * 1000))));
if (this.diameter > 0 && (this.mesh.drawTriangles || this.mesh.fillTriangles)) {
this.g3d.addRenderer (1073741880);
this.g3d.drawFilledCircle (this.colix, this.mesh.fillTriangles ? this.colix : 0, this.diameter, this.pt1i.x, this.pt1i.y, this.pt1i.z);
}break;
case J.shapespecial.Draw.EnumDrawType.CURVE:
case J.shapespecial.Draw.EnumDrawType.LINE_SEGMENT:
break;
case J.shapespecial.Draw.EnumDrawType.ARC:
var nDegreesOffset = (this.vertexCount > 3 ? this.vertices[3].x : 0);
var theta = (this.vertexCount > 3 ? this.vertices[3].y : 360);
if (theta == 0) return;
var fractionalOffset = (this.vertexCount > 3 ? this.vertices[3].z : 0);
this.vTemp.sub2 (this.vertices[1], this.vertices[0]);
this.pt1f.scaleAdd2 (fractionalOffset, this.vTemp, this.vertices[0]);
var mat =  new JU.M3 ().setAA (JU.A4.newVA (this.vTemp, (nDegreesOffset * 3.141592653589793 / 180)));
this.vTemp2.sub2 (this.vertexCount > 2 ? this.vertices[2] : J.shapespecial.Draw.randomPoint (), this.vertices[0]);
this.vTemp2.cross (this.vTemp, this.vTemp2);
this.vTemp2.cross (this.vTemp2, this.vTemp);
this.vTemp2.normalize ();
this.vTemp2.scale (this.dmesh.scale / 2);
mat.rotate (this.vTemp2);
var degrees = theta / 5;
while (Math.abs (degrees) > 5) degrees /= 2;

nPoints = Math.round (theta / degrees) + 1;
while (nPoints < 10) {
degrees /= 2;
nPoints = Math.round (theta / degrees) + 1;
}
mat.setAA (JU.A4.newVA (this.vTemp, (degrees * 3.141592653589793 / 180)));
this.screens = this.vwr.allocTempScreens (nPoints);
var iBase = nPoints - (this.dmesh.scale < 2 ? 3 : 3);
for (var i = 0; i < nPoints; i++) {
if (i == iBase) this.pt0.setT (this.pt1);
this.pt1.scaleAdd2 (1, this.vTemp2, this.pt1f);
if (i == 0) this.pt2.setT (this.pt1);
this.tm.transformPtScr (this.pt1, this.screens[i]);
mat.rotate (this.vTemp2);
}
if (this.dmesh.isVector && !this.dmesh.noHead) {
this.renderArrowHead (this.pt0, this.pt1, 0.3, false, false, this.dmesh.isBarb);
this.tm.transformPtScr (this.pt1f, this.screens[nPoints - 1]);
}this.pt1f.setT (this.pt2);
break;
case J.shapespecial.Draw.EnumDrawType.ARROW:
if (this.vertexCount == 2) {
this.renderArrowHead (this.vertices[0], this.vertices[1], 0, false, true, this.dmesh.isBarb);
break;
}var nHermites = 5;
if (this.controlHermites == null || this.controlHermites.length < nHermites + 1) {
this.controlHermites =  new Array (nHermites + 1);
}JU.GData.getHermiteList (tension, this.vertices[this.vertexCount - 3], this.vertices[this.vertexCount - 2], this.vertices[this.vertexCount - 1], this.vertices[this.vertexCount - 1], this.vertices[this.vertexCount - 1], this.controlHermites, 0, nHermites, true);
this.renderArrowHead (this.controlHermites[nHermites - 2], this.controlHermites[nHermites - 1], 0, false, false, this.dmesh.isBarb);
break;
}
if (this.diameter == 0) this.diameter = 3;
if (isCurved) {
this.g3d.addRenderer (553648147);
for (var i = 0, i0 = 0; i < nPoints - 1; i++) {
this.g3d.fillHermite (tension, this.diameter, this.diameter, this.diameter, this.screens[i0], this.screens[i], this.screens[i + 1], this.screens[i + (i == nPoints - 2 ? 1 : 2)]);
i0 = i;
}
} else if (isSegments) {
for (var i = 0; i < nPoints - 1; i++) this.drawLine (i, i + 1, true, this.vertices[i], this.vertices[i + 1], this.screens[i], this.screens[i + 1]);

}if (isDrawPickMode && !isExport) {
this.renderHandles ();
}}, "~B");
Clazz.defineMethod (c$, "getConnectionPoints", 
 function () {
this.vertexCount = 3;
var dmax = 3.4028235E38;
var i0 = 0;
var j0 = 0;
for (var i = 0; i < 2; i++) for (var j = 2; j < 4; j++) {
var d = this.vertices[i].distance (this.vertices[j]);
if (d < dmax) {
dmax = d;
i0 = i;
j0 = j;
}}

this.pt0.ave (this.vertices[0], this.vertices[1]);
this.pt2.ave (this.vertices[2], this.vertices[3]);
this.pt1.ave (this.pt0, this.pt2);
this.vertices[3] = JU.P3.newP (this.vertices[i0]);
this.vertices[3].add (this.vertices[j0]);
this.vertices[3].scale (0.5);
this.vertices[1] = JU.P3.newP (this.pt1);
this.vertices[0] = JU.P3.newP (this.pt0);
this.vertices[2] = JU.P3.newP (this.pt2);
for (var i = 0; i < 4; i++) this.tm.transformPtScr (this.vertices[i], this.screens[i]);

var f = 4 * this.getArrowScale ();
var endoffset = 0.2;
var offsetside = (this.width == 0 ? 0.1 : this.width);
this.pt0.set (this.screens[0].x, this.screens[0].y, this.screens[0].z);
this.pt1.set (this.screens[1].x, this.screens[1].y, this.screens[1].z);
this.pt2.set (this.screens[3].x, this.screens[3].y, this.screens[3].z);
var dx = (this.screens[1].x - this.screens[0].x) * f;
var dy = (this.screens[1].y - this.screens[0].y) * f;
if (dmax == 0 || JU.Measure.computeTorsion (this.pt2, this.pt0, JU.P3.new3 (this.pt0.x, this.pt0.y, 10000), this.pt1, false) > 0) {
dx = -dx;
dy = -dy;
}this.pt2.set (dy, -dx, 0);
this.pt1.add (this.pt2);
this.tm.unTransformPoint (this.pt1, this.vertices[1]);
this.pt2.scale (offsetside);
this.vTemp.sub2 (this.vertices[1], this.vertices[0]);
this.vTemp.scale (endoffset);
this.vertices[0].add (this.vTemp);
this.vTemp.sub2 (this.vertices[1], this.vertices[2]);
this.vTemp.scale (endoffset);
this.vertices[2].add (this.vTemp);
for (var i = 0; i < 3; i++) {
this.tm.transformPtScr (this.vertices[i], this.screens[i]);
if (offsetside != 0) {
this.screens[i].x += Math.round (this.pt2.x);
this.screens[i].y += Math.round (this.pt2.y);
this.pt1.set (this.screens[i].x, this.screens[i].y, this.screens[i].z);
this.tm.unTransformPoint (this.pt1, this.vertices[i]);
}}
});
Clazz.defineMethod (c$, "drawLineData", 
 function (lineData) {
if (this.diameter == 0) this.diameter = 3;
for (var i = lineData.size (); --i >= 0; ) {
var pts = lineData.get (i);
this.tm.transformPtScr (pts[0], this.pt1i);
this.tm.transformPtScr (pts[1], this.pt2i);
this.drawLine (-1, -2, true, pts[0], pts[1], this.pt1i, this.pt2i);
}
}, "JU.Lst");
Clazz.defineMethod (c$, "renderXyArrow", 
 function (ptXY) {
var ptXYZ = 1 - ptXY;
var arrowPt =  new Array (2);
arrowPt[ptXYZ] = this.pt1;
arrowPt[ptXY] = this.pt0;
this.pt0.set (this.screens[ptXY].x, this.screens[ptXY].y, this.screens[ptXY].z);
this.tm.rotatePoint (this.vertices[ptXYZ], this.pt1);
this.pt1.z *= -1;
var zoomDimension = this.vwr.getScreenDim ();
var scaleFactor = zoomDimension / 20;
this.pt1.scaleAdd2 (this.dmesh.scale * scaleFactor, this.pt1, this.pt0);
if (this.diameter == 0) this.diameter = 1;
this.pt1i.set (Math.round (this.pt0.x), Math.round (this.pt0.y), Math.round (this.pt0.z));
this.pt2i.set (Math.round (this.pt1.x), Math.round (this.pt1.y), Math.round (this.pt1.z));
if (this.diameter < 0) this.g3d.drawDottedLine (this.pt1i, this.pt2i);
 else this.g3d.fillCylinder (2, this.diameter, this.pt1i, this.pt2i);
this.renderArrowHead (this.pt0, this.pt1, 0, true, false, false);
}, "~N");
Clazz.defineMethod (c$, "renderArrowHead", 
 function (pt1, pt2, factor2, isTransformed, withShaft, isBarb) {
if (this.dmesh.noHead) return;
var fScale = this.getArrowScale ();
if (isTransformed) fScale *= 40;
if (factor2 > 0) fScale *= factor2;
this.pt0f.setT (pt1);
this.pt2f.setT (pt2);
var d = this.pt0f.distance (this.pt2f);
if (d == 0) return;
this.vTemp.sub2 (this.pt2f, this.pt0f);
this.vTemp.normalize ();
this.vTemp.scale (fScale / 5);
if (!withShaft) this.pt2f.add (this.vTemp);
this.vTemp.scale (5);
this.pt1f.sub2 (this.pt2f, this.vTemp);
if (isTransformed) {
this.pt1i.set (Math.round (this.pt1f.x), Math.round (this.pt1f.y), Math.round (this.pt1f.z));
this.pt2i.set (Math.round (this.pt2f.x), Math.round (this.pt2f.y), Math.round (this.pt2f.z));
} else {
this.tm.transformPtScr (this.pt2f, this.pt2i);
this.tm.transformPtScr (this.pt1f, this.pt1i);
this.tm.transformPtScr (this.pt0f, this.pt0i);
}if (this.pt2i.z == 1 || this.pt1i.z == 1) return;
var headDiameter;
if (this.diameter > 0) {
headDiameter = this.diameter * 3;
} else {
this.vTemp.set (this.pt2i.x - this.pt1i.x, this.pt2i.y - this.pt1i.y, this.pt2i.z - this.pt1i.z);
headDiameter = Math.round (this.vTemp.length () * .5);
this.diameter = Clazz.doubleToInt (headDiameter / 5);
}if (this.diameter < 1) this.diameter = 1;
if (headDiameter > 2) this.g3d.fillConeScreen (2, headDiameter, this.pt1i, this.pt2i, isBarb);
if (withShaft) this.g3d.fillCylinderScreen3I (4, this.diameter, this.pt0i, this.pt1i, null, null, 0);
}, "JU.T3,JU.T3,~N,~B,~B,~B");
Clazz.defineMethod (c$, "getArrowScale", 
 function () {
var fScale = (this.dmesh.isScaleSet ? this.dmesh.scale : 0);
if (fScale == 0) fScale = this.vwr.getFloat (570425352) * (this.dmesh.connections == null ? 1 : 0.5);
if (fScale <= 0) fScale = 0.5;
return fScale;
});
Clazz.defineMethod (c$, "renderHandles", 
 function () {
var diameter = Math.round (10 * this.imageFontScaling);
switch (this.drawType) {
case J.shapespecial.Draw.EnumDrawType.NONE:
return;
default:
var colixFill = JU.C.getColixTranslucent3 (23, true, 0.5);
this.bsHandles.clearAll ();
this.g3d.addRenderer (1073741880);
for (var i = this.dmesh.pc; --i >= 0; ) {
if (!this.isPolygonDisplayable (i)) continue;
var vertexIndexes = this.dmesh.pis[i];
if (vertexIndexes == null) continue;
for (var j = (this.dmesh.isTriangleSet ? 3 : vertexIndexes.length); --j >= 0; ) {
var k = vertexIndexes[j];
if (this.bsHandles.get (k)) continue;
this.bsHandles.set (k);
this.g3d.drawFilledCircle (23, colixFill, diameter, this.screens[k].x, this.screens[k].y, this.screens[k].z);
}
}
break;
}
});
Clazz.defineMethod (c$, "renderInfo", 
 function () {
if (this.mesh.title == null || this.vwr.getDrawHover () || !this.g3d.setC (this.vwr.getColixBackgroundContrast ())) return;
for (var i = this.dmesh.pc; --i >= 0; ) if (this.isPolygonDisplayable (i)) {
var size = this.vwr.getFloat (570425356);
if (size <= 0) size = 14;
var fid = this.g3d.getFontFid (size * this.imageFontScaling);
this.g3d.setFontFid (fid);
var s = this.mesh.title[i < this.mesh.title.length ? i : this.mesh.title.length - 1];
var pt = 0;
if (s.length > 1 && s.charAt (0) == '>') {
pt = this.dmesh.pis[i].length - 1;
s = s.substring (1);
if (this.drawType === J.shapespecial.Draw.EnumDrawType.ARC) this.pt1f.setT (this.pt2f);
}if (this.drawType !== J.shapespecial.Draw.EnumDrawType.ARC) this.pt1f.setT (this.vertices[this.dmesh.pis[i][pt]]);
this.tm.transformPtScr (this.pt1f, this.pt1i);
var offset = Math.round (5 * this.imageFontScaling);
this.g3d.drawString (s, null, this.pt1i.x + offset, this.pt1i.y - offset, this.pt1i.z, this.pt1i.z, 0);
break;
}
});
});
