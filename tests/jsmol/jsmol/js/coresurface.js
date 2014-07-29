(function(Clazz
,Clazz_doubleToInt
,Clazz_declarePackage
,Clazz_instanceOf
,Clazz_load
,Clazz_instantialize
,Clazz_decorateAsClass
,Clazz_floatToInt
,Clazz_makeConstructor
,Clazz_defineEnumConstant
,Clazz_exceptionOf
,Clazz_newIntArray
,Clazz_defineStatics
,Clazz_newFloatArray
,Clazz_declareType
,Clazz_prepareFields
,Clazz_superConstructor
,Clazz_newByteArray
,Clazz_declareInterface
,Clazz_p0p
,Clazz_pu$h
,Clazz_newShortArray
,Clazz_innerTypeInstance
,Clazz_isClassDefined
,Clazz_prepareCallback
,Clazz_newArray
,Clazz_castNullAs
,Clazz_floatToShort
,Clazz_superCall
,Clazz_decorateAsType
,Clazz_newBooleanArray
,Clazz_newCharArray
,Clazz_implementOf
,Clazz_newDoubleArray
,Clazz_overrideConstructor
,Clazz_clone
,Clazz_doubleToShort
,Clazz_getInheritedLevel
,Clazz_getParamsType
,Clazz_isAF
,Clazz_isAI
,Clazz_isAS
,Clazz_isASS
,Clazz_isAP
,Clazz_isAFloat
,Clazz_isAII
,Clazz_isAFF
,Clazz_isAFFF
,Clazz_tryToSearchAndExecute
,Clazz_getStackTrace
,Clazz_inheritArgs
,Clazz_alert
,Clazz_defineMethod
,Clazz_overrideMethod
,Clazz_declareAnonymous
//,Clazz_checkPrivateMethod
,Clazz_cloneFinals
){
var $t$;
//var c$;
Clazz_declarePackage ("JU");
Clazz_declareInterface (JU, "Triangulator");
Clazz_declarePackage ("JU");
Clazz_load (["JU.Triangulator", "JU.P3i"], "JU.TriangleData", ["JU.AU", "$.BS", "$.P3", "JU.BSUtil"], function () {
c$ = Clazz_declareType (JU, "TriangleData", null, JU.Triangulator);
Clazz_overrideMethod (c$, "intersectPlane", 
function (plane, v, flags) {
if (plane == null) {
v.addLast (JU.TriangleData.fullCubePolygon);
return v;
}var vertices = v.get (0);
if (flags != 0) v.clear ();
var values =  Clazz_newFloatArray (8, 0);
var edgePoints =  new Array (12);
var insideMask = 0;
for (var i = 0; i < 8; i++) {
values[i] = plane.x * vertices[i].x + plane.y * vertices[i].y + plane.z * vertices[i].z + plane.w;
if (values[i] < 0) insideMask |= JU.TriangleData.Pwr2[i];
}
var triangles = JU.TriangleData.triangleTable2[insideMask];
if (triangles == null) return null;
for (var i = 0; i < 24; i += 2) {
var v1 = JU.TriangleData.edgeVertexes[i];
var v2 = JU.TriangleData.edgeVertexes[i + 1];
var result = JU.P3.newP (vertices[v2]);
result.sub (vertices[v1]);
result.scale (values[v1] / (values[v1] - values[v2]));
result.add (vertices[v1]);
edgePoints[i >> 1] = result;
}
if (flags == 0) {
var bsPoints =  new JU.BS ();
v.clear ();
for (var i = 0; i < triangles.length; i++) {
bsPoints.set (triangles[i]);
if (i % 4 == 2) i++;
}
var nPoints = JU.BSUtil.cardinalityOf (bsPoints);
var pts =  new Array (nPoints);
v.addLast (pts);
var list =  Clazz_newIntArray (12, 0);
var ptList = 0;
for (var i = 0; i < triangles.length; i++) {
var pt = triangles[i];
if (bsPoints.get (pt)) {
bsPoints.clear (pt);
pts[ptList] = edgePoints[pt];
list[pt] = ptList++;
}if (i % 4 == 2) i++;
}
var polygons = JU.AU.newInt2 (triangles.length >> 2);
v.addLast (polygons);
for (var i = 0; i < triangles.length; i++) polygons[i >> 2] = [list[triangles[i++]], list[triangles[i++]], list[triangles[i++]], triangles[i]];

return v;
}for (var i = 0; i < triangles.length; i++) {
var pt1 = edgePoints[triangles[i++]];
var pt2 = edgePoints[triangles[i++]];
var pt3 = edgePoints[triangles[i++]];
if ((flags & 1) == 1) v.addLast ([pt1, pt2, pt3]);
if ((flags & 2) == 2) {
var b = triangles[i];
if ((b & 1) == 1) v.addLast ([pt1, pt2]);
if ((b & 2) == 2) v.addLast ([pt2, pt3]);
if ((b & 4) == 4) v.addLast ([pt1, pt3]);
}}
return v;
}, "JU.P4,JU.Lst,~N");
Clazz_defineStatics (c$,
"Pwr2", [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048],
"fullCubePolygon", [[0, 4, 5, 3], [5, 1, 0, 3], [1, 5, 6, 2], [6, 2, 1, 3], [2, 6, 7, 2], [7, 3, 2, 3], [3, 7, 4, 2], [4, 0, 3, 2], [6, 5, 4, 0], [4, 7, 6, 0], [0, 1, 2, 0], [2, 3, 0, 0]]);
c$.cubeVertexOffsets = c$.prototype.cubeVertexOffsets = [JU.P3i.new3 (0, 0, 0), JU.P3i.new3 (1, 0, 0), JU.P3i.new3 (1, 0, 1), JU.P3i.new3 (0, 0, 1), JU.P3i.new3 (0, 1, 0), JU.P3i.new3 (1, 1, 0), JU.P3i.new3 (1, 1, 1), JU.P3i.new3 (0, 1, 1)];
Clazz_defineStatics (c$,
"edgeVertexes", [0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6, 6, 7, 7, 4, 0, 4, 1, 5, 2, 6, 3, 7],
"triangleTable2", [null, [0, 8, 3, 7], [0, 1, 9, 7], [1, 8, 3, 6, 9, 8, 1, 5], [1, 2, 10, 7], [0, 8, 3, 7, 1, 2, 10, 7], [9, 2, 10, 6, 0, 2, 9, 5], [2, 8, 3, 6, 2, 10, 8, 1, 10, 9, 8, 3], [3, 11, 2, 7], [0, 11, 2, 6, 8, 11, 0, 5], [1, 9, 0, 7, 2, 3, 11, 7], [1, 11, 2, 6, 1, 9, 11, 1, 9, 8, 11, 3], [3, 10, 1, 6, 11, 10, 3, 5], [0, 10, 1, 6, 0, 8, 10, 1, 8, 11, 10, 3], [3, 9, 0, 6, 3, 11, 9, 1, 11, 10, 9, 3], [9, 8, 10, 5, 10, 8, 11, 6], [4, 7, 8, 7], [4, 3, 0, 6, 7, 3, 4, 5], [0, 1, 9, 7, 8, 4, 7, 7], [4, 1, 9, 6, 4, 7, 1, 1, 7, 3, 1, 3], [1, 2, 10, 7, 8, 4, 7, 7], [3, 4, 7, 6, 3, 0, 4, 3, 1, 2, 10, 7], [9, 2, 10, 6, 9, 0, 2, 3, 8, 4, 7, 7], [2, 10, 9, 3, 2, 9, 7, 0, 2, 7, 3, 6, 7, 9, 4, 6], [8, 4, 7, 7, 3, 11, 2, 7], [11, 4, 7, 6, 11, 2, 4, 1, 2, 0, 4, 3], [9, 0, 1, 7, 8, 4, 7, 7, 2, 3, 11, 7], [4, 7, 11, 3, 9, 4, 11, 1, 9, 11, 2, 2, 9, 2, 1, 6], [3, 10, 1, 6, 3, 11, 10, 3, 7, 8, 4, 7], [1, 11, 10, 6, 1, 4, 11, 0, 1, 0, 4, 3, 7, 11, 4, 5], [4, 7, 8, 7, 9, 0, 11, 1, 9, 11, 10, 6, 11, 0, 3, 6], [4, 7, 11, 3, 4, 11, 9, 4, 9, 11, 10, 6], [9, 5, 4, 7], [9, 5, 4, 7, 0, 8, 3, 7], [0, 5, 4, 6, 1, 5, 0, 5], [8, 5, 4, 6, 8, 3, 5, 1, 3, 1, 5, 3], [1, 2, 10, 7, 9, 5, 4, 7], [3, 0, 8, 7, 1, 2, 10, 7, 4, 9, 5, 7], [5, 2, 10, 6, 5, 4, 2, 1, 4, 0, 2, 3], [2, 10, 5, 3, 3, 2, 5, 1, 3, 5, 4, 2, 3, 4, 8, 6], [9, 5, 4, 7, 2, 3, 11, 7], [0, 11, 2, 6, 0, 8, 11, 3, 4, 9, 5, 7], [0, 5, 4, 6, 0, 1, 5, 3, 2, 3, 11, 7], [2, 1, 5, 3, 2, 5, 8, 0, 2, 8, 11, 6, 4, 8, 5, 5], [10, 3, 11, 6, 10, 1, 3, 3, 9, 5, 4, 7], [4, 9, 5, 7, 0, 8, 1, 5, 8, 10, 1, 2, 8, 11, 10, 3], [5, 4, 0, 3, 5, 0, 11, 0, 5, 11, 10, 6, 11, 0, 3, 6], [5, 4, 8, 3, 5, 8, 10, 4, 10, 8, 11, 6], [9, 7, 8, 6, 5, 7, 9, 5], [9, 3, 0, 6, 9, 5, 3, 1, 5, 7, 3, 3], [0, 7, 8, 6, 0, 1, 7, 1, 1, 5, 7, 3], [1, 5, 3, 5, 3, 5, 7, 6], [9, 7, 8, 6, 9, 5, 7, 3, 10, 1, 2, 7], [10, 1, 2, 7, 9, 5, 0, 5, 5, 3, 0, 2, 5, 7, 3, 3], [8, 0, 2, 3, 8, 2, 5, 0, 8, 5, 7, 6, 10, 5, 2, 5], [2, 10, 5, 3, 2, 5, 3, 4, 3, 5, 7, 6], [7, 9, 5, 6, 7, 8, 9, 3, 3, 11, 2, 7], [9, 5, 7, 3, 9, 7, 2, 0, 9, 2, 0, 6, 2, 7, 11, 6], [2, 3, 11, 7, 0, 1, 8, 5, 1, 7, 8, 2, 1, 5, 7, 3], [11, 2, 1, 3, 11, 1, 7, 4, 7, 1, 5, 6], [9, 5, 8, 5, 8, 5, 7, 6, 10, 1, 3, 3, 10, 3, 11, 6], [5, 7, 0, 1, 5, 0, 9, 6, 7, 11, 0, 1, 1, 0, 10, 5, 11, 10, 0, 1], [11, 10, 0, 1, 11, 0, 3, 6, 10, 5, 0, 1, 8, 0, 7, 5, 5, 7, 0, 1], [11, 10, 5, 3, 7, 11, 5, 5], [10, 6, 5, 7], [0, 8, 3, 7, 5, 10, 6, 7], [9, 0, 1, 7, 5, 10, 6, 7], [1, 8, 3, 6, 1, 9, 8, 3, 5, 10, 6, 7], [1, 6, 5, 6, 2, 6, 1, 5], [1, 6, 5, 6, 1, 2, 6, 3, 3, 0, 8, 7], [9, 6, 5, 6, 9, 0, 6, 1, 0, 2, 6, 3], [5, 9, 8, 3, 5, 8, 2, 0, 5, 2, 6, 6, 3, 2, 8, 5], [2, 3, 11, 7, 10, 6, 5, 7], [11, 0, 8, 6, 11, 2, 0, 3, 10, 6, 5, 7], [0, 1, 9, 7, 2, 3, 11, 7, 5, 10, 6, 7], [5, 10, 6, 7, 1, 9, 2, 5, 9, 11, 2, 2, 9, 8, 11, 3], [6, 3, 11, 6, 6, 5, 3, 1, 5, 1, 3, 3], [0, 8, 11, 3, 0, 11, 5, 0, 0, 5, 1, 6, 5, 11, 6, 6], [3, 11, 6, 3, 0, 3, 6, 1, 0, 6, 5, 2, 0, 5, 9, 6], [6, 5, 9, 3, 6, 9, 11, 4, 11, 9, 8, 6], [5, 10, 6, 7, 4, 7, 8, 7], [4, 3, 0, 6, 4, 7, 3, 3, 6, 5, 10, 7], [1, 9, 0, 7, 5, 10, 6, 7, 8, 4, 7, 7], [10, 6, 5, 7, 1, 9, 7, 1, 1, 7, 3, 6, 7, 9, 4, 6], [6, 1, 2, 6, 6, 5, 1, 3, 4, 7, 8, 7], [1, 2, 5, 5, 5, 2, 6, 6, 3, 0, 4, 3, 3, 4, 7, 6], [8, 4, 7, 7, 9, 0, 5, 5, 0, 6, 5, 2, 0, 2, 6, 3], [7, 3, 9, 1, 7, 9, 4, 6, 3, 2, 9, 1, 5, 9, 6, 5, 2, 6, 9, 1], [3, 11, 2, 7, 7, 8, 4, 7, 10, 6, 5, 7], [5, 10, 6, 7, 4, 7, 2, 1, 4, 2, 0, 6, 2, 7, 11, 6], [0, 1, 9, 7, 4, 7, 8, 7, 2, 3, 11, 7, 5, 10, 6, 7], [9, 2, 1, 6, 9, 11, 2, 2, 9, 4, 11, 1, 7, 11, 4, 5, 5, 10, 6, 7], [8, 4, 7, 7, 3, 11, 5, 1, 3, 5, 1, 6, 5, 11, 6, 6], [5, 1, 11, 1, 5, 11, 6, 6, 1, 0, 11, 1, 7, 11, 4, 5, 0, 4, 11, 1], [0, 5, 9, 6, 0, 6, 5, 2, 0, 3, 6, 1, 11, 6, 3, 5, 8, 4, 7, 7], [6, 5, 9, 3, 6, 9, 11, 4, 4, 7, 9, 5, 7, 11, 9, 1], [10, 4, 9, 6, 6, 4, 10, 5], [4, 10, 6, 6, 4, 9, 10, 3, 0, 8, 3, 7], [10, 0, 1, 6, 10, 6, 0, 1, 6, 4, 0, 3], [8, 3, 1, 3, 8, 1, 6, 0, 8, 6, 4, 6, 6, 1, 10, 6], [1, 4, 9, 6, 1, 2, 4, 1, 2, 6, 4, 3], [3, 0, 8, 7, 1, 2, 9, 5, 2, 4, 9, 2, 2, 6, 4, 3], [0, 2, 4, 5, 4, 2, 6, 6], [8, 3, 2, 3, 8, 2, 4, 4, 4, 2, 6, 6], [10, 4, 9, 6, 10, 6, 4, 3, 11, 2, 3, 7], [0, 8, 2, 5, 2, 8, 11, 6, 4, 9, 10, 3, 4, 10, 6, 6], [3, 11, 2, 7, 0, 1, 6, 1, 0, 6, 4, 6, 6, 1, 10, 6], [6, 4, 1, 1, 6, 1, 10, 6, 4, 8, 1, 1, 2, 1, 11, 5, 8, 11, 1, 1], [9, 6, 4, 6, 9, 3, 6, 0, 9, 1, 3, 3, 11, 6, 3, 5], [8, 11, 1, 1, 8, 1, 0, 6, 11, 6, 1, 1, 9, 1, 4, 5, 6, 4, 1, 1], [3, 11, 6, 3, 3, 6, 0, 4, 0, 6, 4, 6], [6, 4, 8, 3, 11, 6, 8, 5], [7, 10, 6, 6, 7, 8, 10, 1, 8, 9, 10, 3], [0, 7, 3, 6, 0, 10, 7, 0, 0, 9, 10, 3, 6, 7, 10, 5], [10, 6, 7, 3, 1, 10, 7, 1, 1, 7, 8, 2, 1, 8, 0, 6], [10, 6, 7, 3, 10, 7, 1, 4, 1, 7, 3, 6], [1, 2, 6, 3, 1, 6, 8, 0, 1, 8, 9, 6, 8, 6, 7, 6], [2, 6, 9, 1, 2, 9, 1, 6, 6, 7, 9, 1, 0, 9, 3, 5, 7, 3, 9, 1], [7, 8, 0, 3, 7, 0, 6, 4, 6, 0, 2, 6], [7, 3, 2, 3, 6, 7, 2, 5], [2, 3, 11, 7, 10, 6, 8, 1, 10, 8, 9, 6, 8, 6, 7, 6], [2, 0, 7, 1, 2, 7, 11, 6, 0, 9, 7, 1, 6, 7, 10, 5, 9, 10, 7, 1], [1, 8, 0, 6, 1, 7, 8, 2, 1, 10, 7, 1, 6, 7, 10, 5, 2, 3, 11, 7], [11, 2, 1, 3, 11, 1, 7, 4, 10, 6, 1, 5, 6, 7, 1, 1], [8, 9, 6, 1, 8, 6, 7, 6, 9, 1, 6, 1, 11, 6, 3, 5, 1, 3, 6, 1], [0, 9, 1, 7, 11, 6, 7, 7], [7, 8, 0, 3, 7, 0, 6, 4, 3, 11, 0, 5, 11, 6, 0, 1], [7, 11, 6, 7], [7, 6, 11, 7], [3, 0, 8, 7, 11, 7, 6, 7], [0, 1, 9, 7, 11, 7, 6, 7], [8, 1, 9, 6, 8, 3, 1, 3, 11, 7, 6, 7], [10, 1, 2, 7, 6, 11, 7, 7], [1, 2, 10, 7, 3, 0, 8, 7, 6, 11, 7, 7], [2, 9, 0, 6, 2, 10, 9, 3, 6, 11, 7, 7], [6, 11, 7, 7, 2, 10, 3, 5, 10, 8, 3, 2, 10, 9, 8, 3], [7, 2, 3, 6, 6, 2, 7, 5], [7, 0, 8, 6, 7, 6, 0, 1, 6, 2, 0, 3], [2, 7, 6, 6, 2, 3, 7, 3, 0, 1, 9, 7], [1, 6, 2, 6, 1, 8, 6, 0, 1, 9, 8, 3, 8, 7, 6, 3], [10, 7, 6, 6, 10, 1, 7, 1, 1, 3, 7, 3], [10, 7, 6, 6, 1, 7, 10, 4, 1, 8, 7, 2, 1, 0, 8, 3], [0, 3, 7, 3, 0, 7, 10, 0, 0, 10, 9, 6, 6, 10, 7, 5], [7, 6, 10, 3, 7, 10, 8, 4, 8, 10, 9, 6], [6, 8, 4, 6, 11, 8, 6, 5], [3, 6, 11, 6, 3, 0, 6, 1, 0, 4, 6, 3], [8, 6, 11, 6, 8, 4, 6, 3, 9, 0, 1, 7], [9, 4, 6, 3, 9, 6, 3, 0, 9, 3, 1, 6, 11, 3, 6, 5], [6, 8, 4, 6, 6, 11, 8, 3, 2, 10, 1, 7], [1, 2, 10, 7, 3, 0, 11, 5, 0, 6, 11, 2, 0, 4, 6, 3], [4, 11, 8, 6, 4, 6, 11, 3, 0, 2, 9, 5, 2, 10, 9, 3], [10, 9, 3, 1, 10, 3, 2, 6, 9, 4, 3, 1, 11, 3, 6, 5, 4, 6, 3, 1], [8, 2, 3, 6, 8, 4, 2, 1, 4, 6, 2, 3], [0, 4, 2, 5, 4, 6, 2, 3], [1, 9, 0, 7, 2, 3, 4, 1, 2, 4, 6, 6, 4, 3, 8, 6], [1, 9, 4, 3, 1, 4, 2, 4, 2, 4, 6, 6], [8, 1, 3, 6, 8, 6, 1, 0, 8, 4, 6, 3, 6, 10, 1, 3], [10, 1, 0, 3, 10, 0, 6, 4, 6, 0, 4, 6], [4, 6, 3, 1, 4, 3, 8, 6, 6, 10, 3, 1, 0, 3, 9, 5, 10, 9, 3, 1], [10, 9, 4, 3, 6, 10, 4, 5], [4, 9, 5, 7, 7, 6, 11, 7], [0, 8, 3, 7, 4, 9, 5, 7, 11, 7, 6, 7], [5, 0, 1, 6, 5, 4, 0, 3, 7, 6, 11, 7], [11, 7, 6, 7, 8, 3, 4, 5, 3, 5, 4, 2, 3, 1, 5, 3], [9, 5, 4, 7, 10, 1, 2, 7, 7, 6, 11, 7], [6, 11, 7, 7, 1, 2, 10, 7, 0, 8, 3, 7, 4, 9, 5, 7], [7, 6, 11, 7, 5, 4, 10, 5, 4, 2, 10, 2, 4, 0, 2, 3], [3, 4, 8, 6, 3, 5, 4, 2, 3, 2, 5, 1, 10, 5, 2, 5, 11, 7, 6, 7], [7, 2, 3, 6, 7, 6, 2, 3, 5, 4, 9, 7], [9, 5, 4, 7, 0, 8, 6, 1, 0, 6, 2, 6, 6, 8, 7, 6], [3, 6, 2, 6, 3, 7, 6, 3, 1, 5, 0, 5, 5, 4, 0, 3], [6, 2, 8, 1, 6, 8, 7, 6, 2, 1, 8, 1, 4, 8, 5, 5, 1, 5, 8, 1], [9, 5, 4, 7, 10, 1, 6, 5, 1, 7, 6, 2, 1, 3, 7, 3], [1, 6, 10, 6, 1, 7, 6, 2, 1, 0, 7, 1, 8, 7, 0, 5, 9, 5, 4, 7], [4, 0, 10, 1, 4, 10, 5, 6, 0, 3, 10, 1, 6, 10, 7, 5, 3, 7, 10, 1], [7, 6, 10, 3, 7, 10, 8, 4, 5, 4, 10, 5, 4, 8, 10, 1], [6, 9, 5, 6, 6, 11, 9, 1, 11, 8, 9, 3], [3, 6, 11, 6, 0, 6, 3, 4, 0, 5, 6, 2, 0, 9, 5, 3], [0, 11, 8, 6, 0, 5, 11, 0, 0, 1, 5, 3, 5, 6, 11, 3], [6, 11, 3, 3, 6, 3, 5, 4, 5, 3, 1, 6], [1, 2, 10, 7, 9, 5, 11, 1, 9, 11, 8, 6, 11, 5, 6, 6], [0, 11, 3, 6, 0, 6, 11, 2, 0, 9, 6, 1, 5, 6, 9, 5, 1, 2, 10, 7], [11, 8, 5, 1, 11, 5, 6, 6, 8, 0, 5, 1, 10, 5, 2, 5, 0, 2, 5, 1], [6, 11, 3, 3, 6, 3, 5, 4, 2, 10, 3, 5, 10, 5, 3, 1], [5, 8, 9, 6, 5, 2, 8, 0, 5, 6, 2, 3, 3, 8, 2, 5], [9, 5, 6, 3, 9, 6, 0, 4, 0, 6, 2, 6], [1, 5, 8, 1, 1, 8, 0, 6, 5, 6, 8, 1, 3, 8, 2, 5, 6, 2, 8, 1], [1, 5, 6, 3, 2, 1, 6, 5], [1, 3, 6, 1, 1, 6, 10, 6, 3, 8, 6, 1, 5, 6, 9, 5, 8, 9, 6, 1], [10, 1, 0, 3, 10, 0, 6, 4, 9, 5, 0, 5, 5, 6, 0, 1], [0, 3, 8, 7, 5, 6, 10, 7], [10, 5, 6, 7], [11, 5, 10, 6, 7, 5, 11, 5], [11, 5, 10, 6, 11, 7, 5, 3, 8, 3, 0, 7], [5, 11, 7, 6, 5, 10, 11, 3, 1, 9, 0, 7], [10, 7, 5, 6, 10, 11, 7, 3, 9, 8, 1, 5, 8, 3, 1, 3], [11, 1, 2, 6, 11, 7, 1, 1, 7, 5, 1, 3], [0, 8, 3, 7, 1, 2, 7, 1, 1, 7, 5, 6, 7, 2, 11, 6], [9, 7, 5, 6, 9, 2, 7, 0, 9, 0, 2, 3, 2, 11, 7, 3], [7, 5, 2, 1, 7, 2, 11, 6, 5, 9, 2, 1, 3, 2, 8, 5, 9, 8, 2, 1], [2, 5, 10, 6, 2, 3, 5, 1, 3, 7, 5, 3], [8, 2, 0, 6, 8, 5, 2, 0, 8, 7, 5, 3, 10, 2, 5, 5], [9, 0, 1, 7, 5, 10, 3, 1, 5, 3, 7, 6, 3, 10, 2, 6], [9, 8, 2, 1, 9, 2, 1, 6, 8, 7, 2, 1, 10, 2, 5, 5, 7, 5, 2, 1], [1, 3, 5, 5, 3, 7, 5, 3], [0, 8, 7, 3, 0, 7, 1, 4, 1, 7, 5, 6], [9, 0, 3, 3, 9, 3, 5, 4, 5, 3, 7, 6], [9, 8, 7, 3, 5, 9, 7, 5], [5, 8, 4, 6, 5, 10, 8, 1, 10, 11, 8, 3], [5, 0, 4, 6, 5, 11, 0, 0, 5, 10, 11, 3, 11, 3, 0, 3], [0, 1, 9, 7, 8, 4, 10, 1, 8, 10, 11, 6, 10, 4, 5, 6], [10, 11, 4, 1, 10, 4, 5, 6, 11, 3, 4, 1, 9, 4, 1, 5, 3, 1, 4, 1], [2, 5, 1, 6, 2, 8, 5, 0, 2, 11, 8, 3, 4, 5, 8, 5], [0, 4, 11, 1, 0, 11, 3, 6, 4, 5, 11, 1, 2, 11, 1, 5, 5, 1, 11, 1], [0, 2, 5, 1, 0, 5, 9, 6, 2, 11, 5, 1, 4, 5, 8, 5, 11, 8, 5, 1], [9, 4, 5, 7, 2, 11, 3, 7], [2, 5, 10, 6, 3, 5, 2, 4, 3, 4, 5, 2, 3, 8, 4, 3], [5, 10, 2, 3, 5, 2, 4, 4, 4, 2, 0, 6], [3, 10, 2, 6, 3, 5, 10, 2, 3, 8, 5, 1, 4, 5, 8, 5, 0, 1, 9, 7], [5, 10, 2, 3, 5, 2, 4, 4, 1, 9, 2, 5, 9, 4, 2, 1], [8, 4, 5, 3, 8, 5, 3, 4, 3, 5, 1, 6], [0, 4, 5, 3, 1, 0, 5, 5], [8, 4, 5, 3, 8, 5, 3, 4, 9, 0, 5, 5, 0, 3, 5, 1], [9, 4, 5, 7], [4, 11, 7, 6, 4, 9, 11, 1, 9, 10, 11, 3], [0, 8, 3, 7, 4, 9, 7, 5, 9, 11, 7, 2, 9, 10, 11, 3], [1, 10, 11, 3, 1, 11, 4, 0, 1, 4, 0, 6, 7, 4, 11, 5], [3, 1, 4, 1, 3, 4, 8, 6, 1, 10, 4, 1, 7, 4, 11, 5, 10, 11, 4, 1], [4, 11, 7, 6, 9, 11, 4, 4, 9, 2, 11, 2, 9, 1, 2, 3], [9, 7, 4, 6, 9, 11, 7, 2, 9, 1, 11, 1, 2, 11, 1, 5, 0, 8, 3, 7], [11, 7, 4, 3, 11, 4, 2, 4, 2, 4, 0, 6], [11, 7, 4, 3, 11, 4, 2, 4, 8, 3, 4, 5, 3, 2, 4, 1], [2, 9, 10, 6, 2, 7, 9, 0, 2, 3, 7, 3, 7, 4, 9, 3], [9, 10, 7, 1, 9, 7, 4, 6, 10, 2, 7, 1, 8, 7, 0, 5, 2, 0, 7, 1], [3, 7, 10, 1, 3, 10, 2, 6, 7, 4, 10, 1, 1, 10, 0, 5, 4, 0, 10, 1], [1, 10, 2, 7, 8, 7, 4, 7], [4, 9, 1, 3, 4, 1, 7, 4, 7, 1, 3, 6], [4, 9, 1, 3, 4, 1, 7, 4, 0, 8, 1, 5, 8, 7, 1, 1], [4, 0, 3, 3, 7, 4, 3, 5], [4, 8, 7, 7], [9, 10, 8, 5, 10, 11, 8, 3], [3, 0, 9, 3, 3, 9, 11, 4, 11, 9, 10, 6], [0, 1, 10, 3, 0, 10, 8, 4, 8, 10, 11, 6], [3, 1, 10, 3, 11, 3, 10, 5], [1, 2, 11, 3, 1, 11, 9, 4, 9, 11, 8, 6], [3, 0, 9, 3, 3, 9, 11, 4, 1, 2, 9, 5, 2, 11, 9, 1], [0, 2, 11, 3, 8, 0, 11, 5], [3, 2, 11, 7], [2, 3, 8, 3, 2, 8, 10, 4, 10, 8, 9, 6], [9, 10, 2, 3, 0, 9, 2, 5], [2, 3, 8, 3, 2, 8, 10, 4, 0, 1, 8, 5, 1, 10, 8, 1], [1, 10, 2, 7], [1, 3, 8, 3, 9, 1, 8, 5], [0, 9, 1, 7], [0, 3, 8, 7], null]);
});
Clazz_declarePackage ("J.api");
Clazz_declareInterface (J.api, "VolumeDataInterface");
Clazz_declarePackage ("J.jvxl.api");
Clazz_declareInterface (J.jvxl.api, "VertexDataServer");
Clazz_declarePackage ("J.jvxl.api");
Clazz_load (["J.jvxl.api.VertexDataServer"], "J.jvxl.api.MeshDataServer", null, function () {
Clazz_declareInterface (J.jvxl.api, "MeshDataServer", J.jvxl.api.VertexDataServer);
});
Clazz_declarePackage ("J.shapesurface");
Clazz_load (["J.jvxl.api.MeshDataServer", "J.shape.MeshCollection", "JU.P3i", "$.P4"], "J.shapesurface.Isosurface", ["java.io.BufferedReader", "java.lang.Boolean", "$.Float", "java.util.Hashtable", "JU.A4", "$.AU", "$.BS", "$.CU", "$.Lst", "$.M3", "$.P3", "$.PT", "$.Quat", "$.Rdr", "$.SB", "$.V3", "J.jvxl.data.JvxlCoder", "$.JvxlData", "$.MeshData", "J.jvxl.readers.SurfaceGenerator", "J.shape.Mesh", "J.shapesurface.IsosurfaceMesh", "JU.C", "$.Escape", "$.Logger", "$.MeshSurface", "$.Txt", "JV.JC", "$.Viewer"], function () {
c$ = Clazz_decorateAsClass (function () {
this.isomeshes = null;
this.thisMesh = null;
this.actualID = null;
this.iHaveBitSets = false;
this.explicitContours = false;
this.atomIndex = 0;
this.moNumber = 0;
this.moLinearCombination = null;
this.colorType = 0;
this.defaultColix = 0;
this.meshColix = 0;
this.center = null;
this.scale3d = 0;
this.isPhaseColored = false;
this.isColorExplicit = false;
this.scriptAppendix = "";
this.sg = null;
this.jvxlData = null;
this.withinDistance2 = 0;
this.isWithinNot = false;
this.withinPoints = null;
this.cutoffRange = null;
this.allowMesh = true;
this.script = null;
this.iHaveModelIndex = false;
this.nLCAO = 0;
this.lcaoDir = null;
this.associateNormals = false;
this.ptXY = null;
this.keyXy = null;
Clazz_instantialize (this, arguments);
}, J.shapesurface, "Isosurface", J.shape.MeshCollection, J.jvxl.api.MeshDataServer);
Clazz_prepareFields (c$, function () {
this.isomeshes =  new Array (4);
this.lcaoDir =  new JU.P4 ();
this.ptXY =  new JU.P3i ();
});
Clazz_overrideMethod (c$, "allocMesh", 
function (thisID, m) {
var index = this.meshCount++;
this.meshes = this.isomeshes = JU.AU.ensureLength (this.isomeshes, this.meshCount * 2);
this.currentMesh = this.thisMesh = this.isomeshes[index] = (m == null ?  new J.shapesurface.IsosurfaceMesh (thisID, this.colix, index) : m);
this.currentMesh.index = index;
if (this.sg != null) this.sg.setJvxlData (this.jvxlData = this.thisMesh.jvxlData);
}, "~S,J.shape.Mesh");
Clazz_defineMethod (c$, "initShape", 
function () {
Clazz_superCall (this, J.shapesurface.Isosurface, "initShape", []);
this.myType = "isosurface";
this.newSg ();
});
Clazz_defineMethod (c$, "newSg", 
function () {
this.sg =  new J.jvxl.readers.SurfaceGenerator (this.vwr, this, null, this.jvxlData =  new J.jvxl.data.JvxlData ());
this.sg.getParams ().showTiming = this.vwr.getBoolean (603979934);
this.sg.setVersion ("Jmol " + JV.Viewer.getJmolVersion ());
});
Clazz_defineMethod (c$, "clearSg", 
function () {
this.sg = null;
});
Clazz_overrideMethod (c$, "setProperty", 
function (propertyName, value, bs) {
this.setPropI (propertyName, value, bs);
}, "~S,~O,JU.BS");
Clazz_defineMethod (c$, "setPropI", 
function (propertyName, value, bs) {
if ("cache" === propertyName) {
if (this.currentMesh == null) return;
var id = this.currentMesh.thisID;
var imodel = this.currentMesh.modelIndex;
this.vwr.cachePut ("cache://isosurface_" + id, this.getPropI ("jvxlDataXml"));
this.deleteMeshI (this.currentMesh.index);
this.setPropI ("init", null, null);
this.setPropI ("thisID", id, null);
this.setPropI ("modelIndex", Integer.$valueOf (imodel), null);
this.setPropI ("fileName", "cache://isosurface_" + id, null);
this.setPropI ("readFile", null, null);
this.setPropI ("finalize", "isosurface ID " + JU.PT.esc (id) + (imodel >= 0 ? " modelIndex " + imodel : "") + " /*file*/" + JU.PT.esc ("cache://isosurface_" + id), null);
this.setPropI ("clear", null, null);
return;
}if ("delete" === propertyName) {
this.setPropertySuper (propertyName, value, bs);
if (!this.explicitID) this.nLCAO = this.nUnnamed = 0;
this.currentMesh = this.thisMesh = null;
return;
}if ("remapInherited" === propertyName) {
for (var i = this.meshCount; --i >= 0; ) {
if (this.isomeshes[i] != null && "#inherit;".equals (this.isomeshes[i].colorCommand)) this.isomeshes[i].remapColors (this.vwr, null, NaN);
}
return;
}if ("remapColor" === propertyName) {
if (this.thisMesh != null) this.thisMesh.remapColors (this.vwr, value, this.translucentLevel);
return;
}if ("thisID" === propertyName) {
if (this.actualID != null) value = this.actualID;
this.setPropertySuper ("thisID", value, null);
return;
}if ("params" === propertyName) {
if (this.thisMesh != null) {
this.ensureMeshSource ();
this.thisMesh.checkAllocColixes ();
var data = value;
var colixes = data[0];
var atomMap = null;
if (colixes != null) {
for (var i = 0; i < colixes.length; i++) {
var colix = colixes[i];
var f = 0;
if (f > 0.01) colix = JU.C.getColixTranslucent3 (colix, true, f);
colixes[i] = colix;
}
atomMap =  Clazz_newIntArray (bs.length (), 0);
for (var pt = 0, i = bs.nextSetBit (0); i >= 0; i = bs.nextSetBit (i + 1), pt++) atomMap[i] = pt;

}this.thisMesh.setVertexColixesForAtoms (this.vwr, colixes, atomMap, bs);
this.thisMesh.setVertexColorMap ();
}return;
}if ("atomcolor" === propertyName) {
if (this.thisMesh != null) {
this.ensureMeshSource ();
this.thisMesh.colorVertices (JU.C.getColixO (value), bs, true);
}return;
}if ("pointSize" === propertyName) {
if (this.thisMesh != null) {
this.thisMesh.volumeRenderPointSize = (value).floatValue ();
}return;
}if ("vertexcolor" === propertyName) {
if (this.thisMesh != null) {
this.thisMesh.colorVertices (JU.C.getColixO (value), bs, false);
}return;
}if ("colorPhase" === propertyName) {
var colors = value;
var colix0 = JU.C.getColix ((colors[0]).intValue ());
var colix1 = JU.C.getColix ((colors[1]).intValue ());
var id = (this.thisMesh != null ? this.thisMesh.thisID : JU.Txt.isWild (this.previousMeshID) ? this.previousMeshID : null);
var list = this.getMeshList (id, false);
for (var i = list.size (); --i >= 0; ) this.setColorPhase (list.get (i), colix0, colix1);

return;
}if ("color" === propertyName) {
var color = JU.C.getHexCode (JU.C.getColixO (value));
if (this.thisMesh != null) {
this.setIsoMeshColor (this.thisMesh, color);
} else {
var list = this.getMeshList (JU.Txt.isWild (this.previousMeshID) ? this.previousMeshID : null, false);
for (var i = list.size (); --i >= 0; ) this.setIsoMeshColor (list.get (i), color);

}this.setPropertySuper (propertyName, value, bs);
return;
}if ("nocontour" === propertyName) {
if (this.thisMesh != null) {
this.thisMesh.deleteContours ();
}return;
}if ("fixed" === propertyName) {
this.isFixed = (value).booleanValue ();
this.setMeshI ();
return;
}if ("newObject" === propertyName) {
if (this.thisMesh != null) this.thisMesh.clearType (this.thisMesh.meshType, false);
return;
}if ("moveIsosurface" === propertyName) {
if (this.thisMesh != null) {
this.thisMesh.updateCoordinates (value, null);
this.thisMesh.altVertices = null;
}return;
}if ("refreshTrajectories" === propertyName) {
for (var i = this.meshCount; --i >= 0; ) if (this.meshes[i].connections != null && this.meshes[i].modelIndex == ((value)[0]).intValue ()) (this.meshes[i]).updateCoordinates ((value)[2], (value)[1]);

return;
}if ("modelIndex" === propertyName) {
if (!this.iHaveModelIndex) {
this.modelIndex = (value).intValue ();
this.isFixed = (this.modelIndex < 0);
this.sg.setModelIndex (Math.abs (this.modelIndex));
}return;
}if ("lcaoCartoon" === propertyName || "lonePair" === propertyName || "radical" === propertyName) {
var info = value;
if (!this.explicitID) {
this.setPropertySuper ("thisID", null, null);
}if (!this.sg.setParameter ("lcaoCartoonCenter", info[2])) this.drawLcaoCartoon (info[0], info[1], info[3], ("lonePair" === propertyName ? 2 : "radical" === propertyName ? 1 : 0));
return;
}if ("select" === propertyName) {
if (this.iHaveBitSets) return;
}if ("ignore" === propertyName) {
if (this.iHaveBitSets) return;
}if ("meshcolor" === propertyName) {
var rgb = (value).intValue ();
this.meshColix = JU.C.getColix (rgb);
if (this.thisMesh != null) this.thisMesh.meshColix = this.meshColix;
return;
}if ("offset" === propertyName) {
var offset = JU.P3.newP (value);
if (offset.equals (JV.JC.center)) offset = null;
if (this.thisMesh != null) {
this.thisMesh.rotateTranslate (null, offset, true);
this.thisMesh.altVertices = null;
}return;
}if ("rotate" === propertyName) {
var pt4 = value;
if (this.thisMesh != null) {
this.thisMesh.rotateTranslate (JU.Quat.newP4 (pt4), null, true);
this.thisMesh.altVertices = null;
}return;
}if ("bsDisplay" === propertyName) {
this.bsDisplay = value;
return;
}if ("displayWithin" === propertyName) {
var o = value;
this.displayWithinDistance2 = (o[0]).floatValue ();
this.isDisplayWithinNot = (this.displayWithinDistance2 < 0);
this.displayWithinDistance2 *= this.displayWithinDistance2;
this.displayWithinPoints = o[3];
if (this.displayWithinPoints.size () == 0) this.displayWithinPoints = this.vwr.ms.getAtomPointVector (o[2]);
return;
}if ("finalize" === propertyName) {
if (this.thisMesh != null) {
var cmd = value;
if (cmd != null && !cmd.startsWith ("; isosurface map")) {
this.thisMesh.setDiscreteColixes (this.sg.getParams ().contoursDiscrete, this.sg.getParams ().contourColixes);
this.setJvxlInfo ();
}this.setScriptInfo (cmd);
}this.clearSg ();
return;
}if ("connections" === propertyName) {
if (this.currentMesh != null) {
this.connections = value;
if (this.connections[0] >= 0 && this.connections[0] < this.vwr.getAtomCount ()) this.currentMesh.connections = this.connections;
 else this.connections = this.currentMesh.connections = null;
}return;
}if ("cutoffRange" === propertyName) {
this.cutoffRange = value;
return;
}if ("slab" === propertyName) {
if (Clazz_instanceOf (value, Integer)) {
if (this.thisMesh != null) this.thisMesh.jvxlData.slabValue = (value).intValue ();
return;
}if (this.thisMesh != null) {
var slabInfo = value;
var tok = (slabInfo[0]).intValue ();
switch (tok) {
case 1073742018:
var data = slabInfo[1];
var m = this.getMesh (data[1]);
if (m == null) return;
data[1] = m;
break;
}
this.slabPolygons (slabInfo);
return;
}}if ("cap" === propertyName) {
if (this.thisMesh != null && this.thisMesh.pc != 0) {
this.thisMesh.slabPolygons (value, true);
this.thisMesh.initialize (this.thisMesh.lighting, null, null);
return;
}}if ("map" === propertyName) {
if (this.sg != null) this.sg.getParams ().isMapped = true;
this.setProperty ("squareData", Boolean.FALSE, null);
if (this.thisMesh == null || this.thisMesh.vc == 0) return;
}if ("deleteVdw" === propertyName) {
for (var i = this.meshCount; --i >= 0; ) if (this.isomeshes[i].bsVdw != null && (bs == null || bs.intersects (this.isomeshes[i].bsVdw))) this.deleteMeshI (i);

this.currentMesh = this.thisMesh = null;
return;
}if ("mapColor" === propertyName || "readFile" === propertyName) {
if (value == null) {
value = this.vwr.getBufferedReaderOrErrorMessageFromName (this.sg.getFileName (), null, true);
if (Clazz_instanceOf (value, String)) {
JU.Logger.error ("Isosurface: could not open file " + this.sg.getFileName () + " -- " + value);
return;
}if (!(Clazz_instanceOf (value, java.io.BufferedReader))) try {
value = JU.Rdr.getBufferedReader (value, "ISO-8859-1");
} catch (e) {
if (Clazz_exceptionOf (e, java.io.IOException)) {
} else {
throw e;
}
}
}} else if ("atomIndex" === propertyName) {
this.atomIndex = (value).intValue ();
if (this.thisMesh != null) this.thisMesh.atomIndex = this.atomIndex;
} else if ("center" === propertyName) {
this.center.setT (value);
} else if ("colorRGB" === propertyName) {
var rgb = (value).intValue ();
if (rgb == 1297090050) {
this.colorType = rgb;
} else {
this.colorType = 0;
this.defaultColix = JU.C.getColix (rgb);
}} else if ("contour" === propertyName) {
this.explicitContours = true;
} else if ("functionXY" === propertyName) {
if (this.sg.isStateDataRead ()) this.setScriptInfo (null);
} else if ("init" === propertyName) {
this.newSg ();
} else if ("getSurfaceSets" === propertyName) {
if (this.thisMesh != null) {
this.thisMesh.jvxlData.thisSet = (value).intValue ();
this.thisMesh.calculatedVolume = null;
this.thisMesh.calculatedArea = null;
}} else if ("localName" === propertyName) {
value = this.vwr.getOutputChannel (value, null);
propertyName = "outputChannel";
} else if ("molecularOrbital" === propertyName) {
if (Clazz_instanceOf (value, Integer)) {
this.moNumber = (value).intValue ();
this.moLinearCombination = null;
} else {
this.moLinearCombination = value;
this.moNumber = 0;
}if (!this.isColorExplicit) this.isPhaseColored = true;
} else if ("phase" === propertyName) {
this.isPhaseColored = true;
} else if ("plane" === propertyName) {
} else if ("pocket" === propertyName) {
} else if ("scale3d" === propertyName) {
this.scale3d = (value).floatValue ();
if (this.thisMesh != null) {
this.thisMesh.scale3d = this.thisMesh.jvxlData.scale3d = this.scale3d;
this.thisMesh.altVertices = null;
}} else if ("title" === propertyName) {
if (Clazz_instanceOf (value, String) && "-".equals (value)) value = null;
this.setPropertySuper (propertyName, value, bs);
value = this.title;
} else if ("withinPoints" === propertyName) {
var o = value;
this.withinDistance2 = (o[0]).floatValue ();
this.isWithinNot = (this.withinDistance2 < 0);
this.withinDistance2 *= this.withinDistance2;
this.withinPoints = o[3];
if (this.withinPoints.size () == 0) this.withinPoints = this.vwr.ms.getAtomPointVector (o[2]);
} else if (("nci" === propertyName || "orbital" === propertyName) && this.sg != null) {
this.sg.getParams ().testFlags = (this.vwr.getTestFlag (2) ? 2 : 0);
}if (this.sg != null && this.sg.setProp (propertyName, value, bs)) {
if (this.sg.isValid ()) return;
propertyName = "delete";
}if ("init" === propertyName) {
this.explicitID = false;
this.scriptAppendix = "";
var script = (Clazz_instanceOf (value, String) ? value : null);
var pt = (script == null ? -1 : script.indexOf ("# ID="));
this.actualID = (pt >= 0 ? JU.PT.getQuotedStringAt (script, pt) : null);
this.setPropertySuper ("thisID", "+PREVIOUS_MESH+", null);
if (script != null && !(this.iHaveBitSets = this.getScriptBitSets (script, null))) this.sg.setParameter ("select", bs);
this.initializeIsosurface ();
this.sg.setModelIndex (this.isFixed ? -1 : this.modelIndex);
return;
}if ("clear" === propertyName) {
this.discardTempData (true);
return;
}if ("colorDensity" === propertyName) {
if (value != null && this.currentMesh != null) this.currentMesh.volumeRenderPointSize = (value).floatValue ();
return;
}if (propertyName === "deleteModelAtoms") {
var modelIndex = ((value)[2])[0];
var firstAtomDeleted = ((value)[2])[1];
var nAtomsDeleted = ((value)[2])[2];
for (var i = this.meshCount; --i >= 0; ) {
var m = this.meshes[i];
if (m == null) continue;
if (m.connections != null) {
var iAtom = m.connections[0];
if (iAtom >= firstAtomDeleted + nAtomsDeleted) m.connections[0] = iAtom - nAtomsDeleted;
 else if (iAtom >= firstAtomDeleted) m.connections = null;
}m.connections = null;
if (m.modelIndex == modelIndex) {
this.meshCount--;
if (m === this.currentMesh) this.currentMesh = this.thisMesh = null;
this.meshes = this.isomeshes = JU.AU.deleteElements (this.meshes, i, 1);
} else if (m.modelIndex > modelIndex) {
m.modelIndex--;
if (m.atomIndex >= firstAtomDeleted) m.atomIndex -= nAtomsDeleted;
}}
return;
}this.setPropertySuper (propertyName, value, bs);
}, "~S,~O,JU.BS");
Clazz_defineMethod (c$, "setIsoMeshColor", 
 function (m, color) {
m.jvxlData.baseColor = color;
m.isColorSolid = true;
m.pcs = null;
m.colorEncoder = null;
m.vertexColorMap = null;
}, "J.shapesurface.IsosurfaceMesh,~S");
Clazz_defineMethod (c$, "setColorPhase", 
 function (m, colix0, colix1) {
m.colorPhased = true;
m.colix = m.jvxlData.minColorIndex = colix0;
m.jvxlData.maxColorIndex = colix1;
m.jvxlData.isBicolorMap = true;
m.jvxlData.colorDensity = false;
m.isColorSolid = false;
m.remapColors (this.vwr, null, this.translucentLevel);
}, "J.shapesurface.IsosurfaceMesh,~N,~N");
Clazz_defineMethod (c$, "ensureMeshSource", 
 function () {
var haveColors = (this.thisMesh.vertexSource != null);
if (haveColors) for (var i = this.thisMesh.vc; --i >= 0; ) if (this.thisMesh.vertexSource[i] < 0) {
haveColors = false;
break;
}
if (!haveColors) {
var source = this.thisMesh.vertexSource;
var vertexColixes = this.thisMesh.vcs;
var colix = (this.thisMesh.isColorSolid ? this.thisMesh.colix : 0);
this.setProperty ("init", null, null);
this.setProperty ("map", Boolean.FALSE, null);
this.setProperty ("property",  Clazz_newFloatArray (this.vwr.getAtomCount (), 0), null);
if (colix != 0) {
this.thisMesh.colorCommand = "color isosurface " + JU.C.getHexCode (colix);
this.setProperty ("color", Integer.$valueOf (JU.C.getArgb (colix)), null);
}if (source != null) {
for (var i = this.thisMesh.vc; --i >= 0; ) if (source[i] < 0) source[i] = this.thisMesh.vertexSource[i];

this.thisMesh.vertexSource = source;
this.thisMesh.vcs = vertexColixes;
}}});
Clazz_defineMethod (c$, "slabPolygons", 
function (slabInfo) {
this.thisMesh.slabPolygons (slabInfo, false);
this.thisMesh.reinitializeLightingAndColor (this.vwr);
}, "~A");
Clazz_defineMethod (c$, "setPropertySuper", 
 function (propertyName, value, bs) {
if (propertyName === "thisID" && this.currentMesh != null && this.currentMesh.thisID != null && this.currentMesh.thisID.equals (value)) {
this.checkExplicit (value);
return;
}this.currentMesh = this.thisMesh;
this.setPropMC (propertyName, value, bs);
this.thisMesh = this.currentMesh;
this.jvxlData = (this.thisMesh == null ? null : this.thisMesh.jvxlData);
if (this.sg != null) this.sg.setJvxlData (this.jvxlData);
}, "~S,~O,JU.BS");
Clazz_overrideMethod (c$, "getPropertyData", 
function (property, data) {
if (property === "colorEncoder") {
var mesh = this.getMesh (data[0]);
if (mesh == null || (data[1] = mesh.colorEncoder) == null) return false;
return true;
}if (property === "intersectPlane") {
var mesh = this.getMesh (data[0]);
if (mesh == null) return false;
data[3] = Integer.$valueOf (mesh.modelIndex);
mesh.getIntersection (0, data[1], null, data[2], null, null, null, false, false, 135266319, false);
return true;
}if (property === "getBoundingBox") {
var id = data[0];
var m = this.getMesh (id);
if (m == null || m.vs == null) return false;
data[2] = m.jvxlData.boundingBox;
if (m.mat4 != null) {
var d =  new Array (2);
d[0] = JU.P3.newP (m.jvxlData.boundingBox[0]);
d[1] = JU.P3.newP (m.jvxlData.boundingBox[1]);
var v =  new JU.V3 ();
m.mat4.getTranslation (v);
d[0].add (v);
d[1].add (v);
data[2] = d;
}return true;
}if (property === "unitCell") {
var m = this.getMesh (data[0]);
return (m != null && (data[1] = m.getUnitCell ()) != null);
}if (property === "getCenter") {
var index = (data[1]).intValue ();
if (index == -2147483648) {
var id = data[0];
var m = this.getMesh (id);
if (m == null || m.vs == null) return false;
var p = JU.P3.newP (m.jvxlData.boundingBox[0]);
p.add (m.jvxlData.boundingBox[1]);
p.scale (0.5);
if (m.mat4 != null) {
var v =  new JU.V3 ();
m.mat4.getTranslation (v);
p.add (v);
}data[2] = p;
return true;
}}return this.getPropDataMC (property, data);
}, "~S,~A");
Clazz_overrideMethod (c$, "getProperty", 
function (property, index) {
return this.getPropI (property);
}, "~S,~N");
Clazz_defineMethod (c$, "getPropI", 
function (property) {
var ret = this.getPropMC (property);
if (ret != null) return ret;
if (property === "dataRange") return (this.thisMesh == null || this.jvxlData.jvxlPlane != null && this.thisMesh.colorEncoder == null ? null : [this.jvxlData.mappedDataMin, this.jvxlData.mappedDataMax, (this.jvxlData.isColorReversed ? this.jvxlData.valueMappedToBlue : this.jvxlData.valueMappedToRed), (this.jvxlData.isColorReversed ? this.jvxlData.valueMappedToRed : this.jvxlData.valueMappedToBlue)]);
if (property === "moNumber") return Integer.$valueOf (this.moNumber);
if (property === "moLinearCombination") return this.moLinearCombination;
if (property === "nSets") return Integer.$valueOf (this.thisMesh == null ? 0 : this.thisMesh.nSets);
if (property === "area") return (this.thisMesh == null ? Float.$valueOf (NaN) : this.calculateVolumeOrArea (true));
if (property === "volume") return (this.thisMesh == null ? Float.$valueOf (NaN) : this.calculateVolumeOrArea (false));
if (this.thisMesh == null) return null;
if (property === "cutoff") return Float.$valueOf (this.jvxlData.cutoff);
if (property === "minMaxInfo") return [this.jvxlData.dataMin, this.jvxlData.dataMax];
if (property === "plane") return this.jvxlData.jvxlPlane;
if (property === "contours") return this.thisMesh.getContours ();
if (property === "jvxlDataXml" || property === "jvxlMeshXml") {
var meshData = null;
this.jvxlData.slabInfo = null;
if (property === "jvxlMeshXml" || this.jvxlData.vertexDataOnly || this.thisMesh.bsSlabDisplay != null && this.thisMesh.bsSlabGhost == null) {
meshData =  new J.jvxl.data.MeshData ();
this.fillMeshData (meshData, 1, null);
meshData.polygonColorData = J.shapesurface.Isosurface.getPolygonColorData (meshData.pc, meshData.pcs, meshData.bsSlabDisplay);
} else if (this.thisMesh.bsSlabGhost != null) {
this.jvxlData.slabInfo = this.thisMesh.slabOptions.toString ();
}var sb =  new JU.SB ();
this.getMeshCommand (sb, this.thisMesh.index);
this.thisMesh.setJvxlColorMap (true);
return J.jvxl.data.JvxlCoder.jvxlGetFile (this.jvxlData, meshData, this.title, "", true, 1, sb.toString (), null);
}if (property === "jvxlFileInfo") {
this.thisMesh.setJvxlColorMap (false);
return J.jvxl.data.JvxlCoder.jvxlGetInfo (this.jvxlData);
}if (property === "command") {
var sb =  new JU.SB ();
var list = this.getMeshList (this.previousMeshID, false);
for (var i = list.size (); --i >= 0; ) this.getMeshCommand (sb, i);

return sb.toString ();
}return null;
}, "~S");
Clazz_defineMethod (c$, "calculateVolumeOrArea", 
 function (isArea) {
if (isArea) {
if (this.thisMesh.calculatedArea != null) return this.thisMesh.calculatedArea;
} else {
if (this.thisMesh.calculatedVolume != null) return this.thisMesh.calculatedVolume;
}var meshData =  new J.jvxl.data.MeshData ();
this.fillMeshData (meshData, 1, null);
meshData.nSets = this.thisMesh.nSets;
meshData.vertexSets = this.thisMesh.vertexSets;
if (!isArea && this.thisMesh.jvxlData.colorDensity) {
var f = this.thisMesh.jvxlData.voxelVolume;
f *= (this.thisMesh.bsSlabDisplay == null ? this.thisMesh.vc : this.thisMesh.bsSlabDisplay.cardinality ());
return this.thisMesh.calculatedVolume = Float.$valueOf (f);
}var ret = meshData.calculateVolumeOrArea (this.thisMesh.jvxlData.thisSet, isArea, false);
if (isArea) this.thisMesh.calculatedArea = ret;
 else this.thisMesh.calculatedVolume = ret;
return ret;
}, "~B");
c$.getPolygonColorData = Clazz_defineMethod (c$, "getPolygonColorData", 
function (ccount, colixes, bsSlabDisplay) {
if (colixes == null) return null;
var list1 =  new JU.SB ();
var count = 0;
var colix = 0;
var done = false;
for (var i = 0; i < ccount || (done = true) == true; i++) {
if (!done && bsSlabDisplay != null && !bsSlabDisplay.get (i)) continue;
if (done || colixes[i] != colix) {
if (count != 0) list1.append (" ").appendI (count).append (" ").appendI ((colix == 0 ? 0 : JU.C.getArgb (colix)));
if (done) break;
colix = colixes[i];
count = 1;
} else {
count++;
}}
list1.append ("\n");
return list1.toString ();
}, "~N,~A,JU.BS");
Clazz_overrideMethod (c$, "getShapeState", 
function () {
this.clean ();
var sb =  new JU.SB ();
sb.append ("\n");
for (var i = 0; i < this.meshCount; i++) this.getMeshCommand (sb, i);

return sb.toString ();
});
Clazz_defineMethod (c$, "getMeshCommand", 
 function (sb, i) {
var imesh = this.meshes[i];
if (imesh == null || imesh.scriptCommand == null) return;
var cmd = imesh.scriptCommand;
var modelCount = this.vwr.getModelCount ();
if (modelCount > 1) J.shape.Shape.appendCmd (sb, "frame " + this.vwr.getModelNumberDotted (imesh.modelIndex));
cmd = JU.PT.rep (cmd, ";; isosurface map", " map");
cmd = JU.PT.rep (cmd, "; isosurface map", " map");
cmd = cmd.$replace ('\t', ' ');
cmd = JU.PT.rep (cmd, ";#", "; #");
var pt = cmd.indexOf ("; #");
if (pt >= 0) cmd = cmd.substring (0, pt);
if (imesh.connections != null) cmd += " connect " + JU.Escape.eAI (imesh.connections);
cmd = JU.PT.trim (cmd, ";");
if (imesh.linkedMesh != null) cmd += " LINK";
if (this.myType === "lcaoCartoon" && imesh.atomIndex >= 0) cmd += " ATOMINDEX " + imesh.atomIndex;
J.shape.Shape.appendCmd (sb, cmd);
var id = this.myType + " ID " + JU.PT.esc (imesh.thisID);
if (imesh.jvxlData.thisSet >= 0) J.shape.Shape.appendCmd (sb, id + " set " + (imesh.jvxlData.thisSet + 1));
if (imesh.mat4 != null) J.shape.Shape.appendCmd (sb, id + " move " + JU.Escape.matrixToScript (imesh.mat4));
if (imesh.scale3d != 0) J.shape.Shape.appendCmd (sb, id + " scale3d " + imesh.scale3d);
if (imesh.jvxlData.slabValue != -2147483648) J.shape.Shape.appendCmd (sb, id + " slab " + imesh.jvxlData.slabValue);
if (imesh.slabOptions != null) J.shape.Shape.appendCmd (sb, imesh.slabOptions.toString ());
if (cmd.charAt (0) != '#') {
if (this.allowMesh) J.shape.Shape.appendCmd (sb, imesh.getState (this.myType));
if (!imesh.isColorSolid && imesh.colorType == 0 && JU.C.isColixTranslucent (imesh.colix)) J.shape.Shape.appendCmd (sb, "color " + this.myType + " " + J.shape.Shape.getTranslucentLabel (imesh.colix));
if (imesh.colorCommand != null && imesh.colorType == 0 && !imesh.colorCommand.equals ("#inherit;")) {
J.shape.Shape.appendCmd (sb, imesh.colorCommand);
}var colorArrayed = (imesh.isColorSolid && imesh.pcs != null);
if (imesh.isColorSolid && imesh.colorType == 0 && !colorArrayed) {
J.shape.Shape.appendCmd (sb, J.shape.Shape.getColorCommandUnk (this.myType, imesh.colix, this.translucentAllowed));
} else if (imesh.jvxlData.isBicolorMap && imesh.colorPhased) {
J.shape.Shape.appendCmd (sb, "color isosurface phase " + J.shape.Shape.encodeColor (imesh.jvxlData.minColorIndex) + " " + J.shape.Shape.encodeColor (imesh.jvxlData.maxColorIndex));
}if (imesh.vertexColorMap != null) for (var entry, $entry = imesh.vertexColorMap.entrySet ().iterator (); $entry.hasNext () && ((entry = $entry.next ()) || true);) {
var bs = entry.getValue ();
if (!bs.isEmpty ()) J.shape.Shape.appendCmd (sb, "color " + this.myType + " " + JU.Escape.eBS (bs) + " " + entry.getKey ());
}
}}, "JU.SB,~N");
Clazz_defineMethod (c$, "getScriptBitSets", 
 function (script, bsCmd) {
this.script = script;
var i;
this.iHaveModelIndex = false;
this.modelIndex = -1;
if (script != null && (i = script.indexOf ("MODEL({")) >= 0) {
var j = script.indexOf ("})", i);
if (j > 0) {
var bs = JU.BS.unescape (script.substring (i + 3, j + 1));
this.modelIndex = (bs == null ? -1 : bs.nextSetBit (0));
this.iHaveModelIndex = (this.modelIndex >= 0);
}}if (script == null) return false;
this.getCapSlabInfo (script);
i = script.indexOf ("# ({");
if (i < 0) return false;
var j = script.indexOf ("})", i);
if (j < 0) return false;
var bs = JU.BS.unescape (script.substring (i + 2, j + 2));
if (bsCmd == null) this.sg.setParameter ("select", bs);
 else bsCmd[0] = bs;
if ((i = script.indexOf ("({", j)) < 0) return true;
j = script.indexOf ("})", i);
if (j < 0) return false;
bs = JU.BS.unescape (script.substring (i + 1, j + 1));
if (bsCmd == null) this.sg.setParameter ("ignore", bs);
 else bsCmd[1] = bs;
if ((i = script.indexOf ("/({", j)) == j + 2) {
if ((j = script.indexOf ("})", i)) < 0) return false;
bs = JU.BS.unescape (script.substring (i + 3, j + 1));
if (bsCmd == null) this.vwr.ms.setTrajectoryBs (bs);
 else bsCmd[2] = bs;
}return true;
}, "~S,~A");
Clazz_defineMethod (c$, "getCapSlabInfo", 
function (script) {
var i = script.indexOf ("# SLAB=");
if (i >= 0) this.sg.setParameter ("slab", JU.MeshSurface.getCapSlabObject (JU.PT.getQuotedStringAt (script, i), false));
i = script.indexOf ("# CAP=");
if (i >= 0) this.sg.setParameter ("slab", JU.MeshSurface.getCapSlabObject (JU.PT.getQuotedStringAt (script, i), true));
}, "~S");
Clazz_defineMethod (c$, "initializeIsosurface", 
 function () {
if (!this.iHaveModelIndex) this.modelIndex = this.vwr.am.cmi;
this.atomIndex = -1;
this.bsDisplay = null;
this.center = JU.P3.new3 (3.4028235E38, 3.4028235E38, 3.4028235E38);
this.colix = 5;
this.connections = null;
this.cutoffRange = null;
this.colorType = this.defaultColix = this.meshColix = 0;
this.displayWithinPoints = null;
this.explicitContours = false;
this.isFixed = (this.modelIndex < 0);
this.isPhaseColored = this.isColorExplicit = false;
this.linkedMesh = null;
if (this.modelIndex < 0) this.modelIndex = 0;
this.scale3d = 0;
this.title = null;
this.translucentLevel = 0;
this.withinPoints = null;
this.initState ();
});
Clazz_defineMethod (c$, "initState", 
 function () {
this.associateNormals = true;
this.sg.initState ();
});
Clazz_defineMethod (c$, "setMeshI", 
 function () {
this.thisMesh.visible = true;
if ((this.thisMesh.atomIndex = this.atomIndex) >= 0) this.thisMesh.modelIndex = this.vwr.getAtomModelIndex (this.atomIndex);
 else if (this.isFixed) this.thisMesh.modelIndex = -1;
 else if (this.modelIndex >= 0) this.thisMesh.modelIndex = this.modelIndex;
 else this.thisMesh.modelIndex = this.vwr.am.cmi;
this.thisMesh.scriptCommand = this.script;
this.thisMesh.ptCenter.setT (this.center);
this.thisMesh.scale3d = (this.thisMesh.jvxlData.jvxlPlane == null ? 0 : this.scale3d);
});
Clazz_defineMethod (c$, "discardTempData", 
function (discardAll) {
if (!discardAll) return;
this.title = null;
if (this.thisMesh == null) return;
this.thisMesh.surfaceSet = null;
}, "~B");
Clazz_defineMethod (c$, "getDefaultColix", 
 function () {
if (this.defaultColix != 0) return this.defaultColix;
if (!this.sg.isCubeData ()) return this.colix;
var argb = (this.sg.getCutoff () >= 0 ? -11525984 : -6283184);
return JU.C.getColix (argb);
});
Clazz_defineMethod (c$, "drawLcaoCartoon", 
 function (z, x, rotAxis, nElectrons) {
var lcaoCartoon = this.sg.setLcao ();
var rotRadians = rotAxis.x + rotAxis.y + rotAxis.z;
this.defaultColix = JU.C.getColix (this.sg.getColor (1));
var colixNeg = JU.C.getColix (this.sg.getColor (-1));
var y =  new JU.V3 ();
var isReverse = (lcaoCartoon.length > 0 && lcaoCartoon.charAt (0) == '-');
if (isReverse) lcaoCartoon = lcaoCartoon.substring (1);
var sense = (isReverse ? -1 : 1);
y.cross (z, x);
if (rotRadians != 0) {
var a =  new JU.A4 ();
if (rotAxis.x != 0) a.setVA (x, rotRadians);
 else if (rotAxis.y != 0) a.setVA (y, rotRadians);
 else a.setVA (z, rotRadians);
var m =  new JU.M3 ().setAA (a);
m.rotate (x);
m.rotate (y);
m.rotate (z);
}if (this.thisMesh == null && this.nLCAO == 0) this.nLCAO = this.meshCount;
var id = (this.thisMesh == null ? (nElectrons > 0 ? "lp" : "lcao") + (++this.nLCAO) + "_" + lcaoCartoon : this.thisMesh.thisID);
if (this.thisMesh == null) this.allocMesh (id, null);
if (lcaoCartoon.equals ("px")) {
this.thisMesh.thisID += "a";
var meshA = this.thisMesh;
this.createLcaoLobe (x, sense, nElectrons);
if (nElectrons > 0) return;
this.setProperty ("thisID", id + "b", null);
this.createLcaoLobe (x, -sense, nElectrons);
this.thisMesh.colix = colixNeg;
this.linkedMesh = this.thisMesh.linkedMesh = meshA;
return;
}if (lcaoCartoon.equals ("py")) {
this.thisMesh.thisID += "a";
var meshA = this.thisMesh;
this.createLcaoLobe (y, sense, nElectrons);
if (nElectrons > 0) return;
this.setProperty ("thisID", id + "b", null);
this.createLcaoLobe (y, -sense, nElectrons);
this.thisMesh.colix = colixNeg;
this.linkedMesh = this.thisMesh.linkedMesh = meshA;
return;
}if (lcaoCartoon.equals ("pz")) {
this.thisMesh.thisID += "a";
var meshA = this.thisMesh;
this.createLcaoLobe (z, sense, nElectrons);
if (nElectrons > 0) return;
this.setProperty ("thisID", id + "b", null);
this.createLcaoLobe (z, -sense, nElectrons);
this.thisMesh.colix = colixNeg;
this.linkedMesh = this.thisMesh.linkedMesh = meshA;
return;
}if (lcaoCartoon.equals ("pza") || lcaoCartoon.indexOf ("sp") == 0 || lcaoCartoon.indexOf ("d") == 0 || lcaoCartoon.indexOf ("lp") == 0) {
this.createLcaoLobe (z, sense, nElectrons);
return;
}if (lcaoCartoon.equals ("pzb")) {
this.createLcaoLobe (z, -sense, nElectrons);
return;
}if (lcaoCartoon.equals ("pxa")) {
this.createLcaoLobe (x, sense, nElectrons);
return;
}if (lcaoCartoon.equals ("pxb")) {
this.createLcaoLobe (x, -sense, nElectrons);
return;
}if (lcaoCartoon.equals ("pya")) {
this.createLcaoLobe (y, sense, nElectrons);
return;
}if (lcaoCartoon.equals ("pyb")) {
this.createLcaoLobe (y, -sense, nElectrons);
return;
}if (lcaoCartoon.equals ("spacefill") || lcaoCartoon.equals ("cpk")) {
this.createLcaoLobe (null, 2 * this.vwr.getAtomRadius (this.atomIndex), nElectrons);
return;
}this.createLcaoLobe (null, 1, nElectrons);
return;
}, "JU.V3,JU.V3,JU.V3,~N");
Clazz_defineMethod (c$, "createLcaoLobe", 
 function (lobeAxis, factor, nElectrons) {
this.initState ();
if (JU.Logger.debugging) {
JU.Logger.debug ("creating isosurface ID " + this.thisMesh.thisID);
}if (lobeAxis == null) {
this.setProperty ("sphere", Float.$valueOf (factor / 2), null);
} else {
this.lcaoDir.x = lobeAxis.x * factor;
this.lcaoDir.y = lobeAxis.y * factor;
this.lcaoDir.z = lobeAxis.z * factor;
this.lcaoDir.w = 0.7;
this.setProperty (nElectrons == 2 ? "lp" : nElectrons == 1 ? "rad" : "lobe", this.lcaoDir, null);
}this.thisMesh.colix = this.defaultColix;
this.setScriptInfo (null);
}, "JU.V3,~N,~N");
Clazz_overrideMethod (c$, "invalidateTriangles", 
function () {
this.thisMesh.invalidatePolygons ();
});
Clazz_overrideMethod (c$, "setOutputChannel", 
function (binaryDoc, out) {
binaryDoc.setOutputChannel (out);
}, "javajs.api.GenericBinaryDocument,JU.OC");
Clazz_overrideMethod (c$, "fillMeshData", 
function (meshData, mode, mesh) {
if (meshData == null) {
if (this.thisMesh == null) this.allocMesh (null, null);
if (!this.thisMesh.isMerged) this.thisMesh.clearType (this.myType, this.sg.getIAddGridPoints ());
this.thisMesh.connections = this.connections;
this.thisMesh.colix = this.getDefaultColix ();
this.thisMesh.colorType = this.colorType;
this.thisMesh.meshColix = this.meshColix;
if (this.isPhaseColored || this.thisMesh.jvxlData.isBicolorMap) this.thisMesh.isColorSolid = false;
return;
}if (mesh == null) mesh = this.thisMesh;
if (mesh == null) return;
switch (mode) {
case 1:
meshData.mergeVertexCount0 = mesh.mergeVertexCount0;
meshData.vs = mesh.vs;
meshData.vertexSource = mesh.vertexSource;
meshData.vvs = mesh.vvs;
meshData.vc = mesh.vc;
meshData.vertexIncrement = mesh.vertexIncrement;
meshData.pc = mesh.pc;
meshData.pis = mesh.pis;
meshData.pcs = mesh.pcs;
meshData.bsSlabDisplay = mesh.bsSlabDisplay;
meshData.bsSlabGhost = mesh.bsSlabGhost;
meshData.slabColix = mesh.slabColix;
meshData.slabMeshType = mesh.slabMeshType;
meshData.polygonCount0 = mesh.polygonCount0;
meshData.vertexCount0 = mesh.vertexCount0;
meshData.slabOptions = mesh.slabOptions;
return;
case 2:
if (mesh.vcs == null || mesh.vc > mesh.vcs.length) mesh.vcs =  Clazz_newShortArray (mesh.vc, 0);
meshData.vcs = mesh.vcs;
return;
case 3:
mesh.surfaceSet = meshData.surfaceSet;
mesh.vertexSets = meshData.vertexSets;
mesh.nSets = meshData.nSets;
return;
case 4:
mesh.vs = meshData.vs;
mesh.vvs = meshData.vvs;
mesh.vc = meshData.vc;
mesh.vertexIncrement = meshData.vertexIncrement;
mesh.vertexSource = meshData.vertexSource;
mesh.pc = meshData.pc;
mesh.pis = meshData.pis;
mesh.pcs = meshData.pcs;
mesh.bsSlabDisplay = meshData.bsSlabDisplay;
mesh.bsSlabGhost = meshData.bsSlabGhost;
mesh.slabColix = meshData.slabColix;
mesh.slabMeshType = meshData.slabMeshType;
mesh.polygonCount0 = meshData.polygonCount0;
mesh.vertexCount0 = meshData.vertexCount0;
mesh.mergeVertexCount0 = meshData.mergeVertexCount0;
mesh.slabOptions = meshData.slabOptions;
return;
}
}, "J.jvxl.data.MeshData,~N,J.shapesurface.IsosurfaceMesh");
Clazz_overrideMethod (c$, "notifySurfaceGenerationCompleted", 
function () {
this.setMeshI ();
this.setBsVdw ();
this.thisMesh.insideOut = this.sg.isInsideOut ();
this.thisMesh.vertexSource = this.sg.getVertexSource ();
this.thisMesh.spanningVectors = this.sg.getSpanningVectors ();
this.thisMesh.calculatedArea = null;
this.thisMesh.calculatedVolume = null;
var params = this.sg.getParams ();
if (!this.thisMesh.isMerged) this.thisMesh.initialize (this.sg.isFullyLit () ? 1073741964 : 1073741958, null, this.sg.getPlane ());
if (!params.allowVolumeRender) this.thisMesh.jvxlData.allowVolumeRender = false;
this.thisMesh.setColorsFromJvxlData (this.sg.getParams ().colorRgb);
if (this.thisMesh.jvxlData.slabInfo != null) this.vwr.runScript ("isosurface " + this.thisMesh.jvxlData.slabInfo);
if (this.sg.getParams ().psi_monteCarloCount > 0) this.thisMesh.diameter = -1;
});
Clazz_overrideMethod (c$, "notifySurfaceMappingCompleted", 
function () {
if (!this.thisMesh.isMerged) {
this.thisMesh.initialize (this.sg.isFullyLit () ? 1073741964 : 1073741958, null, this.sg.getPlane ());
this.thisMesh.setJvxlDataRendering ();
}this.setBsVdw ();
this.thisMesh.isColorSolid = false;
this.thisMesh.colorDensity = this.jvxlData.colorDensity;
this.thisMesh.volumeRenderPointSize = this.jvxlData.pointSize;
this.thisMesh.colorEncoder = this.sg.getColorEncoder ();
this.thisMesh.getContours ();
if (this.thisMesh.jvxlData.nContours != 0 && this.thisMesh.jvxlData.nContours != -1) this.explicitContours = true;
if (this.explicitContours && this.thisMesh.jvxlData.jvxlPlane != null) this.thisMesh.havePlanarContours = true;
this.setPropertySuper ("token", Integer.$valueOf (this.explicitContours ? 1073742046 : 1073741938), null);
this.setPropertySuper ("token", Integer.$valueOf (this.explicitContours ? 1073741898 : 1073742039), null);
var slabInfo = this.sg.getSlabInfo ();
if (slabInfo != null) {
this.thisMesh.slabPolygonsList (slabInfo, false);
this.thisMesh.reinitializeLightingAndColor (this.vwr);
}this.thisMesh.setColorCommand ();
});
Clazz_defineMethod (c$, "setBsVdw", 
 function () {
var bs = this.sg.geVdwBitSet ();
if (bs == null) return;
if (this.thisMesh.bsVdw == null) this.thisMesh.bsVdw =  new JU.BS ();
this.thisMesh.bsVdw.or (bs);
});
Clazz_overrideMethod (c$, "calculateGeodesicSurface", 
function (bsSelected, envelopeRadius) {
return this.vwr.calculateSurface (bsSelected, envelopeRadius);
}, "JU.BS,~N");
Clazz_overrideMethod (c$, "getSurfacePointIndexAndFraction", 
function (cutoff, isCutoffAbsolute, x, y, z, offset, vA, vB, valueA, valueB, pointA, edgeVector, isContourType, fReturn) {
return 0;
}, "~N,~B,~N,~N,~N,JU.P3i,~N,~N,~N,~N,JU.T3,JU.V3,~B,~A");
Clazz_overrideMethod (c$, "addVertexCopy", 
function (vertexXYZ, value, assocVertex, asCopy) {
if (this.cutoffRange != null && (value < this.cutoffRange[0] || value > this.cutoffRange[1])) return -1;
return (this.withinPoints != null && !J.shape.Mesh.checkWithin (vertexXYZ, this.withinPoints, this.withinDistance2, this.isWithinNot) ? -1 : this.thisMesh.addVertexCopy (vertexXYZ, value, assocVertex, this.associateNormals, asCopy));
}, "JU.T3,~N,~N,~B");
Clazz_overrideMethod (c$, "addTriangleCheck", 
function (iA, iB, iC, check, check2, isAbsolute, color) {
return (iA < 0 || iB < 0 || iC < 0 || isAbsolute && !J.jvxl.data.MeshData.checkCutoff (iA, iB, iC, this.thisMesh.vvs) ? -1 : this.thisMesh.addTriangleCheck (iA, iB, iC, check, check2, color));
}, "~N,~N,~N,~N,~N,~B,~N");
Clazz_defineMethod (c$, "setScriptInfo", 
function (strCommand) {
var script = (strCommand == null ? this.sg.getScript () : strCommand);
var pt = (script == null ? -1 : script.indexOf ("; isosurface map"));
if (pt == 0) {
if (this.thisMesh.scriptCommand == null) return;
pt = this.thisMesh.scriptCommand.indexOf ("; isosurface map");
if (pt >= 0) this.thisMesh.scriptCommand = this.thisMesh.scriptCommand.substring (0, pt);
this.thisMesh.scriptCommand += script;
return;
}this.thisMesh.title = this.sg.getTitle ();
this.thisMesh.dataType = this.sg.getParams ().dataType;
this.thisMesh.scale3d = this.sg.getParams ().scale3d;
if (script != null) {
if (script.charAt (0) == ' ') {
script = this.myType + " ID " + JU.PT.esc (this.thisMesh.thisID) + script;
pt = script.indexOf ("; isosurface map");
}}if (pt > 0 && this.scriptAppendix.length > 0) this.thisMesh.scriptCommand = script.substring (0, pt) + this.scriptAppendix + script.substring (pt);
 else this.thisMesh.scriptCommand = script + this.scriptAppendix;
if (!this.explicitID && script != null && (pt = script.indexOf ("# ID=")) >= 0) this.thisMesh.thisID = JU.PT.getQuotedStringAt (script, pt);
}, "~S");
Clazz_overrideMethod (c$, "addRequiredFile", 
function (fileName) {
fileName = " # /*file*/\"" + fileName + "\"";
if (this.scriptAppendix.indexOf (fileName) < 0) this.scriptAppendix += fileName;
}, "~S");
Clazz_defineMethod (c$, "setJvxlInfo", 
 function () {
if (this.sg.getJvxlData () !== this.jvxlData || this.sg.getJvxlData () !== this.thisMesh.jvxlData) this.jvxlData = this.thisMesh.jvxlData = this.sg.getJvxlData ();
});
Clazz_overrideMethod (c$, "getShapeDetail", 
function () {
var V =  new JU.Lst ();
for (var i = 0; i < this.meshCount; i++) {
var info =  new java.util.Hashtable ();
var mesh = this.isomeshes[i];
if (mesh == null || mesh.vs == null || mesh.vc == 0 && mesh.pc == 0) continue;
this.addMeshInfo (mesh, info);
V.addLast (info);
}
return V;
});
Clazz_defineMethod (c$, "addMeshInfo", 
function (mesh, info) {
info.put ("ID", (mesh.thisID == null ? "<noid>" : mesh.thisID));
info.put ("visible", Boolean.$valueOf (mesh.visible));
info.put ("vertexCount", Integer.$valueOf (mesh.vc));
if (mesh.calculatedVolume != null) info.put ("volume", mesh.calculatedVolume);
if (mesh.calculatedArea != null) info.put ("area", mesh.calculatedArea);
if (mesh.ptCenter.x != 3.4028235E38) info.put ("center", mesh.ptCenter);
if (mesh.mat4 != null) info.put ("mat4", mesh.mat4);
if (mesh.scale3d != 0) info.put ("scale3d", Float.$valueOf (mesh.scale3d));
info.put ("xyzMin", mesh.jvxlData.boundingBox[0]);
info.put ("xyzMax", mesh.jvxlData.boundingBox[1]);
var s = J.jvxl.data.JvxlCoder.jvxlGetInfo (mesh.jvxlData);
if (s != null) info.put ("jvxlInfo", s.$replace ('\n', ' '));
info.put ("modelIndex", Integer.$valueOf (mesh.modelIndex));
info.put ("color", JU.CU.colorPtFromInt (JU.C.getArgb (mesh.colix), null));
if (mesh.colorEncoder != null) info.put ("colorKey", mesh.colorEncoder.getColorKey ());
if (mesh.title != null) info.put ("title", mesh.title);
if (mesh.jvxlData.contourValues != null || mesh.jvxlData.contourValuesUsed != null) info.put ("contours", mesh.getContourList (this.vwr));
}, "J.shapesurface.IsosurfaceMesh,java.util.Map");
Clazz_overrideMethod (c$, "getPlane", 
function (x) {
return null;
}, "~N");
Clazz_overrideMethod (c$, "getValue", 
function (x, y, z, ptyz) {
return 0;
}, "~N,~N,~N,~N");
Clazz_overrideMethod (c$, "checkObjectHovered", 
function (x, y, bsVisible) {
if (this.keyXy != null && x >= this.keyXy[0] && y >= this.keyXy[1] && x < this.keyXy[2] && y < this.keyXy[3]) {
this.hoverKey (x, y);
return true;
}if (!this.vwr.getDrawHover ()) return false;
var s = this.findValue (x, y, false, bsVisible);
if (s == null) return false;
if (this.gdata.isDisplayAntialiased ()) {
x <<= 1;
y <<= 1;
}this.vwr.hoverOnPt (x, y, s, this.pickedMesh.thisID, this.pickedPt);
return true;
}, "~N,~N,JU.BS");
Clazz_defineMethod (c$, "hoverKey", 
 function (x, y) {
try {
var s;
var f = 1 - 1.0 * (y - this.keyXy[1]) / (this.keyXy[3] - this.keyXy[1]);
if (this.thisMesh.showContourLines) {
var vContours = this.thisMesh.getContours ();
if (vContours == null) {
if (this.thisMesh.jvxlData.contourValues == null) return;
var i = Clazz_doubleToInt (Math.floor (f * this.thisMesh.jvxlData.contourValues.length));
if (i < 0 || i > this.thisMesh.jvxlData.contourValues.length) return;
s = "" + this.thisMesh.jvxlData.contourValues[i];
} else {
var i = Clazz_doubleToInt (Math.floor (f * vContours.length));
if (i < 0 || i > vContours.length) return;
s = "" + (vContours[i].get (2)).floatValue ();
}} else {
var g = this.thisMesh.colorEncoder.quantize (f, true);
f = this.thisMesh.colorEncoder.quantize (f, false);
s = "" + g + " - " + f;
}if (this.gdata.isAntialiased ()) {
x <<= 1;
y <<= 1;
}this.vwr.hoverOnPt (x, y, s, null, null);
} catch (e) {
if (Clazz_exceptionOf (e, Exception)) {
} else {
throw e;
}
}
}, "~N,~N");
Clazz_overrideMethod (c$, "checkObjectClicked", 
function (x, y, action, bsVisible, drawPicking) {
if (!drawPicking) return null;
if (!this.vwr.isBound (action, 18)) return null;
var dmin2 = 100;
if (this.gdata.isAntialiased ()) {
x <<= 1;
y <<= 1;
dmin2 <<= 1;
}var imesh = -1;
var jmaxz = -1;
var jminz = -1;
var maxz = -2147483648;
var minz = 2147483647;
var pickFront = true;
for (var i = 0; i < this.meshCount; i++) {
var m = this.isomeshes[i];
if (!this.isPickable (m, bsVisible)) continue;
var centers = (pickFront ? m.vs : m.getCenters ());
if (centers == null) continue;
for (var j = centers.length; --j >= 0; ) {
var v = centers[j];
if (v == null) continue;
var d2 = this.coordinateInRange (x, y, v, dmin2, this.ptXY);
if (d2 >= 0) {
if (this.ptXY.z < minz) {
if (pickFront) imesh = i;
minz = this.ptXY.z;
jminz = j;
}if (this.ptXY.z > maxz) {
if (!pickFront) imesh = i;
maxz = this.ptXY.z;
jmaxz = j;
}}}
}
if (imesh < 0) return null;
this.pickedMesh = this.isomeshes[imesh];
this.setPropertySuper ("thisID", this.pickedMesh.thisID, null);
var iFace = this.pickedVertex = (pickFront ? jminz : jmaxz);
var ptRet =  new JU.P3 ();
ptRet.setT ((pickFront ? this.pickedMesh.vs[this.pickedVertex] : (this.pickedMesh).centers[iFace]));
this.pickedModel = this.pickedMesh.modelIndex;
var map = this.getPickedPoint (ptRet, this.pickedModel);
this.setStatusPicked (-4, ptRet, map);
return map;
}, "~N,~N,~N,JU.BS,~B");
Clazz_defineMethod (c$, "isPickable", 
 function (m, bsVisible) {
return m.visibilityFlags != 0 && (m.modelIndex < 0 || bsVisible.get (m.modelIndex)) && !JU.C.isColixTranslucent (m.colix);
}, "J.shapesurface.IsosurfaceMesh,JU.BS");
Clazz_defineMethod (c$, "findValue", 
 function (x, y, isPicking, bsVisible) {
var dmin2 = 100;
if (this.gdata.isAntialiased ()) {
x <<= 1;
y <<= 1;
dmin2 <<= 1;
}var pickedVertex = -1;
var pickedContour = null;
var m = null;
for (var i = 0; i < this.meshCount; i++) {
m = this.isomeshes[i];
if (!this.isPickable (m, bsVisible)) continue;
var vs = m.jvxlData.vContours;
var ilast = (m.firstRealVertex < 0 ? 0 : m.firstRealVertex);
var pickedJ = 0;
if (vs != null && vs.length > 0) {
for (var j = 0; j < vs.length; j++) {
var vc = vs[j];
var n = vc.size () - 1;
for (var k = 6; k < n; k++) {
var v = vc.get (k);
var d2 = this.coordinateInRange (x, y, v, dmin2, this.ptXY);
if (d2 >= 0) {
dmin2 = d2;
pickedContour = vc;
pickedJ = j;
this.pickedMesh = m;
this.pickedPt = v;
}}
}
if (pickedContour != null) return pickedContour.get (2).toString () + (JU.Logger.debugging ? " " + pickedJ : "");
} else if (m.jvxlData.jvxlPlane != null && m.vvs != null) {
var vertices = (m.mat4 == null && m.scale3d == 0 ? m.vs : m.getOffsetVertices (m.jvxlData.jvxlPlane));
for (var k = m.vc; --k >= ilast; ) {
var v = vertices[k];
var d2 = this.coordinateInRange (x, y, v, dmin2, this.ptXY);
if (d2 >= 0) {
dmin2 = d2;
pickedVertex = k;
this.pickedMesh = m;
this.pickedPt = v;
}}
if (pickedVertex != -1) break;
} else if (m.vvs != null) {
for (var k = m.vc; --k >= ilast; ) {
var v = m.vs[k];
var d2 = this.coordinateInRange (x, y, v, dmin2, this.ptXY);
if (d2 >= 0) {
dmin2 = d2;
pickedVertex = k;
this.pickedMesh = m;
this.pickedPt = v;
}}
if (pickedVertex != -1) break;
}}
return (pickedVertex == -1 ? null : (JU.Logger.debugging ? "$" + m.thisID + "[" + (pickedVertex + 1) + "] " + m.vs[pickedVertex] + ": " : m.thisID + ": ") + m.vvs[pickedVertex]);
}, "~N,~N,~B,JU.BS");
Clazz_defineMethod (c$, "getCmd", 
function (index) {
var sb =  new JU.SB ().append ("\n");
this.getMeshCommand (sb, index);
return (sb.toString ());
}, "~N");
Clazz_defineStatics (c$,
"MAX_OBJECT_CLICK_DISTANCE_SQUARED", 100);
});
Clazz_declarePackage ("J.jvxl.data");
Clazz_load (null, "J.jvxl.data.JvxlCoder", ["java.lang.Character", "$.Float", "JU.BS", "$.Lst", "$.P3", "$.PT", "$.SB", "$.XmlUtil", "J.jvxl.data.VolumeData", "JU.BSUtil", "$.C", "$.Escape", "$.Logger"], function () {
c$ = Clazz_declareType (J.jvxl.data, "JvxlCoder");
c$.jvxlGetFile = Clazz_defineMethod (c$, "jvxlGetFile", 
function (volumeData, jvxlData, title) {
var counts = volumeData.getVoxelCounts ();
jvxlData.nPointsX = counts[0];
jvxlData.nPointsY = counts[1];
jvxlData.nPointsZ = counts[2];
jvxlData.jvxlVolumeDataXml = volumeData.setVolumetricXml ();
return J.jvxl.data.JvxlCoder.jvxlGetFile (jvxlData, null, title, null, true, 1, null, null);
}, "J.jvxl.data.VolumeData,J.jvxl.data.JvxlData,~A");
c$.jvxlGetFile = Clazz_defineMethod (c$, "jvxlGetFile", 
function (jvxlData, meshData, title, msg, includeHeader, nSurfaces, state, comment) {
return J.jvxl.data.JvxlCoder.jvxlGetFileXml (jvxlData, meshData, title, msg, includeHeader, nSurfaces, state, comment);
}, "J.jvxl.data.JvxlData,J.jvxl.data.MeshData,~A,~S,~B,~N,~S,~S");
c$.jvxlGetFileXml = Clazz_defineMethod (c$, "jvxlGetFileXml", 
 function (jvxlData, meshData, title, msg, includeHeader, nSurfaces, state, comment) {
var data =  new JU.SB ();
if ("TRAILERONLY".equals (msg)) {
JU.XmlUtil.closeTag (data, "jvxlSurfaceSet");
JU.XmlUtil.closeTag (data, "jvxl");
return data.toString ();
}var vertexDataOnly = (meshData != null);
var isHeaderOnly = ("HEADERONLY".equals (msg));
if (includeHeader) {
JU.XmlUtil.openDocument (data);
JU.XmlUtil.openTagAttr (data, "jvxl", ["version", "2.3", "jmolVersion", jvxlData.version, "xmlns", "http://jmol.org/jvxl_schema", "xmlns:cml", "http://www.xml-cml.org/schema"]);
if (jvxlData.jvxlFileTitle != null) JU.XmlUtil.appendCdata (data, "jvxlFileTitle", null, "\n" + jvxlData.jvxlFileTitle);
if (jvxlData.moleculeXml != null) data.append (jvxlData.moleculeXml);
var volumeDataXml = (vertexDataOnly ? null : jvxlData.jvxlVolumeDataXml);
if (volumeDataXml == null) volumeDataXml = ( new J.jvxl.data.VolumeData ()).setVolumetricXml ();
data.append (volumeDataXml);
JU.XmlUtil.openTagAttr (data, "jvxlSurfaceSet", ["count", "" + (nSurfaces > 0 ? nSurfaces : 1)]);
if (isHeaderOnly) return data.toString ();
}var sb;
var type = (vertexDataOnly ? "pmesh" : jvxlData.jvxlPlane == null ? "isosurface" : "plane");
if (jvxlData.jvxlColorData != null && jvxlData.jvxlColorData.length > 0) type = "mapped " + type;
JU.XmlUtil.openTagAttr (data, "jvxlSurface", ["type", type]);
data.append (J.jvxl.data.JvxlCoder.jvxlGetInfoData (jvxlData, vertexDataOnly));
J.jvxl.data.JvxlCoder.jvxlAppendCommandState (data, comment, state);
if (title != null || msg != null && msg.length > 0) {
sb =  new JU.SB ();
if (msg != null && msg.length > 0) sb.append (msg).append ("\n");
if (title != null) for (var i = 0; i < title.length; i++) sb.append (title[i]).appendC ('\n');

JU.XmlUtil.appendCdata (data, "jvxlSurfaceTitle", null, sb.toString ());
}sb =  new JU.SB ();
JU.XmlUtil.openTagAttr (sb, "jvxlSurfaceData", (vertexDataOnly || jvxlData.jvxlPlane == null ? null : jvxlData.mapLattice == null ? ["plane", JU.Escape.eP4 (jvxlData.jvxlPlane)] : ["plane", JU.Escape.eP4 (jvxlData.jvxlPlane), "maplattice", JU.Escape.eP (jvxlData.mapLattice)]));
if (vertexDataOnly) {
J.jvxl.data.JvxlCoder.appendXmlVertexOnlyData (sb, jvxlData, meshData, true);
} else if (jvxlData.jvxlPlane == null) {
if (jvxlData.jvxlEdgeData == null) return "";
J.jvxl.data.JvxlCoder.appendXmlEdgeData (sb, jvxlData);
J.jvxl.data.JvxlCoder.appendXmlColorData (sb, jvxlData.jvxlColorData, true, jvxlData.isJvxlPrecisionColor, jvxlData.valueMappedToRed, jvxlData.valueMappedToBlue);
} else {
J.jvxl.data.JvxlCoder.appendXmlColorData (sb, jvxlData.jvxlColorData, true, jvxlData.isJvxlPrecisionColor, jvxlData.valueMappedToRed, jvxlData.valueMappedToBlue);
}J.jvxl.data.JvxlCoder.appendEncodedBitSetTag (sb, "jvxlInvalidatedVertexData", jvxlData.jvxlExcluded[1], -1, null);
if (jvxlData.excludedVertexCount > 0) {
J.jvxl.data.JvxlCoder.appendEncodedBitSetTag (sb, "jvxlExcludedVertexData", jvxlData.jvxlExcluded[0], jvxlData.excludedVertexCount, null);
J.jvxl.data.JvxlCoder.appendEncodedBitSetTag (sb, "jvxlExcludedPlaneData", jvxlData.jvxlExcluded[2], -1, null);
}J.jvxl.data.JvxlCoder.appendEncodedBitSetTag (sb, "jvxlExcludedTriangleData", jvxlData.jvxlExcluded[3], jvxlData.excludedTriangleCount, null);
JU.XmlUtil.closeTag (sb, "jvxlSurfaceData");
var len = sb.length ();
data.appendSB (sb);
if (jvxlData.vContours != null && jvxlData.vContours.length > 0) {
J.jvxl.data.JvxlCoder.jvxlEncodeContourData (jvxlData.vContours, data);
}if (jvxlData.vertexColorMap != null) {
if (jvxlData.baseColor == null) JU.XmlUtil.openTag (data, "jvxlVertexColorData");
 else JU.XmlUtil.openTagAttr (data, "jvxlVertexColorData", ["baseColor", jvxlData.baseColor]);
for (var entry, $entry = jvxlData.vertexColorMap.entrySet ().iterator (); $entry.hasNext () && ((entry = $entry.next ()) || true);) J.jvxl.data.JvxlCoder.appendEncodedBitSetTag (data, "jvxlColorMap", entry.getValue (), -1, ["color", entry.getKey ()]);

jvxlData.vertexColorMap = null;
JU.XmlUtil.closeTag (data, "jvxlVertexColorData");
}JU.XmlUtil.closeTag (data, "jvxlSurface");
if (includeHeader) {
JU.XmlUtil.closeTag (data, "jvxlSurfaceSet");
JU.XmlUtil.closeTag (data, "jvxl");
}return J.jvxl.data.JvxlCoder.jvxlSetCompressionRatio (data, jvxlData, len);
}, "J.jvxl.data.JvxlData,J.jvxl.data.MeshData,~A,~S,~B,~N,~S,~S");
c$.appendEncodedBitSetTag = Clazz_defineMethod (c$, "appendEncodedBitSetTag", 
 function (sb, name, bs, count, attribs) {
if (count < 0) count = JU.BSUtil.cardinalityOf (bs);
if (count == 0) return;
var sb1 =  new JU.SB ();
sb1.append ("\n ");
J.jvxl.data.JvxlCoder.jvxlEncodeBitSetBuffer (bs, -1, sb1);
JU.XmlUtil.appendTagObj (sb, name, [attribs, "bsEncoding", "base90+35", "count", "" + count, "len", "" + bs.length ()], J.jvxl.data.JvxlCoder.jvxlCompressString (sb1.toString (), true));
}, "JU.SB,~S,JU.BS,~N,~A");
c$.jvxlSetCompressionRatio = Clazz_defineMethod (c$, "jvxlSetCompressionRatio", 
 function (data, jvxlData, len) {
var s = data.toString ();
var r = Clazz_floatToInt (jvxlData.nBytes > 0 ? (jvxlData.nBytes) / len : ((jvxlData.nPointsX * jvxlData.nPointsY * jvxlData.nPointsZ * 13)) / len);
return JU.PT.rep (s, "\"not calculated\"", (r > 0 ? "\"" + r + ":1\"" : "\"?\""));
}, "JU.SB,J.jvxl.data.JvxlData,~N");
c$.appendXmlEdgeData = Clazz_defineMethod (c$, "appendXmlEdgeData", 
 function (sb, jvxlData) {
JU.XmlUtil.appendTagObj (sb, "jvxlEdgeData", ["count", "" + (jvxlData.jvxlEdgeData.length - 1), "encoding", "base90f1", "bsEncoding", "base90+35c", "isXLowToHigh", "" + jvxlData.isXLowToHigh, "data", J.jvxl.data.JvxlCoder.jvxlCompressString (jvxlData.jvxlEdgeData, true)], "\n" + J.jvxl.data.JvxlCoder.jvxlCompressString (jvxlData.jvxlSurfaceData, true));
}, "JU.SB,J.jvxl.data.JvxlData");
c$.jvxlAppendCommandState = Clazz_defineMethod (c$, "jvxlAppendCommandState", 
 function (data, cmd, state) {
if (cmd != null) JU.XmlUtil.appendCdata (data, "jvxlIsosurfaceCommand", null, "\n" + (cmd.indexOf ("#") < 0 ? cmd : cmd.substring (0, cmd.indexOf ("#"))) + "\n");
if (state != null) {
if (state.indexOf ("** XML ** ") >= 0) {
state = JU.PT.split (state, "** XML **")[1].trim ();
JU.XmlUtil.appendTag (data, "jvxlIsosurfaceState", "\n" + state + "\n");
} else {
JU.XmlUtil.appendCdata (data, "jvxlIsosurfaceState", null, "\n" + state);
}}}, "JU.SB,~S,~S");
c$.appendXmlColorData = Clazz_defineMethod (c$, "appendXmlColorData", 
 function (sb, data, isEncoded, isPrecisionColor, value1, value2) {
var n;
if (data == null || (n = data.length - 1) < 0) return;
if (isPrecisionColor) n /= 2;
JU.XmlUtil.appendTagObj (sb, "jvxlColorData", ["count", "" + n, "encoding", (isEncoded ? "base90f" + (isPrecisionColor ? "2" : "1") : "none"), "min", "" + value1, "max", "" + value2, "data", J.jvxl.data.JvxlCoder.jvxlCompressString (data, true)], null);
}, "JU.SB,~S,~B,~B,~N,~N");
c$.jvxlGetInfo = Clazz_defineMethod (c$, "jvxlGetInfo", 
function (jvxlData) {
return J.jvxl.data.JvxlCoder.jvxlGetInfoData (jvxlData, jvxlData.vertexDataOnly);
}, "J.jvxl.data.JvxlData");
c$.jvxlGetInfoData = Clazz_defineMethod (c$, "jvxlGetInfoData", 
function (jvxlData, vertexDataOnly) {
if (jvxlData.jvxlSurfaceData == null) return "";
var attribs =  new JU.Lst ();
var nSurfaceInts = jvxlData.nSurfaceInts;
var bytesUncompressedEdgeData = (vertexDataOnly ? 0 : jvxlData.jvxlEdgeData.length - 1);
var nColorData = (jvxlData.jvxlColorData == null ? -1 : (jvxlData.jvxlColorData.length - 1));
if (!vertexDataOnly) {
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  cutoff", "" + jvxlData.cutoff);
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  isCutoffAbsolute", "" + jvxlData.isCutoffAbsolute);
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  pointsPerAngstrom", "" + jvxlData.pointsPerAngstrom);
var n = jvxlData.jvxlSurfaceData.length + bytesUncompressedEdgeData + nColorData + 1;
if (n > 0) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  nBytesData", "" + n);
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  isXLowToHigh", "" + jvxlData.isXLowToHigh);
if (jvxlData.jvxlPlane == null) {
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  nSurfaceInts", "" + nSurfaceInts);
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  nBytesUncompressedEdgeData", "" + bytesUncompressedEdgeData);
}if (nColorData > 0) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  nBytesUncompressedColorData", "" + nColorData);
}jvxlData.excludedVertexCount = JU.BSUtil.cardinalityOf (jvxlData.jvxlExcluded[0]);
jvxlData.excludedTriangleCount = JU.BSUtil.cardinalityOf (jvxlData.jvxlExcluded[3]);
if (jvxlData.excludedVertexCount > 0) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  nExcludedVertexes", "" + jvxlData.excludedVertexCount);
if (jvxlData.excludedTriangleCount > 0) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  nExcludedTriangles", "" + jvxlData.excludedTriangleCount);
var n = JU.BSUtil.cardinalityOf (jvxlData.jvxlExcluded[1]);
if (n > 0) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  nInvalidatedVertexes", "" + n);
if (jvxlData.slabInfo != null) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  slabInfo", jvxlData.slabInfo);
if (jvxlData.isJvxlPrecisionColor) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  precisionColor", "true");
if (jvxlData.colorDensity) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  colorDensity", "true");
if (!Float.isNaN (jvxlData.pointSize)) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  pointSize", "" + jvxlData.pointSize);
 else if (jvxlData.diameter != 0) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  diameter", "" + jvxlData.diameter);
if (!jvxlData.allowVolumeRender) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  allowVolumeRender", "false");
if (jvxlData.jvxlPlane == null || vertexDataOnly) {
if (jvxlData.isContoured) {
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  contoured", "true");
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  colorMapped", "true");
} else if (jvxlData.isBicolorMap) {
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  bicolorMap", "true");
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  colorNegative", JU.C.getHexCode (jvxlData.minColorIndex));
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  colorPositive", JU.C.getHexCode (jvxlData.maxColorIndex));
} else if (nColorData > 0) {
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  colorMapped", "true");
}if (jvxlData.vContours != null && jvxlData.vContours.length > 0) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  nContourData", "" + jvxlData.vContours.length);
} else {
if (jvxlData.mapLattice != null) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  mapLattice", "" + jvxlData.mapLattice);
if (jvxlData.scale3d != 0) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  scale3d", "" + jvxlData.scale3d);
if (nColorData > 0) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  colorMapped", "true");
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  plane", JU.Escape.eP4 (jvxlData.jvxlPlane));
}if (jvxlData.color != null && jvxlData.color.indexOf ("null") < 0) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  color", jvxlData.color);
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  translucency", "" + jvxlData.translucency);
if (jvxlData.meshColor != null) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  meshColor", jvxlData.meshColor);
if (jvxlData.colorScheme != null) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  colorScheme", jvxlData.colorScheme);
if (jvxlData.rendering != null) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  rendering", jvxlData.rendering);
if (jvxlData.thisSet >= 0) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  set", "" + (jvxlData.thisSet + 1));
if (jvxlData.slabValue != -2147483648) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  slabValue", "" + jvxlData.slabValue);
if (jvxlData.isSlabbable) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  slabbable", "true");
if (jvxlData.nVertexColors > 0) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  nVertexColors", "" + jvxlData.nVertexColors);
var min = (jvxlData.mappedDataMin == 3.4028235E38 ? 0 : jvxlData.mappedDataMin);
var blue = (jvxlData.isColorReversed ? jvxlData.valueMappedToRed : jvxlData.valueMappedToBlue);
var red = (jvxlData.isColorReversed ? jvxlData.valueMappedToBlue : jvxlData.valueMappedToRed);
if (jvxlData.jvxlColorData != null && jvxlData.jvxlColorData.length > 0 && !jvxlData.isBicolorMap) {
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  dataMinimum", "" + min);
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  dataMaximum", "" + jvxlData.mappedDataMax);
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  valueMappedToRed", "" + red);
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  valueMappedToBlue", "" + blue);
}if (jvxlData.isContoured) {
if (jvxlData.contourValues == null || jvxlData.contourColixes == null) {
if (jvxlData.vContours == null) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  nContours", "" + Math.abs (jvxlData.nContours));
} else {
if (jvxlData.jvxlPlane != null) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  contoured", "true");
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  nContours", "" + jvxlData.contourValues.length);
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  contourValues", JU.Escape.eAF (jvxlData.contourValuesUsed == null ? jvxlData.contourValues : jvxlData.contourValuesUsed));
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  contourColors", jvxlData.contourColors);
}}if (jvxlData.insideOut) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  insideOut", "true");
if (jvxlData.vertexDataOnly) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  note", "vertex/face data only");
 else if (jvxlData.isXLowToHigh) J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  note", "progressive JVXL+ -- X values read from low(0) to high(" + (jvxlData.nPointsX - 1) + ")");
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  xyzMin", JU.Escape.eP (jvxlData.boundingBox[0]));
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  xyzMax", JU.Escape.eP (jvxlData.boundingBox[1]));
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  approximateCompressionRatio", "not calculated");
J.jvxl.data.JvxlCoder.addAttrib (attribs, "\n  jmolVersion", jvxlData.version);
var info =  new JU.SB ();
JU.XmlUtil.openTagAttr (info, "jvxlSurfaceInfo", attribs.toArray ( new Array (attribs.size ())));
JU.XmlUtil.closeTag (info, "jvxlSurfaceInfo");
return info.toString ();
}, "J.jvxl.data.JvxlData,~B");
c$.addAttrib = Clazz_defineMethod (c$, "addAttrib", 
 function (attribs, name, value) {
attribs.addLast ([name, value]);
}, "JU.Lst,~S,~S");
c$.jvxlEncodeContourData = Clazz_defineMethod (c$, "jvxlEncodeContourData", 
 function (contours, sb) {
JU.XmlUtil.openTagAttr (sb, "jvxlContourData", ["count", "" + contours.length]);
for (var i = 0; i < contours.length; i++) {
if (contours[i].size () < 6) {
continue;
}var nPolygons = (contours[i].get (0)).intValue ();
var sb1 =  new JU.SB ();
sb1.append ("\n");
var bs = contours[i].get (1);
J.jvxl.data.JvxlCoder.jvxlEncodeBitSetBuffer (bs, nPolygons, sb1);
JU.XmlUtil.appendTagObj (sb, "jvxlContour", ["index", "" + i, "value", "" + contours[i].get (2), "color", JU.Escape.escapeColor ((contours[i].get (4))[0]), "count", "" + bs.length (), "encoding", "base90iff1", "bsEncoding", "base90+35c", "data", J.jvxl.data.JvxlCoder.jvxlCompressString (contours[i].get (5).toString (), true)], J.jvxl.data.JvxlCoder.jvxlCompressString (sb1.toString (), true));
}
JU.XmlUtil.closeTag (sb, "jvxlContourData");
}, "~A,JU.SB");
c$.set3dContourVector = Clazz_defineMethod (c$, "set3dContourVector", 
function (v, polygonIndexes, vertices) {
if (v.size () < 6) return;
var fData = v.get (5);
var bs = v.get (1);
var pt = 0;
var nBuf = fData.length ();
var type = 0;
var c1 = ' ';
var c2 = ' ';
for (var i = bs.nextSetBit (0); i >= 0; i = bs.nextSetBit (i + 1)) {
var vertexIndexes = polygonIndexes[i];
while (pt < nBuf && !Character.isDigit (c1 = fData.charAt (pt++))) {
}
type = c1.charCodeAt (0) - 48;
while (pt < nBuf && Character.isWhitespace (c1 = fData.charAt (pt++))) {
}
while (pt < nBuf && Character.isWhitespace (c2 = fData.charAt (pt++))) {
}
var f1 = J.jvxl.data.JvxlCoder.jvxlFractionFromCharacter (c1.charCodeAt (0), 35, 90, 0);
var f2 = J.jvxl.data.JvxlCoder.jvxlFractionFromCharacter (c2.charCodeAt (0), 35, 90, 0);
var i1;
var i2;
var i3;
var i4;
if ((type & 1) == 0) {
i1 = vertexIndexes[1];
i2 = i3 = vertexIndexes[2];
i4 = vertexIndexes[0];
} else {
i1 = vertexIndexes[0];
i2 = vertexIndexes[1];
if ((type & 2) != 0) {
i3 = i2;
i4 = vertexIndexes[2];
} else {
i3 = vertexIndexes[2];
i4 = i1;
}}v.addLast (J.jvxl.data.JvxlCoder.getContourPoint (vertices, i1, i2, f1));
v.addLast (J.jvxl.data.JvxlCoder.getContourPoint (vertices, i3, i4, f2));
}
}, "JU.Lst,~A,~A");
c$.getContourPoint = Clazz_defineMethod (c$, "getContourPoint", 
 function (vertices, i, j, f) {
var pt =  new JU.P3 ();
pt.sub2 (vertices[j], vertices[i]);
pt.scaleAdd2 (f, pt, vertices[i]);
return pt;
}, "~A,~N,~N,~N");
c$.appendContourTriangleIntersection = Clazz_defineMethod (c$, "appendContourTriangleIntersection", 
function (type, f1, f2, fData) {
fData.appendI (type);
fData.appendC (J.jvxl.data.JvxlCoder.jvxlFractionAsCharacter (f1));
fData.appendC (J.jvxl.data.JvxlCoder.jvxlFractionAsCharacter (f2));
}, "~N,~N,~N,JU.SB");
c$.jvxlCreateColorData = Clazz_defineMethod (c$, "jvxlCreateColorData", 
function (jvxlData, vertexValues) {
if (vertexValues == null) {
jvxlData.jvxlColorData = "";
return;
}var writePrecisionColor = jvxlData.isJvxlPrecisionColor;
var doTruncate = jvxlData.isTruncated;
var colorFractionBase = jvxlData.colorFractionBase;
var colorFractionRange = jvxlData.colorFractionRange;
var valueBlue = jvxlData.valueMappedToBlue;
var valueRed = jvxlData.valueMappedToRed;
var vertexCount = (jvxlData.saveVertexCount > 0 ? jvxlData.saveVertexCount : jvxlData.vertexCount);
if (vertexCount > vertexValues.length) System.out.println ("JVXLCODER ERROR");
var min = jvxlData.mappedDataMin;
var max = jvxlData.mappedDataMax;
var list1 =  new JU.SB ();
var list2 =  new JU.SB ();
if (vertexValues.length < vertexCount) System.out.println ("JVXLCOLOR OHOHO");
for (var i = 0; i < vertexCount; i++) {
var value = vertexValues[i];
if (Float.isNaN (value)) value = min;
if (doTruncate) value = (value > 0 ? 0.999 : -0.999);
if (writePrecisionColor) J.jvxl.data.JvxlCoder.jvxlAppendCharacter2 (value, min, max, colorFractionBase, colorFractionRange, list1, list2);
 else list1.appendC (J.jvxl.data.JvxlCoder.jvxlValueAsCharacter (value, valueRed, valueBlue, colorFractionBase, colorFractionRange));
}
jvxlData.jvxlColorData = list1.appendSB (list2).appendC ('\n').toString ();
}, "J.jvxl.data.JvxlData,~A");
c$.appendXmlVertexOnlyData = Clazz_defineMethod (c$, "appendXmlVertexOnlyData", 
 function (sb, jvxlData, meshData, escapeXml) {
var vertexIdNew =  Clazz_newIntArray (meshData.vc, 0);
if (J.jvxl.data.JvxlCoder.appendXmlTriangleData (sb, meshData.pis, meshData.pc, meshData.bsSlabDisplay, vertexIdNew, escapeXml)) J.jvxl.data.JvxlCoder.appendXmlVertexData (sb, jvxlData, vertexIdNew, meshData.vs, meshData.vvs, meshData.vc, meshData.polygonColorData, meshData.pc, meshData.bsSlabDisplay, jvxlData.vertexColors, jvxlData.jvxlColorData.length > 0, escapeXml);
}, "JU.SB,J.jvxl.data.JvxlData,J.jvxl.data.MeshData,~B");
c$.appendXmlTriangleData = Clazz_defineMethod (c$, "appendXmlTriangleData", 
 function (sb, triangles, nData, bsSlabDisplay, vertexIdNew, escapeXml) {
var list1 =  new JU.SB ();
var list2 =  new JU.SB ();
var ilast = 1;
var p = 0;
var inew = 0;
var addPlus = false;
var nTri = 0;
var removeSlabbed = (bsSlabDisplay != null);
for (var i = 0; i < nData; ) {
if (triangles[i] == null || (removeSlabbed && !bsSlabDisplay.get (i))) {
i++;
continue;
}var idata = triangles[i][p];
if (vertexIdNew[idata] > 0) {
idata = vertexIdNew[idata];
} else {
idata = vertexIdNew[idata] = ++inew;
}var diff = idata - ilast;
ilast = idata;
if (diff == 0) {
list1.appendC ('!');
addPlus = false;
} else if (diff > 32) {
if (addPlus) list1.appendC ('+');
list1.appendI (diff);
addPlus = true;
} else if (diff < -32) {
list1.appendI (diff);
addPlus = true;
} else {
list1.appendC (String.fromCharCode (92 + diff));
addPlus = false;
}if (++p % 3 == 0) {
list2.appendI (triangles[i][3]);
p = 0;
i++;
nTri++;
}}
if (list1.length () == 0) return true;
JU.XmlUtil.appendTagObj (sb, "jvxlTriangleData", ["count", "" + nTri, "encoding", "jvxltdiff", "data", J.jvxl.data.JvxlCoder.jvxlCompressString (list1.toString (), escapeXml)], null);
JU.XmlUtil.appendTagObj (sb, "jvxlTriangleEdgeData", ["count", "" + nTri, "encoding", "jvxlsc", "data", J.jvxl.data.JvxlCoder.jvxlCompressString (list2.toString (), escapeXml)], null);
return true;
}, "JU.SB,~A,~N,JU.BS,~A,~B");
c$.appendXmlVertexData = Clazz_defineMethod (c$, "appendXmlVertexData", 
 function (sb, jvxlData, vertexIdNew, vertices, vertexValues, vertexCount, polygonColorData, polygonCount, bsSlabDisplay, vertexColors, addColorData, escapeXml) {
var colorFractionBase = jvxlData.colorFractionBase;
var colorFractionRange = jvxlData.colorFractionRange;
var p;
var min = jvxlData.boundingBox[0];
var max = jvxlData.boundingBox[1];
var list1 =  new JU.SB ();
var list2 =  new JU.SB ();
var vertexIdOld = null;
var removeSlabbed = (bsSlabDisplay != null);
if (polygonCount > 0) {
if (removeSlabbed) polygonCount = bsSlabDisplay.cardinality ();
removeSlabbed = false;
vertexIdOld =  Clazz_newIntArray (vertexCount, 0);
for (var i = 0; i < vertexCount; i++) if (vertexIdNew[i] > 0) vertexIdOld[vertexIdNew[i] - 1] = i;

}var n = 0;
for (var i = 0; i < vertexCount; i++) if (!removeSlabbed || bsSlabDisplay.get (i)) {
n++;
p = vertices[(polygonCount == 0 ? i : vertexIdOld[i])];
J.jvxl.data.JvxlCoder.jvxlAppendCharacter2 (p.x, min.x, max.x, colorFractionBase, colorFractionRange, list1, list2);
J.jvxl.data.JvxlCoder.jvxlAppendCharacter2 (p.y, min.y, max.y, colorFractionBase, colorFractionRange, list1, list2);
J.jvxl.data.JvxlCoder.jvxlAppendCharacter2 (p.z, min.z, max.z, colorFractionBase, colorFractionRange, list1, list2);
}
list1.appendSB (list2);
JU.XmlUtil.appendTagObj (sb, "jvxlVertexData", ["count", "" + n, "min", JU.Escape.eP (min), "max", JU.Escape.eP (max), "encoding", "base90xyz2", "data", J.jvxl.data.JvxlCoder.jvxlCompressString (list1.toString (), escapeXml)], null);
if (polygonColorData != null) JU.XmlUtil.appendTagObj (sb, "jvxlPolygonColorData", ["encoding", "jvxlnc", "count", "" + polygonCount], "\n" + polygonColorData);
if (!addColorData) return;
list1 =  new JU.SB ();
list2 =  new JU.SB ();
if (vertexColors == null) {
for (var i = 0; i < vertexCount; i++) if (!removeSlabbed || bsSlabDisplay.get (i)) {
var value = vertexValues[polygonCount == 0 ? i : vertexIdOld[i]];
J.jvxl.data.JvxlCoder.jvxlAppendCharacter2 (value, jvxlData.mappedDataMin, jvxlData.mappedDataMax, colorFractionBase, colorFractionRange, list1, list2);
}
} else {
var lastColor = 0;
list1.appendI (n).append (" ");
for (var i = 0; i < vertexCount; i++) if (!removeSlabbed || bsSlabDisplay.get (i)) {
var c = vertexColors[polygonCount == 0 ? i : vertexIdOld[i]];
if (c == lastColor) c = 0;
 else lastColor = c;
list1.appendI (c);
list1.append (" ");
}
}J.jvxl.data.JvxlCoder.appendXmlColorData (sb, list1.appendSB (list2).append ("\n").toString (), (vertexColors == null), true, jvxlData.valueMappedToRed, jvxlData.valueMappedToBlue);
}, "JU.SB,J.jvxl.data.JvxlData,~A,~A,~A,~N,~S,~N,JU.BS,~A,~B,~B");
c$.jvxlFractionAsCharacter = Clazz_defineMethod (c$, "jvxlFractionAsCharacter", 
function (fraction) {
return J.jvxl.data.JvxlCoder.jvxlFractionAsCharacterRange (fraction, 35, 90);
}, "~N");
c$.jvxlFractionAsCharacterRange = Clazz_defineMethod (c$, "jvxlFractionAsCharacterRange", 
function (fraction, base, range) {
if (fraction > 0.9999) fraction = 0.9999;
 else if (Float.isNaN (fraction)) fraction = 1.0001;
var ich = Clazz_doubleToInt (Math.floor (fraction * range + base));
if (ich < base) return String.fromCharCode (base);
if (ich == 92) return '!';
return String.fromCharCode (ich);
}, "~N,~N,~N");
c$.jvxlAppendCharacter2 = Clazz_defineMethod (c$, "jvxlAppendCharacter2", 
 function (value, min, max, base, range, list1, list2) {
var fraction = (min == max ? value : (value - min) / (max - min));
var ch1 = J.jvxl.data.JvxlCoder.jvxlFractionAsCharacterRange (fraction, base, range);
list1.appendC (ch1);
fraction -= J.jvxl.data.JvxlCoder.jvxlFractionFromCharacter (ch1.charCodeAt (0), base, range, 0);
list2.appendC (J.jvxl.data.JvxlCoder.jvxlFractionAsCharacterRange (fraction * range, base, range));
}, "~N,~N,~N,~N,~N,JU.SB,JU.SB");
c$.jvxlFractionFromCharacter = Clazz_defineMethod (c$, "jvxlFractionFromCharacter", 
function (ich, base, range, fracOffset) {
if (ich == base + range) return NaN;
if (ich < base) ich = 92;
var fraction = (ich - base + fracOffset) / range;
if (fraction < 0) return 0;
if (fraction > 1) return 0.999999;
return fraction;
}, "~N,~N,~N,~N");
c$.jvxlFractionFromCharacter2 = Clazz_defineMethod (c$, "jvxlFractionFromCharacter2", 
function (ich1, ich2, base, range) {
var fraction = J.jvxl.data.JvxlCoder.jvxlFractionFromCharacter (ich1, base, range, 0);
var remains = J.jvxl.data.JvxlCoder.jvxlFractionFromCharacter (ich2, base, range, 0.5);
return fraction + remains / range;
}, "~N,~N,~N,~N");
c$.jvxlValueAsCharacter = Clazz_defineMethod (c$, "jvxlValueAsCharacter", 
function (value, min, max, base, range) {
var fraction = (min == max ? value : (value - min) / (max - min));
return J.jvxl.data.JvxlCoder.jvxlFractionAsCharacterRange (fraction, base, range);
}, "~N,~N,~N,~N,~N");
c$.jvxlValueFromCharacter2 = Clazz_defineMethod (c$, "jvxlValueFromCharacter2", 
function (ich, ich2, min, max, base, range) {
var fraction = J.jvxl.data.JvxlCoder.jvxlFractionFromCharacter2 (ich, ich2, base, range);
return (max == min ? fraction : min + fraction * (max - min));
}, "~N,~N,~N,~N,~N,~N");
c$.jvxlEncodeBitSet0 = Clazz_defineMethod (c$, "jvxlEncodeBitSet0", 
function (bs, nPoints, sb) {
var dataCount = 0;
var prevCount = -1;
var nPrev = 0;
if (nPoints < 0) nPoints = bs.length ();
var n = 0;
var isset = false;
var lastPoint = nPoints - 1;
for (var i = 0; i < nPoints; ++i) {
if (isset == bs.get (i)) {
dataCount++;
} else {
if (dataCount == prevCount && i != lastPoint) {
nPrev++;
} else {
if (nPrev > 0) {
sb.appendC (' ').appendI (-nPrev);
nPrev = 0;
n++;
}sb.appendC (' ').appendI (dataCount);
n++;
prevCount = dataCount;
}dataCount = 1;
isset = !isset;
}}
sb.appendC (' ').appendI (dataCount).appendC ('\n');
return n;
}, "JU.BS,~N,JU.SB");
c$.jvxlEncodeBitSet = Clazz_defineMethod (c$, "jvxlEncodeBitSet", 
function (bs) {
var sb =  new JU.SB ();
J.jvxl.data.JvxlCoder.jvxlEncodeBitSetBuffer (bs, -1, sb);
return sb.toString ();
}, "JU.BS");
c$.jvxlEncodeBitSetBuffer = Clazz_defineMethod (c$, "jvxlEncodeBitSetBuffer", 
function (bs, nPoints, sb) {
var dataCount = 0;
var n = 0;
var isset = false;
if (nPoints < 0) nPoints = bs.length ();
if (nPoints == 0) return 0;
sb.append ("-");
for (var i = 0; i < nPoints; ++i) {
if (isset == bs.get (i)) {
dataCount++;
} else {
J.jvxl.data.JvxlCoder.jvxlAppendEncodedNumber (sb, dataCount, 35, 90);
n++;
dataCount = 1;
isset = !isset;
}}
J.jvxl.data.JvxlCoder.jvxlAppendEncodedNumber (sb, dataCount, 35, 90);
sb.appendC ('\n');
return n;
}, "JU.BS,~N,JU.SB");
c$.jvxlAppendEncodedNumber = Clazz_defineMethod (c$, "jvxlAppendEncodedNumber", 
function (sb, n, base, range) {
var isInRange = (n < range);
if (n == 0) sb.appendC (String.fromCharCode (base));
 else if (!isInRange) sb.appendC (String.fromCharCode (base + range));
while (n > 0) {
var n1 = Clazz_doubleToInt (n / range);
var x = base + n - n1 * range;
if (x == 92) x = 33;
sb.appendC (String.fromCharCode (x));
n = n1;
}
if (!isInRange) sb.append (" ");
}, "JU.SB,~N,~N,~N");
c$.jvxlDecodeBitSetRange = Clazz_defineMethod (c$, "jvxlDecodeBitSetRange", 
function (data, base, range) {
var bs =  new JU.BS ();
var dataCount = 0;
var ptr = 0;
var isset = false;
var next =  Clazz_newIntArray (1, 0);
while ((dataCount = J.jvxl.data.JvxlCoder.jvxlParseEncodedInt (data, base, range, next)) != -2147483648) {
if (isset) bs.setBits (ptr, ptr + dataCount);
ptr += dataCount;
isset = !isset;
}
return bs;
}, "~S,~N,~N");
c$.jvxlParseEncodedInt = Clazz_defineMethod (c$, "jvxlParseEncodedInt", 
function (str, offset, base, next) {
var digitSeen = false;
var value = 0;
var ich = next[0];
var ichMax = str.length;
if (ich < 0) return -2147483648;
while (ich < ichMax && Character.isWhitespace (str.charAt (ich))) ++ich;

if (ich >= ichMax) return -2147483648;
var factor = 1;
var isLong = (str.charCodeAt (ich) == (offset + base));
if (isLong) ich++;
while (ich < ichMax && !Character.isWhitespace (str.charAt (ich))) {
var i = str.charCodeAt (ich);
if (i < offset) i = 92;
value += (i - offset) * factor;
digitSeen = true;
++ich;
if (!isLong) break;
factor *= base;
}
if (!digitSeen) value = -2147483648;
next[0] = ich;
return value;
}, "~S,~N,~N,~A");
c$.jvxlDecodeBitSet = Clazz_defineMethod (c$, "jvxlDecodeBitSet", 
function (data) {
if (data.startsWith ("-")) return J.jvxl.data.JvxlCoder.jvxlDecodeBitSetRange (J.jvxl.data.JvxlCoder.jvxlDecompressString (data.substring (1)), 35, 90);
var bs =  new JU.BS ();
var dataCount = 0;
var lastCount = 0;
var nPrev = 0;
var ptr = 0;
var isset = false;
var next =  Clazz_newIntArray (1, 0);
while (true) {
dataCount = (nPrev++ < 0 ? dataCount : JU.PT.parseIntNext (data, next));
if (dataCount == -2147483648) break;
if (dataCount < 0) {
nPrev = dataCount;
dataCount = lastCount;
continue;
}if (isset) bs.setBits (ptr, ptr + dataCount);
ptr += dataCount;
lastCount = dataCount;
isset = !isset;
}
return bs;
}, "~S");
c$.jvxlCompressString = Clazz_defineMethod (c$, "jvxlCompressString", 
function (data, escapeXml) {
if (data.indexOf ("~") >= 0) return data;
var dataOut =  new JU.SB ();
var chLast = '\u0000';
var escaped = false;
var lastEscaped = false;
var nLast = 0;
var n = data.length;
for (var i = 0; i <= n; i++) {
var ch = (i == n ? '\0' : data.charAt (i));
switch (ch) {
case '\n':
case '\r':
continue;
case '&':
case '<':
escaped = escapeXml;
break;
default:
escaped = false;
}
if (ch == chLast) {
++nLast;
ch = '\0';
} else if (nLast > 0 || lastEscaped) {
if (nLast < 4 && !lastEscaped || chLast == ' ' || chLast == '\t') {
while (--nLast >= 0) dataOut.appendC (chLast);

} else {
if (lastEscaped) lastEscaped = false;
 else dataOut.appendC ('~');
dataOut.appendI (nLast);
dataOut.appendC (' ');
}nLast = 0;
}if (ch != '\0') {
if (escaped) {
lastEscaped = true;
escaped = false;
dataOut.appendC ('~');
chLast = ch;
ch = String.fromCharCode (ch.charCodeAt (0) - 1);
} else {
chLast = ch;
}dataOut.appendC (ch);
}}
return dataOut.toString ();
}, "~S,~B");
c$.jvxlDecompressString = Clazz_defineMethod (c$, "jvxlDecompressString", 
function (data) {
if (data.indexOf ("~") < 0) return data;
var dataOut =  new JU.SB ();
var chLast = '\u0000';
var next =  Clazz_newIntArray (1, 0);
for (var i = 0; i < data.length; i++) {
var ch = data.charAt (i);
if (ch == '~') {
next[0] = ++i;
switch (ch = data.charAt (i)) {
case ';':
case '%':
next[0]++;
dataOut.appendC (chLast = (ch = String.fromCharCode (ch.charCodeAt (0) + 1)));
case '1':
case '2':
case '3':
case '4':
case '5':
case '6':
case '7':
case '8':
case '9':
var nChar = JU.PT.parseIntNext (data, next);
for (var c = 0; c < nChar; c++) dataOut.appendC (chLast);

i = next[0];
continue;
case '~':
--i;
break;
default:
JU.Logger.error ("Error uncompressing string " + data.substring (0, i) + "?");
}
}dataOut.appendC (ch);
chLast = ch;
}
return dataOut.toString ();
}, "~S");
c$.jvxlCreateHeaderWithoutTitleOrAtoms = Clazz_defineMethod (c$, "jvxlCreateHeaderWithoutTitleOrAtoms", 
function (v, bs) {
J.jvxl.data.JvxlCoder.jvxlCreateHeader (v, bs);
}, "J.jvxl.data.VolumeData,JU.SB");
c$.jvxlCreateHeader = Clazz_defineMethod (c$, "jvxlCreateHeader", 
function (v, sb) {
v.setVolumetricXml ();
if (sb.length () == 0) sb.append ("Line 1\nLine 2\n");
}, "J.jvxl.data.VolumeData,JU.SB");
Clazz_defineStatics (c$,
"JVXL_VERSION1", "2.0",
"JVXL_VERSION_XML", "2.3",
"CONTOUR_NPOLYGONS", 0,
"CONTOUR_BITSET", 1,
"CONTOUR_VALUE", 2,
"CONTOUR_COLIX", 3,
"CONTOUR_COLOR", 4,
"CONTOUR_FDATA", 5,
"CONTOUR_POINTS", 6,
"defaultEdgeFractionBase", 35,
"defaultEdgeFractionRange", 90,
"defaultColorFractionBase", 35,
"defaultColorFractionRange", 90);
});
Clazz_declarePackage ("J.jvxl.data");
Clazz_load (["J.api.VolumeDataInterface", "JU.M3", "$.P3", "$.V3"], "J.jvxl.data.VolumeData", ["java.lang.Float", "java.util.Hashtable", "JU.SB", "JU.Escape", "$.Logger"], function () {
c$ = Clazz_decorateAsClass (function () {
this.sr = null;
this.doIterate = true;
this.volumetricOrigin = null;
this.origin = null;
this.volumetricVectors = null;
this.voxelCounts = null;
this.nPoints = 0;
this.voxelData = null;
this.voxelMap = null;
this.volumetricVectorLengths = null;
this.maxVectorLength = 0;
this.minToPlaneDistance = 0;
this.yzCount = 0;
this.unitVolumetricVectors = null;
this.volumetricMatrix = null;
this.inverseMatrix = null;
this.thePlane = null;
this.thePlaneNormalMag = 0;
this.ptXyzTemp = null;
this.xmlData = null;
this.mappingPlane = null;
this.mappingPlaneNormalMag = 0;
this.minGrid = 0;
this.maxGrid = 0;
this.voxelVolume = 0;
this.spanningVectors = null;
this.isPeriodic = false;
this.isSquared = false;
this.edgeVector = null;
this.ptTemp = null;
Clazz_instantialize (this, arguments);
}, J.jvxl.data, "VolumeData", null, J.api.VolumeDataInterface);
Clazz_prepareFields (c$, function () {
this.volumetricOrigin =  new JU.P3 ();
this.origin =  Clazz_newFloatArray (3, 0);
this.volumetricVectors =  new Array (3);
this.voxelCounts =  Clazz_newIntArray (3, 0);
this.volumetricVectorLengths =  Clazz_newFloatArray (3, 0);
this.unitVolumetricVectors =  new Array (3);
this.volumetricMatrix =  new JU.M3 ();
this.inverseMatrix =  new JU.M3 ();
this.ptXyzTemp =  new JU.P3 ();
this.edgeVector =  new JU.V3 ();
this.ptTemp =  new JU.P3 ();
});
Clazz_overrideMethod (c$, "getVoxelData", 
function () {
return this.voxelData;
});
Clazz_overrideMethod (c$, "setVoxelDataAsArray", 
function (voxelData) {
this.voxelData = voxelData;
if (voxelData != null) this.sr = null;
}, "~A");
Clazz_defineMethod (c$, "hasPlane", 
function () {
return (this.thePlane != null);
});
Clazz_makeConstructor (c$, 
function () {
this.volumetricVectors[0] =  new JU.V3 ();
this.volumetricVectors[1] =  new JU.V3 ();
this.volumetricVectors[2] =  new JU.V3 ();
this.unitVolumetricVectors[0] =  new JU.V3 ();
this.unitVolumetricVectors[1] =  new JU.V3 ();
this.unitVolumetricVectors[2] =  new JU.V3 ();
});
Clazz_defineMethod (c$, "setMappingPlane", 
function (plane) {
this.mappingPlane = plane;
if (plane == null) return;
this.mappingPlaneNormalMag = Math.sqrt (plane.x * plane.x + plane.y * plane.y + plane.z * plane.z);
}, "JU.P4");
Clazz_defineMethod (c$, "distanceToMappingPlane", 
function (pt) {
return (this.mappingPlane.x * pt.x + this.mappingPlane.y * pt.y + this.mappingPlane.z * pt.z + this.mappingPlane.w) / this.mappingPlaneNormalMag;
}, "JU.T3");
Clazz_overrideMethod (c$, "setVolumetricOrigin", 
function (x, y, z) {
this.volumetricOrigin.set (x, y, z);
}, "~N,~N,~N");
Clazz_overrideMethod (c$, "getOriginFloat", 
function () {
return this.origin;
});
Clazz_defineMethod (c$, "getSpanningVectors", 
function () {
return this.spanningVectors;
});
Clazz_defineMethod (c$, "getYzCount", 
function () {
this.minGrid = this.volumetricVectors[0].length ();
this.minGrid = Math.min (this.minGrid, this.volumetricVectors[1].length ());
this.minGrid = Math.min (this.minGrid, this.volumetricVectors[2].length ());
this.maxGrid = this.volumetricVectors[0].length ();
this.maxGrid = Math.max (this.maxGrid, this.volumetricVectors[1].length ());
this.maxGrid = Math.max (this.maxGrid, this.volumetricVectors[2].length ());
this.nPoints = this.voxelCounts[0] * this.voxelCounts[1] * this.voxelCounts[2];
return this.yzCount = this.voxelCounts[1] * this.voxelCounts[2];
});
Clazz_overrideMethod (c$, "getVolumetricVectorLengths", 
function () {
return this.volumetricVectorLengths;
});
Clazz_overrideMethod (c$, "setVolumetricVector", 
function (i, x, y, z) {
this.volumetricVectors[i].x = x;
this.volumetricVectors[i].y = y;
this.volumetricVectors[i].z = z;
this.setUnitVectors ();
}, "~N,~N,~N,~N");
Clazz_overrideMethod (c$, "getVoxelCounts", 
function () {
return this.voxelCounts;
});
Clazz_overrideMethod (c$, "setVoxelCounts", 
function (nPointsX, nPointsY, nPointsZ) {
this.voxelCounts[0] = nPointsX;
this.voxelCounts[1] = nPointsY;
this.voxelCounts[2] = nPointsZ;
return nPointsX * nPointsY * nPointsZ;
}, "~N,~N,~N");
Clazz_defineMethod (c$, "getVoxelDataAt", 
function (pt) {
var ix = Clazz_doubleToInt (pt / this.yzCount);
pt -= ix * this.yzCount;
var iy = Clazz_doubleToInt (pt / this.voxelCounts[2]);
var iz = pt - iy * this.voxelCounts[2];
return this.voxelData[ix][iy][iz];
}, "~N");
Clazz_defineMethod (c$, "getPointIndex", 
function (x, y, z) {
return x * this.yzCount + y * this.voxelCounts[2] + z;
}, "~N,~N,~N");
Clazz_defineMethod (c$, "getPoint", 
function (ipt, pt) {
var ix = Clazz_doubleToInt (ipt / this.yzCount);
ipt -= ix * this.yzCount;
var iy = Clazz_doubleToInt (ipt / this.voxelCounts[2]);
var iz = ipt - iy * this.voxelCounts[2];
this.voxelPtToXYZ (ix, iy, iz, pt);
}, "~N,JU.P3");
Clazz_defineMethod (c$, "setVoxelData", 
function (pt, value) {
var ix = Clazz_doubleToInt (pt / this.yzCount);
pt -= ix * this.yzCount;
var iy = Clazz_doubleToInt (pt / this.voxelCounts[2]);
var iz = pt - iy * this.voxelCounts[2];
this.voxelData[ix][iy][iz] = value;
}, "~N,~N");
Clazz_defineMethod (c$, "setVoxelMap", 
function () {
this.voxelMap =  new java.util.Hashtable ();
this.getYzCount ();
});
Clazz_defineMethod (c$, "setMatrix", 
 function () {
for (var i = 0; i < 3; i++) this.volumetricMatrix.setColumnV (i, this.volumetricVectors[i]);

try {
this.inverseMatrix.invertM (this.volumetricMatrix);
} catch (e) {
if (Clazz_exceptionOf (e, Exception)) {
JU.Logger.error ("VolumeData error setting matrix -- bad unit vectors? ");
return false;
} else {
throw e;
}
}
return true;
});
Clazz_overrideMethod (c$, "transform", 
function (v1, v2) {
this.volumetricMatrix.rotate2 (v1, v2);
}, "JU.V3,JU.V3");
Clazz_overrideMethod (c$, "setPlaneParameters", 
function (plane) {
this.thePlane = plane;
this.thePlaneNormalMag = Math.sqrt (plane.x * plane.x + plane.y * plane.y + plane.z * plane.z);
}, "JU.P4");
Clazz_overrideMethod (c$, "calcVoxelPlaneDistance", 
function (x, y, z) {
this.voxelPtToXYZ (x, y, z, this.ptXyzTemp);
return this.distancePointToPlane (this.ptXyzTemp);
}, "~N,~N,~N");
Clazz_defineMethod (c$, "getToPlaneParameter", 
function () {
return (Math.sqrt (this.thePlane.x * this.thePlane.x + this.thePlane.y * this.thePlane.y + this.thePlane.z * this.thePlane.z) * this.minToPlaneDistance);
});
Clazz_defineMethod (c$, "isNearPlane", 
function (x, y, z, toPlaneParameter) {
this.voxelPtToXYZ (x, y, z, this.ptXyzTemp);
return ((this.thePlane.x * this.ptXyzTemp.x + this.thePlane.y * this.ptXyzTemp.y + this.thePlane.z * this.ptXyzTemp.z + this.thePlane.w) < toPlaneParameter);
}, "~N,~N,~N,~N");
Clazz_overrideMethod (c$, "distancePointToPlane", 
function (pt) {
return (this.thePlane.x * pt.x + this.thePlane.y * pt.y + this.thePlane.z * pt.z + this.thePlane.w) / this.thePlaneNormalMag;
}, "JU.T3");
Clazz_overrideMethod (c$, "voxelPtToXYZ", 
function (x, y, z, pt) {
pt.scaleAdd2 (x, this.volumetricVectors[0], this.volumetricOrigin);
pt.scaleAdd2 (y, this.volumetricVectors[1], pt);
pt.scaleAdd2 (z, this.volumetricVectors[2], pt);
}, "~N,~N,~N,JU.T3");
Clazz_defineMethod (c$, "setUnitVectors", 
function () {
this.maxVectorLength = 0;
this.voxelVolume = 1;
for (var i = 0; i < 3; i++) {
var d = this.volumetricVectorLengths[i] = this.volumetricVectors[i].length ();
if (d == 0) return false;
if (d > this.maxVectorLength) this.maxVectorLength = d;
this.voxelVolume *= d;
this.unitVolumetricVectors[i].setT (this.volumetricVectors[i]);
this.unitVolumetricVectors[i].normalize ();
}
this.minToPlaneDistance = this.maxVectorLength * 2;
this.origin[0] = this.volumetricOrigin.x;
this.origin[1] = this.volumetricOrigin.y;
this.origin[2] = this.volumetricOrigin.z;
this.spanningVectors =  new Array (4);
this.spanningVectors[0] = JU.V3.newV (this.volumetricOrigin);
for (var i = 0; i < 3; i++) {
var v = this.spanningVectors[i + 1] =  new JU.V3 ();
v.scaleAdd2 (this.voxelCounts[i] - 1, this.volumetricVectors[i], v);
}
return this.setMatrix ();
});
Clazz_overrideMethod (c$, "xyzToVoxelPt", 
function (x, y, z, pt3i) {
this.ptXyzTemp.set (x, y, z);
this.ptXyzTemp.sub (this.volumetricOrigin);
this.inverseMatrix.rotate (this.ptXyzTemp);
pt3i.set (Math.round (this.ptXyzTemp.x), Math.round (this.ptXyzTemp.y), Math.round (this.ptXyzTemp.z));
}, "~N,~N,~N,JU.T3i");
Clazz_overrideMethod (c$, "lookupInterpolatedVoxelValue", 
function (point, getSource) {
if (this.mappingPlane != null) return this.distanceToMappingPlane (point);
if (this.sr != null) {
var v = this.sr.getValueAtPoint (point, getSource);
return (this.isSquared ? v * v : v);
}this.ptXyzTemp.sub2 (point, this.volumetricOrigin);
this.inverseMatrix.rotate (this.ptXyzTemp);
var iMax;
var xLower = this.indexLower (this.ptXyzTemp.x, iMax = this.voxelCounts[0] - 1);
var xUpper = this.indexUpper (this.ptXyzTemp.x, xLower, iMax);
var yLower = this.indexLower (this.ptXyzTemp.y, iMax = this.voxelCounts[1] - 1);
var yUpper = this.indexUpper (this.ptXyzTemp.y, yLower, iMax);
var zLower = this.indexLower (this.ptXyzTemp.z, iMax = this.voxelCounts[2] - 1);
var zUpper = this.indexUpper (this.ptXyzTemp.z, zLower, iMax);
var v1 = J.jvxl.data.VolumeData.getFractional2DValue (this.mantissa (this.ptXyzTemp.x - xLower), this.mantissa (this.ptXyzTemp.y - yLower), this.getVoxelValue (xLower, yLower, zLower), this.getVoxelValue (xUpper, yLower, zLower), this.getVoxelValue (xLower, yUpper, zLower), this.getVoxelValue (xUpper, yUpper, zLower));
var v2 = J.jvxl.data.VolumeData.getFractional2DValue (this.mantissa (this.ptXyzTemp.x - xLower), this.mantissa (this.ptXyzTemp.y - yLower), this.getVoxelValue (xLower, yLower, zUpper), this.getVoxelValue (xUpper, yLower, zUpper), this.getVoxelValue (xLower, yUpper, zUpper), this.getVoxelValue (xUpper, yUpper, zUpper));
return v1 + this.mantissa (this.ptXyzTemp.z - zLower) * (v2 - v1);
}, "JU.T3,~B");
Clazz_defineMethod (c$, "mantissa", 
 function (f) {
return (this.isPeriodic ? f - Math.floor (f) : f);
}, "~N");
Clazz_defineMethod (c$, "getVoxelValue", 
function (x, y, z) {
if (this.voxelMap == null) return this.voxelData[x][y][z];
var f = this.voxelMap.get (Integer.$valueOf (this.getPointIndex (x, y, z)));
return (f == null ? NaN : f.floatValue ());
}, "~N,~N,~N");
c$.getFractional2DValue = Clazz_defineMethod (c$, "getFractional2DValue", 
function (fx, fy, x11, x12, x21, x22) {
var v1 = x11 + fx * (x12 - x11);
var v2 = x21 + fx * (x22 - x21);
return v1 + fy * (v2 - v1);
}, "~N,~N,~N,~N,~N,~N");
Clazz_defineMethod (c$, "indexLower", 
 function (x, xMax) {
if (this.isPeriodic && xMax > 0) {
while (x < 0) x += xMax;

while (x >= xMax) x -= xMax;

return Clazz_doubleToInt (Math.floor (x));
}if (x < 0) return 0;
var floor = Clazz_doubleToInt (Math.floor (x));
return (floor > xMax ? xMax : floor);
}, "~N,~N");
Clazz_defineMethod (c$, "indexUpper", 
 function (x, xLower, xMax) {
return (!this.isPeriodic && x < 0 || xLower == xMax ? xLower : xLower + 1);
}, "~N,~N,~N");
Clazz_defineMethod (c$, "offsetCenter", 
function (center) {
var pt =  new JU.P3 ();
pt.scaleAdd2 ((this.voxelCounts[0] - 1) / 2, this.volumetricVectors[0], pt);
pt.scaleAdd2 ((this.voxelCounts[1] - 1) / 2, this.volumetricVectors[1], pt);
pt.scaleAdd2 ((this.voxelCounts[2] - 1) / 2, this.volumetricVectors[2], pt);
this.volumetricOrigin.sub2 (center, pt);
}, "JU.P3");
Clazz_overrideMethod (c$, "setDataDistanceToPlane", 
function (plane) {
this.setPlaneParameters (plane);
var nx = this.voxelCounts[0];
var ny = this.voxelCounts[1];
var nz = this.voxelCounts[2];
this.voxelData =  Clazz_newFloatArray (nx, ny, nz, 0);
for (var x = 0; x < nx; x++) for (var y = 0; y < ny; y++) for (var z = 0; z < nz; z++) this.voxelData[x][y][z] = this.calcVoxelPlaneDistance (x, y, z);



}, "JU.P4");
Clazz_overrideMethod (c$, "filterData", 
function (isSquared, invertCutoff) {
var doInvert = (!Float.isNaN (invertCutoff));
if (this.sr != null) {
this.isSquared = isSquared;
return;
}var nx = this.voxelCounts[0];
var ny = this.voxelCounts[1];
var nz = this.voxelCounts[2];
if (isSquared) for (var x = 0; x < nx; x++) for (var y = 0; y < ny; y++) for (var z = 0; z < nz; z++) this.voxelData[x][y][z] *= this.voxelData[x][y][z];



if (doInvert) for (var x = 0; x < nx; x++) for (var y = 0; y < ny; y++) for (var z = 0; z < nz; z++) this.voxelData[x][y][z] = invertCutoff - this.voxelData[x][y][z];



}, "~B,~N");
Clazz_overrideMethod (c$, "capData", 
function (plane, cutoff) {
if (this.voxelData == null) return;
var nx = this.voxelCounts[0];
var ny = this.voxelCounts[1];
var nz = this.voxelCounts[2];
var normal = JU.V3.new3 (plane.x, plane.y, plane.z);
normal.normalize ();
var f = 1;
for (var x = 0; x < nx; x++) for (var y = 0; y < ny; y++) for (var z = 0; z < nz; z++) {
var value = this.voxelData[x][y][z] - cutoff;
this.voxelPtToXYZ (x, y, z, this.ptXyzTemp);
var d = (this.ptXyzTemp.x * normal.x + this.ptXyzTemp.y * normal.y + this.ptXyzTemp.z * normal.z + plane.w - cutoff) / f;
if (d >= 0 || d > value) this.voxelData[x][y][z] = d;
}


}, "JU.P4,~N");
Clazz_defineMethod (c$, "setVolumetricXml", 
function () {
var sb =  new JU.SB ();
if (this.voxelCounts[0] == 0) {
sb.append ("<jvxlVolumeData>\n");
} else {
sb.append ("<jvxlVolumeData origin=\"" + JU.Escape.eP (this.volumetricOrigin) + "\">\n");
for (var i = 0; i < 3; i++) sb.append ("<jvxlVolumeVector type=\"" + i + "\" count=\"" + this.voxelCounts[i] + "\" vector=\"" + JU.Escape.eP (this.volumetricVectors[i]) + "\"></jvxlVolumeVector>\n");

}sb.append ("</jvxlVolumeData>\n");
return this.xmlData = sb.toString ();
});
Clazz_defineMethod (c$, "setVoxelMapValue", 
function (x, y, z, v) {
if (this.voxelMap == null) return;
this.voxelMap.put (Integer.$valueOf (this.getPointIndex (x, y, z)), Float.$valueOf (v));
}, "~N,~N,~N,~N");
Clazz_defineMethod (c$, "calculateFractionalPoint", 
function (cutoff, pointA, pointB, valueA, valueB, pt) {
var d = (valueB - valueA);
var fraction = (cutoff - valueA) / d;
this.edgeVector.sub2 (pointB, pointA);
pt.scaleAdd2 (fraction, this.edgeVector, pointA);
if (this.sr == null || !this.doIterate || valueB == valueA || fraction < 0.01 || fraction > 0.99 || (this.edgeVector.length ()) < 0.01) return cutoff;
var n = 0;
this.ptTemp.setT (pt);
var v = this.lookupInterpolatedVoxelValue (this.ptTemp, false);
var v0 = NaN;
while (++n < 10) {
var fnew = (v - valueA) / d;
if (fnew < 0 || fnew > 1) break;
var diff = (cutoff - v) / d / 2;
fraction += diff;
if (fraction < 0 || fraction > 1) break;
pt.setT (this.ptTemp);
v0 = v;
if (Math.abs (diff) < 0.005) break;
this.ptTemp.scaleAdd2 (diff, this.edgeVector, pt);
v = this.lookupInterpolatedVoxelValue (this.ptTemp, false);
}
return v0;
}, "~N,JU.P3,JU.P3,~N,~N,JU.P3");
});
Clazz_declarePackage ("J.jvxl.data");
Clazz_load (null, "J.jvxl.data.JvxlData", ["java.lang.Float", "JU.SB", "J.jvxl.data.JvxlCoder"], function () {
c$ = Clazz_decorateAsClass (function () {
this.wasJvxl = false;
this.wasCubic = false;
this.jvxlFileTitle = null;
this.jvxlFileMessage = null;
this.jvxlSurfaceData = null;
this.jvxlEdgeData = null;
this.jvxlColorData = null;
this.jvxlVolumeDataXml = null;
this.jvxlExcluded = null;
this.jvxlPlane = null;
this.isJvxlPrecisionColor = false;
this.jvxlDataIsColorMapped = false;
this.jvxlDataIs2dContour = false;
this.jvxlDataIsColorDensity = false;
this.isColorReversed = false;
this.thisSet = -2147483648;
this.edgeFractionBase = 35;
this.edgeFractionRange = 90;
this.colorFractionBase = 35;
this.colorFractionRange = 90;
this.dataXYReversed = false;
this.insideOut = false;
this.isXLowToHigh = false;
this.isContoured = false;
this.isBicolorMap = false;
this.isTruncated = false;
this.isCutoffAbsolute = false;
this.vertexDataOnly = false;
this.mappedDataMin = 0;
this.mappedDataMax = 0;
this.valueMappedToRed = 0;
this.valueMappedToBlue = 0;
this.cutoff = 0;
this.pointsPerAngstrom = 0;
this.nPointsX = 0;
this.nPointsY = 0;
this.nPointsZ = 0;
this.nBytes = 0;
this.nContours = 0;
this.nEdges = 0;
this.nSurfaceInts = 0;
this.vertexCount = 0;
this.vContours = null;
this.contourColixes = null;
this.contourColors = null;
this.contourValues = null;
this.contourValuesUsed = null;
this.scale3d = 0;
this.minColorIndex = -1;
this.maxColorIndex = 0;
this.title = null;
this.version = null;
this.boundingBox = null;
this.excludedTriangleCount = 0;
this.excludedVertexCount = 0;
this.colorDensity = false;
this.pointSize = 0;
this.moleculeXml = null;
this.dataMin = 0;
this.dataMax = 0;
this.saveVertexCount = 0;
this.vertexColorMap = null;
this.nVertexColors = 0;
this.vertexColors = null;
this.color = null;
this.meshColor = null;
this.translucency = 0;
this.colorScheme = null;
this.rendering = null;
this.slabValue = -2147483648;
this.isSlabbable = false;
this.diameter = 0;
this.slabInfo = null;
this.allowVolumeRender = false;
this.voxelVolume = 0;
this.mapLattice = null;
this.baseColor = null;
Clazz_instantialize (this, arguments);
}, J.jvxl.data, "JvxlData");
Clazz_prepareFields (c$, function () {
this.jvxlExcluded =  new Array (4);
});
Clazz_makeConstructor (c$, 
function () {
});
Clazz_defineMethod (c$, "clear", 
function () {
this.allowVolumeRender = true;
this.jvxlSurfaceData = "";
this.jvxlEdgeData = "";
this.jvxlColorData = "";
this.jvxlVolumeDataXml = "";
this.color = null;
this.colorScheme = null;
this.colorDensity = false;
this.pointSize = NaN;
this.contourValues = null;
this.contourValuesUsed = null;
this.contourColixes = null;
this.contourColors = null;
this.isSlabbable = false;
this.mapLattice = null;
this.meshColor = null;
this.nPointsX = 0;
this.nVertexColors = 0;
this.slabInfo = null;
this.slabValue = -2147483648;
this.thisSet = -2147483648;
this.rendering = null;
this.translucency = 0;
this.vContours = null;
this.vertexColorMap = null;
this.vertexColors = null;
this.voxelVolume = 0;
});
Clazz_defineMethod (c$, "setSurfaceInfo", 
function (thePlane, mapLattice, nSurfaceInts, surfaceData) {
this.jvxlSurfaceData = surfaceData;
if (this.jvxlSurfaceData.indexOf ("--") == 0) this.jvxlSurfaceData = this.jvxlSurfaceData.substring (2);
this.jvxlPlane = thePlane;
this.mapLattice = mapLattice;
this.nSurfaceInts = nSurfaceInts;
}, "JU.P4,JU.P3,~N,~S");
Clazz_defineMethod (c$, "setSurfaceInfoFromBitSet", 
function (bs, thePlane) {
this.setSurfaceInfoFromBitSetPts (bs, thePlane, null);
}, "JU.BS,JU.P4");
Clazz_defineMethod (c$, "setSurfaceInfoFromBitSetPts", 
function (bs, thePlane, mapLattice) {
var sb =  new JU.SB ();
var nSurfaceInts = (thePlane != null ? 0 : J.jvxl.data.JvxlCoder.jvxlEncodeBitSetBuffer (bs, this.nPointsX * this.nPointsY * this.nPointsZ, sb));
this.setSurfaceInfo (thePlane, mapLattice, nSurfaceInts, sb.toString ());
}, "JU.BS,JU.P4,JU.P3");
Clazz_defineMethod (c$, "jvxlUpdateInfo", 
function (title, nBytes) {
this.title = title;
this.nBytes = nBytes;
}, "~A,~N");
c$.updateSurfaceData = Clazz_defineMethod (c$, "updateSurfaceData", 
function (edgeData, vertexValues, vertexCount, vertexIncrement, isNaN) {
if (edgeData.length == 0) return "";
var chars = edgeData.toCharArray ();
for (var i = 0, ipt = 0; i < vertexCount; i += vertexIncrement, ipt++) if (Float.isNaN (vertexValues[i])) chars[ipt] = isNaN;

return String.copyValueOf (chars);
}, "~S,~A,~N,~N,~S");
});
Clazz_declarePackage ("J.jvxl.data");
Clazz_load (["JU.MeshSurface"], "J.jvxl.data.MeshData", ["java.lang.Float", "java.util.Arrays", "JU.AU", "$.BS", "$.V3"], function () {
c$ = Clazz_decorateAsClass (function () {
this.setsSuccessful = false;
this.vertexIncrement = 1;
this.polygonColorData = null;
if (!Clazz_isClassDefined ("J.jvxl.data.MeshData.SSet")) {
J.jvxl.data.MeshData.$MeshData$SSet$ ();
}
if (!Clazz_isClassDefined ("J.jvxl.data.MeshData.SortSet")) {
J.jvxl.data.MeshData.$MeshData$SortSet$ ();
}
Clazz_instantialize (this, arguments);
}, J.jvxl.data, "MeshData", JU.MeshSurface);
Clazz_defineMethod (c$, "addVertexCopy", 
function (vertex, value, assocVertex, asCopy) {
if (assocVertex < 0) this.vertexIncrement = -assocVertex;
return this.addVCVal (vertex, value, asCopy);
}, "JU.T3,~N,~N,~B");
Clazz_defineMethod (c$, "getSurfaceSet", 
function () {
return (this.surfaceSet == null ? this.getSurfaceSetForLevel (0) : this.surfaceSet);
});
Clazz_defineMethod (c$, "getSurfaceSetForLevel", 
 function (level) {
if (level == 0) {
this.surfaceSet =  new Array (100);
this.nSets = 0;
}this.setsSuccessful = true;
for (var i = 0; i < this.pc; i++) if (this.pis[i] != null) {
if (this.bsSlabDisplay != null && !this.bsSlabDisplay.get (i)) continue;
var p = this.pis[i];
var pt0 = this.findSet (p[0]);
var pt1 = this.findSet (p[1]);
var pt2 = this.findSet (p[2]);
if (pt0 < 0 && pt1 < 0 && pt2 < 0) {
this.createSet (p[0], p[1], p[2]);
continue;
}if (pt0 == pt1 && pt1 == pt2) continue;
if (pt0 >= 0) {
this.surfaceSet[pt0].set (p[1]);
this.surfaceSet[pt0].set (p[2]);
if (pt1 >= 0 && pt1 != pt0) this.mergeSets (pt0, pt1);
if (pt2 >= 0 && pt2 != pt0 && pt2 != pt1) this.mergeSets (pt0, pt2);
continue;
}if (pt1 >= 0) {
this.surfaceSet[pt1].set (p[0]);
this.surfaceSet[pt1].set (p[2]);
if (pt2 >= 0 && pt2 != pt1) this.mergeSets (pt1, pt2);
continue;
}this.surfaceSet[pt2].set (p[0]);
this.surfaceSet[pt2].set (p[1]);
}
var n = 0;
for (var i = 0; i < this.nSets; i++) if (this.surfaceSet[i] != null) n++;

var temp =  new Array (this.surfaceSet.length);
n = 0;
for (var i = 0; i < this.nSets; i++) if (this.surfaceSet[i] != null) temp[n++] = this.surfaceSet[i];

this.nSets = n;
this.surfaceSet = temp;
if (!this.setsSuccessful && level < 2) this.getSurfaceSetForLevel (level + 1);
if (level == 0) {
this.sortSurfaceSets ();
this.setVertexSets (false);
}return this.surfaceSet;
}, "~N");
Clazz_defineMethod (c$, "sortSurfaceSets", 
 function () {
var sets =  new Array (this.nSets);
for (var i = 0; i < this.nSets; i++) sets[i] = Clazz_innerTypeInstance (J.jvxl.data.MeshData.SSet, this, null, this.surfaceSet[i]);

java.util.Arrays.sort (sets, Clazz_innerTypeInstance (J.jvxl.data.MeshData.SortSet, this, null));
for (var i = 0; i < this.nSets; i++) this.surfaceSet[i] = sets[i].bs;

});
Clazz_defineMethod (c$, "setVertexSets", 
function (onlyIfNull) {
if (this.surfaceSet == null) return;
var nNull = 0;
for (var i = 0; i < this.nSets; i++) {
if (this.surfaceSet[i] != null && this.surfaceSet[i].cardinality () == 0) this.surfaceSet[i] = null;
if (this.surfaceSet[i] == null) nNull++;
}
if (nNull > 0) {
var bsNew =  new Array (this.nSets - nNull);
for (var i = 0, n = 0; i < this.nSets; i++) if (this.surfaceSet[i] != null) bsNew[n++] = this.surfaceSet[i];

this.surfaceSet = bsNew;
this.nSets -= nNull;
} else if (onlyIfNull) {
return;
}this.vertexSets =  Clazz_newIntArray (this.vc, 0);
for (var i = 0; i < this.nSets; i++) for (var j = this.surfaceSet[i].nextSetBit (0); j >= 0; j = this.surfaceSet[i].nextSetBit (j + 1)) this.vertexSets[j] = i;


}, "~B");
Clazz_defineMethod (c$, "findSet", 
 function (vertex) {
for (var i = 0; i < this.nSets; i++) if (this.surfaceSet[i] != null && this.surfaceSet[i].get (vertex)) return i;

return -1;
}, "~N");
Clazz_defineMethod (c$, "createSet", 
 function (v1, v2, v3) {
var i;
for (i = 0; i < this.nSets; i++) if (this.surfaceSet[i] == null) break;

if (i == this.surfaceSet.length) this.surfaceSet = JU.AU.ensureLength (this.surfaceSet, this.surfaceSet.length + 100);
this.surfaceSet[i] =  new JU.BS ();
this.surfaceSet[i].set (v1);
this.surfaceSet[i].set (v2);
this.surfaceSet[i].set (v3);
if (i == this.nSets) this.nSets++;
}, "~N,~N,~N");
Clazz_defineMethod (c$, "mergeSets", 
 function (a, b) {
this.surfaceSet[a].or (this.surfaceSet[b]);
this.surfaceSet[b] = null;
}, "~N,~N");
Clazz_defineMethod (c$, "invalidateSurfaceSet", 
function (i) {
for (var j = this.surfaceSet[i].nextSetBit (0); j >= 0; j = this.surfaceSet[i].nextSetBit (j + 1)) this.vvs[j] = NaN;

this.surfaceSet[i] = null;
}, "~N");
c$.checkCutoff = Clazz_defineMethod (c$, "checkCutoff", 
function (iA, iB, iC, vertexValues) {
if (iA < 0 || iB < 0 || iC < 0) return false;
var val1 = vertexValues[iA];
var val2 = vertexValues[iB];
var val3 = vertexValues[iC];
return (val1 >= 0 && val2 >= 0 && val3 >= 0 || val1 <= 0 && val2 <= 0 && val3 <= 0);
}, "~N,~N,~N,~A");
Clazz_defineMethod (c$, "calculateVolumeOrArea", 
function (thisSet, isArea, getSets) {
if (getSets || this.nSets == 0) this.getSurfaceSet ();
var justOne = (thisSet >= -1);
var n = (justOne || this.nSets == 0 ? 1 : this.nSets);
var v =  Clazz_newDoubleArray (n, 0);
var vAB =  new JU.V3 ();
var vAC =  new JU.V3 ();
var vTemp =  new JU.V3 ();
for (var i = this.pc; --i >= 0; ) {
if (!this.setABC (i)) continue;
var iSet = (this.nSets == 0 ? 0 : this.vertexSets[this.iA]);
if (thisSet >= 0 && iSet != thisSet) continue;
if (isArea) {
vAB.sub2 (this.vs[this.iB], this.vs[this.iA]);
vAC.sub2 (this.vs[this.iC], this.vs[this.iA]);
vTemp.cross (vAB, vAC);
v[justOne ? 0 : iSet] += vTemp.length ();
} else {
vAB.setT (this.vs[this.iB]);
vAC.setT (this.vs[this.iC]);
vTemp.cross (vAB, vAC);
vAC.setT (this.vs[this.iA]);
v[justOne ? 0 : iSet] += vAC.dot (vTemp);
}}
var factor = (isArea ? 2 : 6);
for (var i = 0; i < n; i++) v[i] /= factor;

if (justOne) return Float.$valueOf (v[0]);
return v;
}, "~N,~B,~B");
Clazz_defineMethod (c$, "updateInvalidatedVertices", 
function (bs) {
bs.clearAll ();
for (var i = 0; i < this.vc; i += this.vertexIncrement) if (Float.isNaN (this.vvs[i])) bs.set (i);

}, "JU.BS");
Clazz_defineMethod (c$, "invalidateVertices", 
function (bsInvalid) {
for (var i = bsInvalid.nextSetBit (0); i >= 0; i = bsInvalid.nextSetBit (i + 1)) this.vvs[i] = NaN;

}, "JU.BS");
c$.$MeshData$SSet$ = function () {
Clazz_pu$h(self.c$);
c$ = Clazz_decorateAsClass (function () {
Clazz_prepareCallback (this, arguments);
this.bs = null;
this.n = 0;
Clazz_instantialize (this, arguments);
}, J.jvxl.data.MeshData, "SSet");
Clazz_makeConstructor (c$, 
function (a) {
this.bs = a;
this.n = a.cardinality ();
}, "JU.BS");
c$ = Clazz_p0p ();
};
c$.$MeshData$SortSet$ = function () {
Clazz_pu$h(self.c$);
c$ = Clazz_decorateAsClass (function () {
Clazz_prepareCallback (this, arguments);
Clazz_instantialize (this, arguments);
}, J.jvxl.data.MeshData, "SortSet", null, java.util.Comparator);
Clazz_overrideMethod (c$, "compare", 
function (a, b) {
return (a.n > b.n ? -1 : a.n < b.n ? 1 : 0);
}, "J.jvxl.data.MeshData.SSet,J.jvxl.data.MeshData.SSet");
c$ = Clazz_p0p ();
};
Clazz_defineStatics (c$,
"MODE_GET_VERTICES", 1,
"MODE_GET_COLOR_INDEXES", 2,
"MODE_PUT_SETS", 3,
"MODE_PUT_VERTICES", 4);
});
Clazz_declarePackage ("J.jvxl.readers");
Clazz_load (null, "J.jvxl.readers.XmlReader", ["java.lang.Character", "JU.P3", "$.SB", "$.XmlUtil", "JU.Escape"], function () {
c$ = Clazz_decorateAsClass (function () {
this.br = null;
this.line = null;
Clazz_instantialize (this, arguments);
}, J.jvxl.readers, "XmlReader");
Clazz_defineMethod (c$, "getLine", 
function () {
return this.line;
});
Clazz_makeConstructor (c$, 
function (br) {
this.br = br;
}, "java.io.BufferedReader");
Clazz_defineMethod (c$, "toTag", 
function (name) {
this.skipTo ("<" + name);
if (this.line == null) return "";
var i = this.line.indexOf ("<" + name) + name.length + 1;
if (i == this.line.length) return this.line;
if (this.line.charAt (i) == ' ' || this.line.charAt (i) == '>') return this.line;
this.line = null;
return this.toTag (name);
}, "~S");
Clazz_defineMethod (c$, "skipTag", 
function (name) {
this.skipTo ("</" + name + ">");
}, "~S");
Clazz_defineMethod (c$, "getXmlData", 
function (name, data, withTag, allowSelfCloseOption) {
var closer = "</" + name + ">";
var tag = "<" + name;
if (data == null) {
var sb =  new JU.SB ();
try {
if (this.line == null) this.line = this.br.readLine ();
while (this.line.indexOf (tag) < 0) {
this.line = this.br.readLine ();
}
} catch (e) {
if (Clazz_exceptionOf (e, Exception)) {
return null;
} else {
throw e;
}
}
sb.append (this.line);
var selfClosed = false;
var pt = this.line.indexOf ("/>");
var pt1 = this.line.indexOf (">");
if (pt1 < 0 || pt == pt1 - 1) selfClosed = allowSelfCloseOption;
while (this.line.indexOf (closer) < 0 && (!selfClosed || this.line.indexOf ("/>") < 0)) sb.append (this.line = this.br.readLine ());

data = sb.toString ();
}return J.jvxl.readers.XmlReader.extractTag (data, tag, closer, withTag);
}, "~S,~S,~B,~B");
c$.extractTagOnly = Clazz_defineMethod (c$, "extractTagOnly", 
function (data, tag) {
return J.jvxl.readers.XmlReader.extractTag (data, "<" + tag + ">", "</" + tag + ">", false);
}, "~S,~S");
c$.extractTag = Clazz_defineMethod (c$, "extractTag", 
 function (data, tag, closer, withTag) {
var pt1 = data.indexOf (tag);
if (pt1 < 0) return "";
var pt2 = data.indexOf (closer, pt1);
if (pt2 < 0) {
pt2 = data.indexOf ("/>", pt1);
closer = "/>";
}if (pt2 < 0) return "";
if (withTag) {
pt2 += closer.length;
return data.substring (pt1, pt2);
}var quoted = false;
for (; pt1 < pt2; pt1++) {
var ch;
if ((ch = data.charAt (pt1)) == '"') quoted = !quoted;
 else if (quoted && ch == '\\') pt1++;
 else if (!quoted && (ch == '>' || ch == '/')) break;
}
if (pt1 >= pt2) return "";
while (Character.isWhitespace (data.charAt (++pt1))) {
}
return JU.XmlUtil.unwrapCdata (data.substring (pt1, pt2));
}, "~S,~S,~S,~B");
c$.getXmlAttrib = Clazz_defineMethod (c$, "getXmlAttrib", 
function (data, what) {
var nexta =  Clazz_newIntArray (1, 0);
var pt = J.jvxl.readers.XmlReader.setNext (data, what, nexta, 1);
if (pt < 2 || (pt = J.jvxl.readers.XmlReader.setNext (data, "\"", nexta, 0)) < 2) return "";
var pt1 = J.jvxl.readers.XmlReader.setNext (data, "\"", nexta, -1);
return (pt1 <= 0 ? "" : data.substring (pt, pt1));
}, "~S,~S");
Clazz_defineMethod (c$, "getXmlPoint", 
function (data, key) {
var spt = J.jvxl.readers.XmlReader.getXmlAttrib (data, key).$replace ('(', '{').$replace (')', '}');
var value = JU.Escape.uP (spt);
if (Clazz_instanceOf (value, JU.P3)) return value;
return  new JU.P3 ();
}, "~S,~S");
c$.setNext = Clazz_defineMethod (c$, "setNext", 
 function (data, what, next, offset) {
var ipt = next[0];
if (ipt < 0 || (ipt = data.indexOf (what, next[0])) < 0) return -1;
ipt += what.length;
next[0] = ipt + offset;
if (offset > 0 && ipt < data.length && data.charAt (ipt) != '=') return J.jvxl.readers.XmlReader.setNext (data, what, next, offset);
return next[0];
}, "~S,~S,~A,~N");
Clazz_defineMethod (c$, "skipTo", 
 function (key) {
if (this.line == null) this.line = this.br.readLine ();
while (this.line != null && this.line.indexOf (key) < 0) this.line = this.br.readLine ();

}, "~S");
Clazz_defineMethod (c$, "isNext", 
function (name) {
if (this.line == null || this.line.indexOf ("</") >= 0 && this.line.indexOf ("</") == this.line.indexOf ("<")) this.line = this.br.readLine ();
return (this.line.indexOf ("<" + name) >= 0);
}, "~S");
});
Clazz_declarePackage ("J.jvxl.readers");
Clazz_load (["JU.P3", "$.V3"], "J.jvxl.readers.SurfaceGenerator", ["java.lang.Float", "java.util.Map", "JU.AU", "$.BS", "$.P4", "$.PT", "$.Rdr", "J.io.JmolBinary", "J.jvxl.data.JvxlCoder", "$.JvxlData", "$.MeshData", "$.VolumeData", "J.jvxl.readers.Parameters", "$.SurfaceReader", "JU.Logger", "$.Measure"], function () {
c$ = Clazz_decorateAsClass (function () {
this.jvxlData = null;
this.meshData = null;
this.params = null;
this.volumeData = null;
this.meshDataServer = null;
this.atomDataServer = null;
this.marchingSquares = null;
this.version = null;
this.$isValid = true;
this.fileType = null;
this.out = null;
this.surfaceReader = null;
this.colorPtr = 0;
this.readerData = null;
this.vAC = null;
this.vAB = null;
this.vNorm = null;
this.ptRef = null;
this.bsVdw = null;
Clazz_instantialize (this, arguments);
}, J.jvxl.readers, "SurfaceGenerator");
Clazz_prepareFields (c$, function () {
this.vAC =  new JU.V3 ();
this.vAB =  new JU.V3 ();
this.vNorm =  new JU.V3 ();
this.ptRef = JU.P3.new3 (0, 0, 1e15);
});
Clazz_defineMethod (c$, "isValid", 
function () {
return this.$isValid;
});
Clazz_defineMethod (c$, "getFileType", 
function () {
return this.fileType;
});
Clazz_defineMethod (c$, "setVersion", 
function (version) {
this.version = version;
}, "~S");
Clazz_makeConstructor (c$, 
function (atomDataServer, meshDataServer, meshData, jvxlData) {
this.atomDataServer = atomDataServer;
this.meshDataServer = meshDataServer;
this.params =  new J.jvxl.readers.Parameters ();
this.meshData = (meshData == null ?  new J.jvxl.data.MeshData () : meshData);
this.jvxlData = (jvxlData == null ?  new J.jvxl.data.JvxlData () : jvxlData);
this.volumeData =  new J.jvxl.data.VolumeData ();
this.initializeIsosurface ();
}, "J.atomdata.AtomDataServer,J.jvxl.api.MeshDataServer,J.jvxl.data.MeshData,J.jvxl.data.JvxlData");
Clazz_defineMethod (c$, "isStateDataRead", 
function () {
return this.params.state == 2;
});
Clazz_defineMethod (c$, "getFileName", 
function () {
return this.params.fileName;
});
Clazz_defineMethod (c$, "getMeshDataServer", 
function () {
return this.meshDataServer;
});
Clazz_defineMethod (c$, "getAtomDataServer", 
function () {
return this.atomDataServer;
});
Clazz_defineMethod (c$, "getColorEncoder", 
function () {
return this.params.colorEncoder;
});
Clazz_defineMethod (c$, "getVertexSource", 
function () {
return this.params.vertexSource;
});
Clazz_defineMethod (c$, "setJvxlData", 
function (jvxlData) {
this.jvxlData = jvxlData;
if (jvxlData != null) jvxlData.version = this.version;
}, "J.jvxl.data.JvxlData");
Clazz_defineMethod (c$, "getJvxlData", 
function () {
return this.jvxlData;
});
Clazz_defineMethod (c$, "getMeshData", 
function () {
return this.meshData;
});
Clazz_defineMethod (c$, "setMarchingSquares", 
function (marchingSquares) {
this.marchingSquares = marchingSquares;
}, "J.jvxl.calc.MarchingSquares");
Clazz_defineMethod (c$, "getMarchingSquares", 
function () {
return this.marchingSquares;
});
Clazz_defineMethod (c$, "getParams", 
function () {
return this.params;
});
Clazz_defineMethod (c$, "getScript", 
function () {
return this.params.script;
});
Clazz_defineMethod (c$, "getTitle", 
function () {
return this.params.title;
});
Clazz_defineMethod (c$, "getBsSelected", 
function () {
return this.params.bsSelected;
});
Clazz_defineMethod (c$, "getBsIgnore", 
function () {
return this.params.bsIgnore;
});
Clazz_defineMethod (c$, "getVolumeData", 
function () {
return this.volumeData;
});
Clazz_defineMethod (c$, "getPlane", 
function () {
return this.params.thePlane;
});
Clazz_defineMethod (c$, "getColor", 
function (which) {
switch (which) {
case -1:
return this.params.colorNeg;
case 1:
return this.params.colorPos;
}
return 0;
}, "~N");
Clazz_defineMethod (c$, "setModelIndex", 
function (modelIndex) {
this.params.modelIndex = modelIndex;
}, "~N");
Clazz_defineMethod (c$, "getIAddGridPoints", 
function () {
return this.params.iAddGridPoints;
});
Clazz_defineMethod (c$, "getIsPositiveOnly", 
function () {
return this.params.isPositiveOnly;
});
Clazz_defineMethod (c$, "isInsideOut", 
function () {
return this.params.insideOut != this.params.dataXYReversed;
});
Clazz_defineMethod (c$, "getCutoff", 
function () {
return this.params.cutoff;
});
Clazz_defineMethod (c$, "getMoData", 
function () {
return this.params.moData;
});
Clazz_defineMethod (c$, "isCubeData", 
function () {
return this.jvxlData.wasCubic;
});
Clazz_defineMethod (c$, "setParameter", 
function (propertyName, value) {
return this.setProp (propertyName, value, null);
}, "~S,~O");
Clazz_defineMethod (c$, "setProp", 
function (propertyName, value, bs) {
if ("debug" === propertyName) {
var TF = (value).booleanValue ();
this.params.logMessages = TF;
this.params.logCube = TF;
return true;
}if ("init" === propertyName) {
this.initializeIsosurface ();
if (Clazz_instanceOf (value, J.jvxl.readers.Parameters)) {
this.params = value;
} else {
this.params.script = value;
if (this.params.script != null && this.params.script.indexOf (";#") >= 0) {
this.params.script = JU.PT.rep (this.params.script, ";#", "; #");
}}return false;
}if ("map" === propertyName) {
this.params.resetForMapping ((value).booleanValue ());
if (this.surfaceReader != null) this.surfaceReader.minMax = null;
return true;
}if ("finalize" === propertyName) {
this.initializeIsosurface ();
return true;
}if ("clear" === propertyName) {
if (this.surfaceReader != null) this.surfaceReader.discardTempData (true);
return false;
}if ("fileIndex" === propertyName) {
this.params.fileIndex = (value).intValue ();
if (this.params.fileIndex < 0) this.params.fileIndex = 0;
this.params.readAllData = false;
return true;
}if ("blockData" === propertyName) {
this.params.blockCubeData = (value).booleanValue ();
return true;
}if ("withinPoints" === propertyName) {
this.params.boundingBox = (value)[1];
return true;
}if ("boundingBox" === propertyName) {
var pts = value;
this.params.boundingBox = [JU.P3.newP (pts[0]), JU.P3.newP (pts[pts.length - 1])];
return true;
}if ("func" === propertyName) {
this.params.func = value;
return true;
}if ("intersection" === propertyName) {
this.params.intersection = value;
return true;
}if ("bsSolvent" === propertyName) {
this.params.bsSolvent = value;
return true;
}if ("select" === propertyName) {
this.params.bsSelected = value;
return true;
}if ("ignore" === propertyName) {
this.params.bsIgnore = value;
return true;
}if ("propertySmoothing" === propertyName) {
this.params.propertySmoothing = (value).booleanValue ();
return true;
}if ("propertyDistanceMax" === propertyName) {
this.params.propertyDistanceMax = (value).floatValue ();
return true;
}if ("propertySmoothingPower" === propertyName) {
this.params.propertySmoothingPower = (value).intValue ();
return true;
}if ("title" === propertyName) {
if (value == null) {
this.params.title = null;
return true;
} else if (JU.PT.isAS (value)) {
this.params.title = value;
for (var i = 0; i < this.params.title.length; i++) if (this.params.title[i].length > 0) JU.Logger.info (this.params.title[i]);

}return true;
}if ("sigma" === propertyName) {
this.params.cutoff = this.params.sigma = (value).floatValue ();
this.params.isPositiveOnly = false;
this.params.cutoffAutomatic = false;
return true;
}if ("cutoff" === propertyName) {
this.params.cutoff = (value).floatValue ();
this.params.isPositiveOnly = false;
this.params.cutoffAutomatic = false;
return true;
}if ("parameters" === propertyName) {
this.params.parameters = JU.AU.ensureLengthA (value, 2);
if (this.params.parameters.length > 0 && this.params.parameters[0] != 0) this.params.cutoff = this.params.parameters[0];
return true;
}if ("cutoffPositive" === propertyName) {
this.params.cutoff = (value).floatValue ();
this.params.isPositiveOnly = true;
return true;
}if ("cap" === propertyName || "slab" === propertyName) {
if (value != null) this.params.addSlabInfo (value);
return true;
}if ("scale" === propertyName) {
this.params.scale = (value).floatValue ();
return true;
}if ("scale3d" === propertyName) {
this.params.scale3d = (value).floatValue ();
return true;
}if ("angstroms" === propertyName) {
this.params.isAngstroms = true;
return true;
}if ("resolution" === propertyName) {
var resolution = (value).floatValue ();
this.params.resolution = (resolution > 0 ? resolution : 3.4028235E38);
return true;
}if ("downsample" === propertyName) {
var rate = (value).intValue ();
this.params.downsampleFactor = (rate >= 0 ? rate : 0);
return true;
}if ("anisotropy" === propertyName) {
if ((this.params.dataType & 32) == 0) this.params.setAnisotropy (value);
return true;
}if ("eccentricity" === propertyName) {
this.params.setEccentricity (value);
return true;
}if ("addHydrogens" === propertyName) {
this.params.addHydrogens = (value).booleanValue ();
return true;
}if ("squareData" === propertyName) {
this.params.isSquared = (value == null ? false : (value).booleanValue ());
return true;
}if ("squareLinear" === propertyName) {
this.params.isSquaredLinear = (value == null ? false : (value).booleanValue ());
return true;
}if ("gridPoints" === propertyName) {
this.params.iAddGridPoints = true;
return true;
}if ("atomIndex" === propertyName) {
this.params.atomIndex = (value).intValue ();
return true;
}if ("insideOut" === propertyName) {
this.params.insideOut = true;
return true;
}if ("sign" === propertyName) {
this.params.isCutoffAbsolute = true;
this.params.colorBySign = true;
this.colorPtr = 0;
return true;
}if ("colorRGB" === propertyName) {
var rgb = (value).intValue ();
this.params.colorRgb = this.params.colorPos = this.params.colorPosLCAO = rgb;
if (this.colorPtr++ == 0) {
this.params.colorNeg = this.params.colorNegLCAO = rgb;
} else {
this.params.colorRgb = 2147483647;
}return true;
}if ("monteCarloCount" === propertyName) {
this.params.psi_monteCarloCount = (value).intValue ();
return true;
}if ("rangeAll" === propertyName) {
this.params.rangeAll = true;
return true;
}if ("rangeSelected" === propertyName) {
this.params.rangeSelected = true;
return true;
}if ("red" === propertyName) {
this.params.valueMappedToRed = (value).floatValue ();
return true;
}if ("blue" === propertyName) {
this.params.valueMappedToBlue = (value).floatValue ();
if (this.params.valueMappedToRed > this.params.valueMappedToBlue) {
var f = this.params.valueMappedToRed;
this.params.valueMappedToRed = this.params.valueMappedToBlue;
this.params.valueMappedToBlue = f;
this.params.isColorReversed = !this.params.isColorReversed;
}this.params.rangeDefined = true;
this.params.rangeAll = false;
return true;
}if ("reverseColor" === propertyName) {
this.params.isColorReversed = true;
return true;
}if ("setColorScheme" === propertyName) {
this.getSurfaceSets ();
this.params.colorBySets = true;
this.mapSurface ();
return true;
}if ("center" === propertyName) {
this.params.center.setT (value);
return true;
}if ("volumeData" === propertyName) {
this.params.volumeData = value;
return true;
}if ("origin" === propertyName) {
this.params.origin = value;
return true;
}if ("step" === propertyName) {
this.params.steps = value;
return true;
}if ("point" === propertyName) {
this.params.points = value;
return true;
}if ("withinDistance" === propertyName) {
this.params.distance = (value).floatValue ();
return true;
}if ("withinPoint" === propertyName) {
this.params.point = value;
return true;
}if ("progressive" === propertyName) {
this.params.isXLowToHigh = true;
return true;
}if ("phase" === propertyName) {
var color = value;
this.params.isCutoffAbsolute = true;
this.params.colorBySign = true;
this.params.colorByPhase = true;
this.params.colorPhase = J.jvxl.readers.SurfaceReader.getColorPhaseIndex (color);
if (this.params.colorPhase < 0) {
JU.Logger.warn (" invalid color phase: " + color);
this.params.colorPhase = 0;
}this.params.colorByPhase = this.params.colorPhase != 0;
if (this.params.state >= 2) {
this.params.dataType = this.params.surfaceType;
this.params.state = 3;
this.params.isBicolorMap = true;
this.surfaceReader.applyColorScale ();
}return true;
}if ("radius" === propertyName) {
JU.Logger.info ("solvent probe radius set to " + value);
this.params.atomRadiusData = value;
return true;
}if ("envelopeRadius" === propertyName) {
this.params.envelopeRadius = (value).floatValue ();
return true;
}if ("cavityRadius" === propertyName) {
this.params.cavityRadius = (value).floatValue ();
return true;
}if ("cavity" === propertyName) {
this.params.isCavity = true;
return true;
}if ("doFullMolecular" === propertyName) {
this.params.doFullMolecular = true;
return true;
}if ("pocket" === propertyName) {
this.params.pocket = value;
this.params.fullyLit = this.params.pocket.booleanValue ();
return true;
}if ("minset" === propertyName) {
this.params.minSet = (value).intValue ();
return true;
}if ("maxset" === propertyName) {
this.params.maxSet = (value).intValue ();
return true;
}if ("plane" === propertyName) {
this.params.setPlane (value);
return true;
}if ("contour" === propertyName) {
this.params.isContoured = true;
var n;
if (JU.PT.isAF (value)) {
this.params.contoursDiscrete = value;
this.params.nContours = this.params.contoursDiscrete.length;
} else if (Clazz_instanceOf (value, JU.P3)) {
var pt = this.params.contourIncrements = value;
var from = pt.x;
var to = pt.y;
var step = pt.z;
if (step <= 0) step = 1;
n = 0;
for (var p = from; p <= to + step / 10; p += step, n++) {
}
this.params.contoursDiscrete =  Clazz_newFloatArray (n, 0);
var p = from;
for (var i = 0; i < n; i++, p += step) {
this.params.contoursDiscrete[i] = p;
}
this.params.nContours = n;
} else {
n = (value).intValue ();
if (n == 0) this.params.nContours = 9;
 else if (n > 0) this.params.nContours = n;
 else this.params.thisContour = -n;
}return true;
}if ("colorDiscrete" === propertyName) {
this.params.contourColixes = value;
return true;
}if ("colorDensity" === propertyName) {
this.params.colorDensity = true;
if (value != null) this.params.pointSize = (value).floatValue ();
return false;
}if ("fullPlane" === propertyName) {
this.params.contourFromZero = !(value).booleanValue ();
return true;
}if ("mapLattice" === propertyName) {
this.params.mapLattice = value;
return true;
}if ("property" === propertyName) {
this.params.dataType = 1206;
this.params.theProperty = value;
this.mapSurface ();
return true;
}if ("sphere" === propertyName) {
this.params.setSphere ((value).floatValue (), false);
this.readerData = Float.$valueOf (this.params.distance);
this.surfaceReader = this.newReader ("IsoShapeReader");
this.generateSurface ();
return true;
}if ("geodesic" === propertyName) {
this.params.setSphere ((value).floatValue (), true);
this.readerData = Float.$valueOf (this.params.distance);
this.surfaceReader = this.newReader ("IsoShapeReader");
this.generateSurface ();
return true;
}if ("ellipsoid" === propertyName) {
if (Clazz_instanceOf (value, JU.P4)) this.params.setEllipsoidP4 (value);
 else if (JU.PT.isAF (value)) this.params.setEllipsoidAF (value);
 else return true;
this.readerData = Float.$valueOf (this.params.distance);
this.surfaceReader = this.newReader ("IsoShapeReader");
this.generateSurface ();
return true;
}if ("ellipsoid3" === propertyName) {
this.params.setEllipsoidAF (value);
this.readerData = Float.$valueOf (this.params.distance);
this.surfaceReader = this.newReader ("IsoShapeReader");
this.generateSurface ();
return true;
}if ("lp" === propertyName) {
this.params.setLp (value);
this.readerData = [3, 2, 0, 15, 0];
this.surfaceReader = this.newReader ("IsoShapeReader");
this.generateSurface ();
return true;
}if ("rad" === propertyName) {
this.params.setRadical (value);
this.readerData = [3, 2, 0, 15, 0];
this.surfaceReader = this.newReader ("IsoShapeReader");
this.generateSurface ();
return true;
}if ("lobe" === propertyName) {
this.params.setLobe (value);
this.readerData = [3, 2, 0, 15, 0];
this.surfaceReader = this.newReader ("IsoShapeReader");
this.generateSurface ();
return true;
}if ("hydrogenOrbital" === propertyName) {
if (!this.params.setAtomicOrbital (value)) {
this.$isValid = false;
return true;
}this.readerData = [this.params.psi_n, this.params.psi_l, this.params.psi_m, this.params.psi_Znuc, this.params.psi_monteCarloCount];
this.surfaceReader = this.newReader ("IsoShapeReader");
this.processState ();
return true;
}if ("functionXY" === propertyName) {
this.params.setFunctionXY (value);
if (this.params.isContoured) this.volumeData.setPlaneParameters (JU.P4.new4 (0, 0, 1, 0));
if ((this.params.functionInfo.get (0)).indexOf ("_xyz") >= 0) this.getFunctionZfromXY ();
this.processState ();
return true;
}if ("functionXYZ" === propertyName) {
this.params.setFunctionXYZ (value);
this.processState ();
return true;
}if ("lcaoType" === propertyName) {
this.params.setLcao (value, this.colorPtr);
return true;
}if ("lcaoCartoonCenter" === propertyName) {
if (++this.params.state != 2) return true;
if (this.params.center.x == 3.4028235E38) this.params.center.setT (value);
return false;
}if ("molecular" === propertyName || "solvent" === propertyName || "sasurface" === propertyName || "nomap" === propertyName) {
this.params.setSolvent (propertyName, (value).floatValue ());
JU.Logger.info (this.params.calculationType);
this.processState ();
return true;
}if ("moData" === propertyName) {
this.params.moData = value;
return true;
}if ("mepCalcType" === propertyName) {
this.params.mep_calcType = (value).intValue ();
return true;
}if ("mep" === propertyName) {
this.params.setMep (value, false);
this.processState ();
return true;
}if ("mlp" === propertyName) {
this.params.setMep (value, true);
this.processState ();
return true;
}if ("nci" === propertyName) {
var isPromolecular = (value).booleanValue ();
this.params.setNci (isPromolecular);
if (isPromolecular) this.processState ();
return true;
}if ("calculationType" === propertyName) {
this.params.calculationType = value;
return true;
}if ("charges" === propertyName) {
this.params.theProperty = value;
return true;
}if ("randomSeed" === propertyName) {
this.params.randomSeed = (value).intValue ();
return true;
}if ("molecularOrbital" === propertyName) {
var iMo = 0;
var linearCombination = null;
if (Clazz_instanceOf (value, Integer)) {
iMo = (value).intValue ();
} else {
linearCombination = value;
}this.params.setMO (iMo, linearCombination);
JU.Logger.info (this.params.calculationType);
this.processState ();
return true;
}if ("fileType" === propertyName) {
this.fileType = value;
return true;
}if ("fileName" === propertyName) {
this.params.fileName = value;
return true;
}if ("outputChannel" === propertyName) {
this.out = value;
return true;
}if ("readFile" === propertyName) {
if ((this.surfaceReader = this.setFileData (value)) == null) {
JU.Logger.error ("Could not set the surface data");
return true;
}this.surfaceReader.setOutputChannel (this.out);
this.generateSurface ();
return true;
}if ("getSurfaceSets" === propertyName) {
this.getSurfaceSets ();
return true;
}if ("mapColor" === propertyName) {
if ((this.surfaceReader = this.setFileData (value)) == null) {
JU.Logger.error ("Could not set the mapping data");
return true;
}this.surfaceReader.setOutputChannel (this.out);
this.mapSurface ();
return true;
}if ("periodic" === propertyName) {
this.params.isPeriodic = true;
}return false;
}, "~S,~O,JU.BS");
Clazz_defineMethod (c$, "newReader", 
 function (name) {
var sr = J.jvxl.readers.SurfaceGenerator.getInterface (name);
if (sr != null) sr.init (this);
return sr;
}, "~S");
Clazz_defineMethod (c$, "newReaderBr", 
 function (name, br) {
var sr = J.jvxl.readers.SurfaceGenerator.getInterface (name);
if (sr != null) sr.init2 (this, br);
return sr;
}, "~S,java.io.BufferedReader");
c$.getInterface = Clazz_defineMethod (c$, "getInterface", 
 function (name) {
try {
var x = Clazz._4Name ("J.jvxl.readers." + name);
return (x == null ? null : x.newInstance ());
} catch (e) {
if (Clazz_exceptionOf (e, Exception)) {
JU.Logger.error ("Interface.java Error creating instance for " + name + ": \n" + e.toString ());
return null;
} else {
throw e;
}
}
}, "~S");
Clazz_defineMethod (c$, "getSurfaceSets", 
 function () {
if (this.meshDataServer == null) {
this.meshData.getSurfaceSet ();
} else {
this.meshDataServer.fillMeshData (this.meshData, 1, null);
this.meshData.getSurfaceSet ();
this.meshDataServer.fillMeshData (this.meshData, 3, null);
}});
Clazz_defineMethod (c$, "processState", 
 function () {
if (this.params.state == 1 && this.params.thePlane != null) this.params.state++;
if (this.params.state >= 2) {
this.mapSurface ();
} else {
this.generateSurface ();
}});
Clazz_defineMethod (c$, "setReader", 
 function () {
this.readerData = null;
if (this.surfaceReader != null) return !this.surfaceReader.vertexDataOnly;
switch (this.params.dataType) {
case 1205:
this.surfaceReader = this.newReader ("IsoPlaneReader");
break;
case 1206:
this.surfaceReader = this.newReader ("AtomPropertyMapper");
break;
case 1328:
case 1329:
this.readerData = (this.params.dataType == 1328 ? "Mep" : "Mlp");
if (this.params.state == 3) {
this.surfaceReader = this.newReader ("AtomPropertyMapper");
} else {
this.surfaceReader = this.newReader ("Iso" + this.readerData + "Reader");
}break;
case 1333:
this.surfaceReader = this.newReader ("IsoIntersectReader");
break;
case 1195:
case 1203:
case 1196:
this.surfaceReader = this.newReader ("IsoSolventReader");
break;
case 1844:
case 1837:
this.surfaceReader = this.newReader ("IsoMOReader");
break;
case 8:
this.surfaceReader = this.newReader ("IsoFxyReader");
break;
case 9:
this.surfaceReader = this.newReader ("IsoFxyzReader");
break;
}
JU.Logger.info ("Using surface reader " + this.surfaceReader);
return true;
});
Clazz_defineMethod (c$, "generateSurface", 
 function () {
if (++this.params.state != 2) return;
this.setReader ();
var haveMeshDataServer = (this.meshDataServer != null);
if (this.params.colorBySign) this.params.isBicolorMap = true;
if (this.surfaceReader == null) {
JU.Logger.error ("surfaceReader is null for " + this.params.dataType);
return;
}if (!this.surfaceReader.createIsosurface (false)) {
JU.Logger.error ("Could not create isosurface");
this.params.cutoff = NaN;
this.surfaceReader.closeReader ();
return;
}if (this.params.pocket != null && haveMeshDataServer) this.surfaceReader.selectPocket (!this.params.pocket.booleanValue ());
if (this.params.minSet > 0) this.surfaceReader.excludeMinimumSet ();
if (this.params.maxSet > 0) this.surfaceReader.excludeMaximumSet ();
if (this.params.slabInfo != null) this.surfaceReader.slabIsosurface (this.params.slabInfo);
if (haveMeshDataServer) this.meshDataServer.notifySurfaceGenerationCompleted ();
if (this.jvxlData.thisSet >= 0) this.getSurfaceSets ();
if (this.jvxlData.jvxlDataIs2dContour) {
this.surfaceReader.colorIsosurface ();
this.params.state = 3;
}if (this.jvxlData.jvxlDataIsColorDensity) {
this.params.state = 3;
}if (this.params.colorBySign || this.params.isBicolorMap) {
this.params.state = 3;
this.surfaceReader.applyColorScale ();
}if (this.jvxlData.vertexColorMap != null) {
this.jvxlData.vertexColorMap = null;
this.surfaceReader.hasColorData = false;
}this.surfaceReader.jvxlUpdateInfo ();
this.setMarchingSquares (this.surfaceReader.marchingSquares);
this.surfaceReader.discardTempData (false);
this.params.mappedDataMin = 3.4028235E38;
this.surfaceReader.closeReader ();
if (this.params.state != 3 && (this.surfaceReader.hasColorData || this.params.colorDensity)) {
this.params.state = 3;
this.colorIsosurface ();
}this.surfaceReader = null;
});
Clazz_defineMethod (c$, "mapSurface", 
 function () {
if (this.params.state == 1 && this.params.thePlane != null) this.params.state++;
if (++this.params.state < 3) return;
if (!this.setReader ()) return;
if (this.params.isPeriodic) this.volumeData.isPeriodic = true;
if (this.params.thePlane != null) {
var isSquared = this.params.isSquared;
this.params.isSquared = false;
this.params.cutoff = 0;
this.volumeData.setMappingPlane (this.params.thePlane);
this.surfaceReader.createIsosurface (!this.params.isPeriodic);
this.volumeData.setMappingPlane (null);
if (this.meshDataServer != null) this.meshDataServer.notifySurfaceGenerationCompleted ();
if (this.params.dataType == 1205) {
this.surfaceReader.discardTempData (true);
return;
}this.params.isSquared = isSquared;
this.params.mappedDataMin = 3.4028235E38;
this.surfaceReader.readVolumeData (true);
if (this.params.mapLattice != null) this.volumeData.isPeriodic = true;
} else if (!this.params.colorBySets && !this.params.colorDensity) {
this.surfaceReader.readAndSetVolumeParameters (true);
this.params.mappedDataMin = 3.4028235E38;
this.surfaceReader.readVolumeData (true);
}this.colorIsosurface ();
this.surfaceReader.closeReader ();
this.surfaceReader = null;
});
Clazz_defineMethod (c$, "getSlabInfo", 
function () {
return this.params.slabInfo;
});
Clazz_defineMethod (c$, "colorIsosurface", 
function () {
this.surfaceReader.colorIsosurface ();
this.surfaceReader.jvxlUpdateInfo ();
this.surfaceReader.updateTriangles ();
this.surfaceReader.discardTempData (true);
if (this.meshDataServer != null) this.meshDataServer.notifySurfaceMappingCompleted ();
});
Clazz_defineMethod (c$, "getProperty", 
function (property, index) {
if (property === "jvxlFileData") return J.jvxl.data.JvxlCoder.jvxlGetFile (this.jvxlData, null, this.params.title, "", true, index, null, null);
if (property === "jvxlFileInfo") return J.jvxl.data.JvxlCoder.jvxlGetInfo (this.jvxlData);
return null;
}, "~S,~N");
Clazz_defineMethod (c$, "setFileData", 
 function (value) {
var fileType = this.fileType;
this.fileType = null;
if (Clazz_instanceOf (value, J.jvxl.data.VolumeData)) {
this.volumeData = value;
return this.newReader ("VolumeDataReader");
}if (Clazz_instanceOf (value, java.util.Map)) {
var map = value;
if (map.containsKey ("__pymolSurfaceData__")) {
this.readerData = map;
return this.newReaderBr ("PyMOLMeshReader", null);
}this.volumeData = map.get ("volumeData");
return this.newReader ("VolumeDataReader");
}var data = null;
if (Clazz_instanceOf (value, String)) {
data = value;
value = JU.Rdr.getBR (value);
}var br = value;
if (fileType == null) fileType = J.io.JmolBinary.determineSurfaceFileType (br);
if (fileType != null && fileType.startsWith ("UPPSALA")) {
var fname = this.params.fileName;
fname = fname.substring (0, fname.indexOf ("/", 10));
fname += JU.PT.getQuotedStringAt (fileType, fileType.indexOf ("A HREF") + 1);
this.params.fileName = fname;
value = this.atomDataServer.getBufferedInputStream (fname);
if (value == null) {
JU.Logger.error ("Isosurface: could not open file " + fname);
return null;
}try {
br = JU.Rdr.getBufferedReader (value, null);
} catch (e) {
if (Clazz_exceptionOf (e, Exception)) {
} else {
throw e;
}
}
fileType = J.io.JmolBinary.determineSurfaceFileType (br);
}if (fileType == null) fileType = "UNKNOWN";
JU.Logger.info ("data file type was determined to be " + fileType);
if (fileType.equals ("Jvxl+")) return this.newReaderBr ("JvxlReader", br);
this.readerData = [this.params.fileName, data];
if ("MRC DELPHI DSN6".indexOf (fileType.toUpperCase ()) >= 0) {
try {
br.close ();
} catch (e) {
if (Clazz_exceptionOf (e, java.io.IOException)) {
} else {
throw e;
}
}
br = null;
fileType += "Binary";
}return this.newReaderBr (fileType + "Reader", br);
}, "~O");
Clazz_defineMethod (c$, "getReaderData", 
function () {
var o = this.readerData;
this.readerData = null;
return o;
});
Clazz_defineMethod (c$, "initializeIsosurface", 
function () {
this.params.initialize ();
this.colorPtr = 0;
this.surfaceReader = null;
this.marchingSquares = null;
this.initState ();
});
Clazz_defineMethod (c$, "initState", 
function () {
this.params.state = 1;
this.params.dataType = this.params.surfaceType = 0;
});
Clazz_defineMethod (c$, "setLcao", 
function () {
this.params.colorPos = this.params.colorPosLCAO;
this.params.colorNeg = this.params.colorNegLCAO;
return this.params.lcaoType;
});
Clazz_defineMethod (c$, "getFunctionZfromXY", 
 function () {
var origin = this.params.functionInfo.get (1);
var counts =  Clazz_newIntArray (3, 0);
var nearest =  Clazz_newIntArray (3, 0);
var vectors =  new Array (3);
for (var i = 0; i < 3; i++) {
var info = this.params.functionInfo.get (i + 2);
counts[i] = Math.abs (Clazz_floatToInt (info.x));
vectors[i] = JU.V3.new3 (info.y, info.z, info.w);
}
var nx = counts[0];
var ny = counts[1];
var pt =  new JU.P3 ();
var pta =  new JU.P3 ();
var ptb =  new JU.P3 ();
var ptc =  new JU.P3 ();
var data = this.params.functionInfo.get (5);
var data2 =  Clazz_newFloatArray (nx, ny, 0);
var d;
for (var i = 0; i < nx; i++) for (var j = 0; j < ny; j++) {
pt.scaleAdd2 (i, vectors[0], origin);
pt.scaleAdd2 (j, vectors[1], pt);
var dist = J.jvxl.readers.SurfaceGenerator.findNearestThreePoints (pt.x, pt.y, data, nearest);
pta.set ((d = data[nearest[0]])[0], d[1], d[2]);
if (dist < 0.00001) {
pt.z = d[2];
} else {
ptb.set ((d = data[nearest[1]])[0], d[1], d[2]);
ptc.set ((d = data[nearest[2]])[0], d[1], d[2]);
pt.z = this.distanceVerticalToPlane (pt.x, pt.y, pta, ptb, ptc);
}data2[i][j] = pt.z;
}

this.params.functionInfo.set (5, data2);
});
Clazz_defineMethod (c$, "distanceVerticalToPlane", 
 function (x, y, pta, ptb, ptc) {
var d = JU.Measure.getDirectedNormalThroughPoints (pta, ptb, ptc, this.ptRef, this.vNorm, this.vAB, this.vAC);
return (this.vNorm.x * x + this.vNorm.y * y + d) / -this.vNorm.z;
}, "~N,~N,JU.P3,JU.P3,JU.P3");
c$.findNearestThreePoints = Clazz_defineMethod (c$, "findNearestThreePoints", 
 function (x, y, xyz, result) {
var d;
var dist1;
var dist2;
var dist3;
var i1;
var i2;
var i3;
i1 = i2 = i3 = -1;
dist1 = dist2 = dist3 = 3.4028235E38;
for (var i = xyz.length; --i >= 0; ) {
d = (d = xyz[i][0] - x) * d + (d = xyz[i][1] - y) * d;
if (d < dist1) {
dist3 = dist2;
dist2 = dist1;
dist1 = d;
i3 = i2;
i2 = i1;
i1 = i;
} else if (d < dist2) {
dist3 = dist2;
dist2 = d;
i3 = i2;
i2 = i;
} else if (d < dist3) {
dist3 = d;
i3 = i;
}}
result[0] = i1;
result[1] = i2;
result[2] = i3;
return dist1;
}, "~N,~N,~A,~A");
Clazz_defineMethod (c$, "addRequiredFile", 
function (fileName) {
if (this.meshDataServer == null) return;
this.meshDataServer.addRequiredFile (fileName);
}, "~S");
Clazz_defineMethod (c$, "log", 
function (msg) {
if (this.atomDataServer == null) System.out.println (msg);
 else this.atomDataServer.log (msg);
}, "~S");
Clazz_defineMethod (c$, "setOutputChannel", 
function (binaryDoc, out) {
if (this.meshDataServer == null) return;
this.meshDataServer.setOutputChannel (binaryDoc, out);
}, "javajs.api.GenericBinaryDocument,JU.OC");
Clazz_defineMethod (c$, "isFullyLit", 
function () {
return (this.params.thePlane != null || this.params.fullyLit);
});
Clazz_defineMethod (c$, "geVdwBitSet", 
function () {
return this.bsVdw;
});
Clazz_defineMethod (c$, "fillAtomData", 
function (atomData, mode) {
if ((mode & 2) != 0 && atomData.bsSelected != null) {
if (this.bsVdw == null) this.bsVdw =  new JU.BS ();
this.bsVdw.or (atomData.bsSelected);
}this.atomDataServer.fillAtomData (atomData, mode);
}, "J.atomdata.AtomData,~N");
Clazz_defineMethod (c$, "getSpanningVectors", 
function () {
return this.surfaceReader.getSpanningVectors ();
});
});
Clazz_declarePackage ("J.jvxl.readers");
Clazz_load (null, "J.jvxl.readers.Parameters", ["java.lang.Float", "java.util.Hashtable", "JU.A4", "$.Lst", "$.M3", "$.P3", "$.P4", "$.V3", "JU.Escape", "$.Logger"], function () {
c$ = Clazz_decorateAsClass (function () {
this.state = 0;
this.testFlags = 0;
this.logMessages = false;
this.logCompression = false;
this.logCube = false;
this.isSilent = false;
this.assocCutoff = 0.3;
this.dataType = 0;
this.surfaceType = 0;
this.calculationType = "";
this.atomRadiusData = null;
this.addHydrogens = false;
this.solventRadius = 0;
this.solventExtendedAtomRadius = 0;
this.propertySmoothing = false;
this.propertySmoothingPower = 4;
this.envelopeRadius = 0;
this.cavityRadius = 0;
this.isCavity = false;
this.pocket = null;
this.minSet = 0;
this.slabInfo = null;
this.slabPlaneOffset = NaN;
this.theProperty = null;
this.solvent_ptsPerAngstrom = 4;
this.solvent_gridMax = 60;
this.plane_ptsPerAngstrom = 4;
this.plane_gridMax = 81;
this.colorBySign = false;
this.colorByPhase = false;
this.colorBySets = false;
this.colorRgb = 0;
this.colorNeg = 0;
this.colorPos = 0;
this.colorPosLCAO = 0;
this.colorNegLCAO = 0;
this.colorPhase = 0;
this.colorDensity = false;
this.iAddGridPoints = false;
this.atomIndex = 0;
this.isAngstroms = false;
this.scale = 0;
this.scale3d = 0;
this.anisotropy = null;
this.isAnisotropic = false;
this.eccentricityMatrix = null;
this.eccentricityMatrixInverse = null;
this.isEccentric = false;
this.eccentricityScale = 0;
this.eccentricityRatio = 0;
this.aniosU = null;
this.anisoB = null;
this.lcaoType = null;
this.functionInfo = null;
this.psi_n = 2;
this.psi_l = 1;
this.psi_m = 1;
this.psi_Znuc = 1;
this.psi_ptsPerAngstrom = 5;
this.psi_monteCarloCount = 0;
this.mep_gridMax = 40;
this.mep_ptsPerAngstrom = 3;
this.mep_marginAngstroms = 1;
this.mep_calcType = -1;
this.qmOrbitalType = 0;
this.qmOrbitalCount = 0;
this.moData = null;
this.qm_gridMax = 80;
this.qm_ptsPerAngstrom = 10;
this.qm_marginAngstroms = 1;
this.qm_nAtoms = 0;
this.qm_moNumber = 2147483647;
this.qm_moLinearCombination = null;
this.center = null;
this.point = null;
this.distance = 0;
this.allowVolumeRender = false;
this.script = null;
this.bsSelected = null;
this.bsIgnore = null;
this.bsSolvent = null;
this.func = null;
this.title = null;
this.blockCubeData = false;
this.readAllData = false;
this.fileIndex = -1;
this.fileName = null;
this.modelIndex = -1;
this.isXLowToHigh = false;
this.insideOut = false;
this.dataXYReversed = false;
this.cutoff = 3.4028235E38;
this.sigma = 3.4028235E38;
this.cutoffAutomatic = true;
this.isCutoffAbsolute = false;
this.isPositiveOnly = false;
this.rangeAll = false;
this.rangeSelected = false;
this.rangeDefined = false;
this.valueMappedToRed = 0;
this.valueMappedToBlue = 0;
this.mappedDataMin = 0;
this.mappedDataMax = 0;
this.isColorReversed = false;
this.isBicolorMap = false;
this.isSquared = false;
this.isSquaredLinear = false;
this.thePlane = null;
this.isContoured = false;
this.nContours = 0;
this.thisContour = 0;
this.contourFromZero = false;
this.parameters = null;
this.resolution = 0;
this.downsampleFactor = 0;
this.maxSet = 0;
this.contoursDiscrete = null;
this.contourColixes = null;
this.contourIncrements = null;
this.boundingBox = null;
this.bsExcluded = null;
this.contourType = 0;
this.colorSchemeTranslucent = false;
this.colorEncoder = null;
this.usePropertyForColorRange = true;
this.isPeriodic = false;
this.doFullMolecular = false;
this.propertyDistanceMax = 2147483647;
this.randomSeed = 0;
this.fullyLit = false;
this.vertexSource = null;
this.intersection = null;
this.origin = null;
this.steps = null;
this.points = null;
this.volumeData = null;
this.contactPair = null;
this.mapLattice = null;
this.isMapped = false;
this.showTiming = false;
this.pointSize = 0;
Clazz_instantialize (this, arguments);
}, J.jvxl.readers, "Parameters");
Clazz_prepareFields (c$, function () {
this.anisotropy =  Clazz_newFloatArray (3, 0);
});
Clazz_defineMethod (c$, "initialize", 
function () {
this.addHydrogens = false;
this.allowVolumeRender = true;
this.atomRadiusData = null;
this.atomIndex = -1;
this.blockCubeData = false;
this.boundingBox = null;
this.bsExcluded =  new Array (4);
this.bsIgnore = null;
this.bsSelected = null;
this.bsSolvent = null;
this.calculationType = "";
this.center = JU.P3.new3 (3.4028235E38, 3.4028235E38, 3.4028235E38);
this.colorBySign = this.colorByPhase = this.colorBySets = false;
this.colorDensity = false;
this.colorEncoder = null;
this.colorNeg = -65536;
this.colorNegLCAO = -8388480;
this.colorPos = -16776961;
this.colorPosLCAO = -23296;
this.colorRgb = -2147483648;
this.colorSchemeTranslucent = false;
this.contactPair = null;
this.contourIncrements = null;
this.contoursDiscrete = null;
this.contourColixes = null;
this.contourFromZero = true;
this.cutoff = 3.4028235E38;
this.cutoffAutomatic = true;
this.dataXYReversed = false;
this.distance = 3.4028235E38;
this.doFullMolecular = false;
this.envelopeRadius = 10;
this.fileIndex = 1;
this.readAllData = true;
this.fileName = "";
this.fullyLit = false;
this.func = null;
this.functionInfo = null;
this.iAddGridPoints = false;
this.insideOut = false;
this.intersection = null;
this.isAngstroms = false;
this.isBicolorMap = this.isCutoffAbsolute = this.isPositiveOnly = false;
this.isCavity = false;
this.isColorReversed = false;
this.isSquared = false;
this.isSquaredLinear = false;
this.isContoured = false;
this.isEccentric = this.isAnisotropic = false;
this.isMapped = false;
this.isPeriodic = false;
this.isSilent = false;
this.mapLattice = null;
this.logCube = this.logCompression = false;
this.logMessages = JU.Logger.debugging;
this.mappedDataMin = 3.4028235E38;
this.mep_calcType = -1;
this.minSet = 0;
this.modelIndex = -1;
this.nContours = 0;
this.pocket = null;
this.pointSize = NaN;
this.propertyDistanceMax = 2147483647;
this.propertySmoothing = false;
this.propertySmoothingPower = 4;
this.rangeDefined = false;
this.rangeAll = false;
this.rangeSelected = false;
this.resolution = 3.4028235E38;
this.scale = NaN;
this.scale3d = 0;
this.sigma = NaN;
this.slabInfo = null;
this.solventExtendedAtomRadius = 0;
this.state = 1;
this.testFlags = 0;
this.thePlane = null;
this.theProperty = null;
this.thisContour = -1;
this.title = null;
this.usePropertyForColorRange = true;
this.vertexSource = null;
this.volumeData = null;
});
Clazz_defineMethod (c$, "setAnisotropy", 
function (pt) {
this.anisotropy[0] = pt.x;
this.anisotropy[1] = pt.y;
this.anisotropy[2] = pt.z;
this.isAnisotropic = true;
if (this.center.x == 3.4028235E38) this.center.set (0, 0, 0);
}, "JU.P3");
Clazz_defineMethod (c$, "setEccentricity", 
function (info) {
var ecc = JU.V3.new3 (info.x, info.y, info.z);
var c = (this.scale > 0 ? this.scale : info.w < 0 ? 1 : ecc.length ());
var fab_c = Math.abs (info.w);
ecc.normalize ();
var z = JU.V3.new3 (0, 0, 1);
ecc.add (z);
ecc.normalize ();
if (Float.isNaN (ecc.x)) ecc.set (1, 0, 0);
this.eccentricityMatrixInverse =  new JU.M3 ();
this.eccentricityMatrixInverse.invertM (this.eccentricityMatrix =  new JU.M3 ().setAA (JU.A4.newVA (ecc, 3.141592653589793)));
this.isEccentric = this.isAnisotropic = true;
this.eccentricityScale = c;
this.eccentricityRatio = fab_c;
if (fab_c > 1) this.eccentricityScale *= fab_c;
this.anisotropy[0] = fab_c * c;
this.anisotropy[1] = fab_c * c;
this.anisotropy[2] = c;
if (this.center.x == 3.4028235E38) this.center.set (0, 0, 0);
}, "JU.P4");
Clazz_defineMethod (c$, "setPlane", 
function (plane) {
this.thePlane = plane;
if (this.thePlane.x == 0 && this.thePlane.y == 0 && this.thePlane.z == 0) this.thePlane.z = 1;
this.isContoured = true;
}, "JU.P4");
Clazz_defineMethod (c$, "setSphere", 
function (radius, isGeodesic) {
this.dataType = (isGeodesic ? 74 : 65);
this.distance = radius;
this.setEccentricity (JU.P4.new4 (0, 0, 1, 1));
this.cutoff = 1.4E-45;
this.isCutoffAbsolute = false;
this.isSilent = !this.logMessages;
this.script = this.getScriptParams () + " SPHERE " + radius + ";";
}, "~N,~B");
Clazz_defineMethod (c$, "setEllipsoidP4", 
function (v) {
this.dataType = 66;
this.distance = 1;
this.setEccentricity (v);
this.cutoff = 1.4E-45;
this.isCutoffAbsolute = false;
this.isSilent = !this.logMessages;
}, "JU.P4");
Clazz_defineMethod (c$, "setEllipsoidAF", 
function (bList) {
this.anisoB = bList;
this.dataType = 67;
this.distance = 0.3 * (Float.isNaN (this.scale) ? 1 : this.scale);
this.cutoff = 1.4E-45;
this.isCutoffAbsolute = false;
this.isSilent = !this.logMessages;
if (this.center.x == 3.4028235E38) this.center.set (0, 0, 0);
if (this.resolution == 3.4028235E38) this.resolution = 6;
}, "~A");
Clazz_defineMethod (c$, "setLobe", 
function (v) {
this.dataType = 68;
this.setEccentricity (v);
if (this.cutoff == 3.4028235E38) {
this.cutoff = 0.14;
if (this.isSquared) this.cutoff = this.cutoff * this.cutoff;
}this.isSilent = !this.logMessages;
this.script = this.getScriptParams () + " LOBE {" + v.x + " " + v.y + " " + v.z + " " + v.w + "};";
}, "JU.P4");
Clazz_defineMethod (c$, "getScriptParams", 
 function () {
return " center " + JU.Escape.eP (this.center) + (Float.isNaN (this.scale) ? "" : " scale " + this.scale);
});
Clazz_defineMethod (c$, "setLp", 
function (v) {
this.dataType = 70;
this.setEccentricity (v);
if (this.cutoff == 3.4028235E38) {
this.cutoff = 0.14;
if (this.isSquared) this.cutoff = this.cutoff * this.cutoff;
}this.isSilent = !this.logMessages;
this.script = " center " + JU.Escape.eP (this.center) + (Float.isNaN (this.scale) ? "" : " scale " + this.scale) + " LP {" + v.x + " " + v.y + " " + v.z + " " + v.w + "};";
}, "JU.P4");
Clazz_defineMethod (c$, "setRadical", 
function (v) {
this.dataType = 71;
this.setEccentricity (v);
if (this.cutoff == 3.4028235E38) {
this.cutoff = 0.14;
if (this.isSquared) this.cutoff = this.cutoff * this.cutoff;
}this.isSilent = !this.logMessages;
this.script = " center " + JU.Escape.eP (this.center) + (Float.isNaN (this.scale) ? "" : " scale " + this.scale) + " RAD {" + v.x + " " + v.y + " " + v.z + " " + v.w + "};";
}, "JU.P4");
Clazz_defineMethod (c$, "setLcao", 
function (type, colorPtr) {
this.lcaoType = type;
if (colorPtr == 1) this.colorPosLCAO = this.colorNegLCAO;
this.isSilent = !this.logMessages;
}, "~S,~N");
Clazz_defineMethod (c$, "setSolvent", 
function (propertyName, radius) {
this.isEccentric = this.isAnisotropic = false;
this.solventRadius = Math.abs (radius);
this.dataType = (this.intersection != null ? 1333 : "nomap" === propertyName ? 1205 : "molecular" === propertyName ? 1203 : "sasurface" === propertyName || this.solventRadius == 0 ? 1196 : 1195);
if (this.state < 2 && (this.cutoffAutomatic || !this.colorDensity) && (this.intersection == null || this.cutoff == 3.4028235E38)) this.cutoff = 0.0;
switch (this.dataType) {
case 1333:
this.calculationType = "VDW intersection";
break;
case 1205:
this.calculationType = "unmapped plane";
break;
case 1203:
this.calculationType = "molecular surface with radius " + this.solventRadius;
break;
case 1195:
this.calculationType = "solvent-excluded surface with radius " + this.solventRadius;
break;
case 1196:
this.calculationType = "solvent-accessible surface with radius " + this.solventRadius;
break;
}
switch (this.dataType) {
case 1205:
this.solventExtendedAtomRadius = this.solventRadius;
this.solventRadius = 0;
this.isContoured = false;
break;
case 1203:
this.solventExtendedAtomRadius = 0;
break;
case 1195:
this.solventExtendedAtomRadius = 0;
if (this.bsIgnore == null) this.bsIgnore = this.bsSolvent;
 else if (this.bsSolvent != null) this.bsIgnore.or (this.bsSolvent);
break;
case 1196:
this.solventExtendedAtomRadius = this.solventRadius;
this.solventRadius = 0;
if (this.bsIgnore == null) this.bsIgnore = this.bsSolvent;
 else if (this.bsSolvent != null) this.bsIgnore.or (this.bsSolvent);
break;
}
}, "~S,~N");
Clazz_defineMethod (c$, "setFunctionXY", 
function (value) {
this.dataType = 8;
this.functionInfo = value;
this.cutoff = 1.4E-45;
this.isEccentric = this.isAnisotropic = false;
}, "JU.Lst");
Clazz_defineMethod (c$, "setFunctionXYZ", 
function (value) {
this.dataType = 9;
this.functionInfo = value;
if (this.cutoff == 3.4028235E38) this.cutoff = 1.4E-45;
this.isEccentric = this.isAnisotropic = false;
}, "JU.Lst");
Clazz_defineMethod (c$, "setAtomicOrbital", 
function (nlmZprs) {
this.dataType = 1294;
this.setEccentricity (JU.P4.new4 (0, 0, 1, 1));
this.psi_n = Clazz_floatToInt (nlmZprs[0]);
this.psi_l = Clazz_floatToInt (nlmZprs[1]);
this.psi_m = Clazz_floatToInt (nlmZprs[2]);
this.psi_Znuc = nlmZprs[3];
this.psi_monteCarloCount = Clazz_floatToInt (nlmZprs[4]);
this.distance = nlmZprs[5];
if (this.distance != 0 || this.thePlane != null) this.allowVolumeRender = false;
this.randomSeed = Clazz_floatToInt (nlmZprs[6]);
this.psi_ptsPerAngstrom = 10;
if (this.cutoff == 3.4028235E38 || this.cutoff == 0.14) {
this.cutoff = (this.psi_monteCarloCount > 0 ? 0 : 0.04);
if (this.isSquared) this.cutoff = this.cutoff * this.cutoff;
}this.isCutoffAbsolute = true;
if (this.state < 2 && this.thePlane == null && this.colorBySign) this.isBicolorMap = true;
return (this.psi_Znuc > 0 && Math.abs (this.psi_m) <= this.psi_l && this.psi_l < this.psi_n);
}, "~A");
Clazz_defineMethod (c$, "setMep", 
function (charges, isMLP) {
this.dataType = (isMLP ? 1329 : 1328);
this.theProperty = charges;
this.usePropertyForColorRange = false;
this.isEccentric = this.isAnisotropic = false;
if (this.cutoff == 3.4028235E38) {
this.cutoff = 0.1;
if (this.isSquared) this.cutoff = this.cutoff * this.cutoff;
}this.isCutoffAbsolute = (this.cutoff > 0 && !this.isPositiveOnly);
this.contourFromZero = false;
if (this.state >= 2 || this.thePlane != null) {
if (!this.rangeDefined && !this.rangeAll) {
this.valueMappedToRed = -0.1;
this.valueMappedToBlue = 0.1;
this.rangeDefined = true;
}} else {
this.colorBySign = true;
this.isBicolorMap = true;
}}, "~A,~B");
Clazz_defineMethod (c$, "setNci", 
function (isPromolecular) {
this.fullyLit = true;
this.qm_gridMax = 200;
if (isPromolecular) this.dataType = 1844;
this.qm_marginAngstroms = 2;
this.qmOrbitalType = (isPromolecular ? 3 : 4);
if (isPromolecular) {
if (this.parameters == null || this.parameters.length < 2) this.parameters = [this.cutoff, 2];
}if (this.cutoff == 3.4028235E38 || this.cutoff == 0) this.cutoff = 0.3;
if (this.isSquared) this.cutoff *= this.cutoff;
if (this.title == null) this.title =  new Array (0);
this.moData =  new java.util.Hashtable ();
}, "~B");
Clazz_defineMethod (c$, "setMO", 
function (iMo, linearCombination) {
this.qm_moLinearCombination = linearCombination;
this.qm_moNumber = (linearCombination == null ? Math.abs (iMo) : Clazz_floatToInt (linearCombination[1]));
this.qmOrbitalType = (this.moData.containsKey ("haveVolumeData") ? 5 : this.moData.containsKey ("gaussians") ? 1 : this.moData.containsKey ("slaters") ? 2 : 0);
var isElectronDensity = (iMo <= 0 && linearCombination == null);
if (this.qmOrbitalType == 0) {
JU.Logger.error ("MO ERROR: No basis functions found in file for MO calculation. (GAUSSIAN 'gfprint' keyword may be missing?)");
this.title = ["no basis functions found in file"];
} else {
var mos = this.moData.get ("mos");
this.qmOrbitalCount = mos.size ();
this.calculationType = this.moData.get ("calculationType");
this.calculationType = "Molecular orbital #" + this.qm_moNumber + "/" + this.qmOrbitalCount + " " + (this.calculationType == null ? "" : this.calculationType);
if (!isElectronDensity) {
if (this.title == null) {
this.title =  new Array (5);
this.title[0] = "%F";
this.title[1] = "Model %M  MO %I/%N %T";
this.title[2] = "?Energy = %E %U";
this.title[3] = "?Symmetry = %S";
this.title[4] = "?Occupancy = %O";
}}}this.dataType = 1837;
if (this.cutoff == 3.4028235E38) {
this.cutoff = (isElectronDensity ? 0.01 : 0.05);
}if (this.isSquared || this.isSquaredLinear) this.cutoff = this.cutoff * this.cutoff;
this.isEccentric = this.isAnisotropic = false;
this.isCutoffAbsolute = (this.cutoff > 0 && !this.isPositiveOnly);
if (this.state >= 2 || this.thePlane != null) return;
this.colorBySign = true;
if (this.colorByPhase && this.colorPhase == 0) this.colorByPhase = false;
this.isBicolorMap = true;
}, "~N,~A");
Clazz_defineMethod (c$, "setMapRanges", 
function (surfaceReader, haveData) {
if (!this.colorDensity) if (this.colorByPhase || this.colorBySign || (this.thePlane != null || this.isBicolorMap) && !this.isContoured) {
this.mappedDataMin = -1;
this.mappedDataMax = 1;
}if (this.mappedDataMin == 3.4028235E38 || this.mappedDataMin == this.mappedDataMax) {
var minMax = surfaceReader.getMinMaxMappedValues (haveData);
this.mappedDataMin = minMax[0];
this.mappedDataMax = minMax[1];
}if (this.mappedDataMin == 0 && this.mappedDataMax == 0) {
this.mappedDataMin = -1;
this.mappedDataMax = 1;
}if (!this.rangeDefined) {
this.valueMappedToRed = this.mappedDataMin;
this.valueMappedToBlue = this.mappedDataMax;
}}, "J.jvxl.readers.SurfaceReader,~B");
Clazz_defineMethod (c$, "resetForMapping", 
function (haveSurface) {
if (!haveSurface) this.state = 2;
this.isMapped = true;
this.qmOrbitalType = 0;
this.parameters = null;
this.colorDensity = false;
this.mappedDataMin = 3.4028235E38;
this.intersection = null;
this.func = null;
this.points = null;
this.origin = null;
this.steps = null;
this.volumeData = null;
this.center.x = 3.4028235E38;
this.isAnisotropic = false;
}, "~B");
Clazz_defineMethod (c$, "addSlabInfo", 
function (slabObject) {
if (this.slabInfo == null) this.slabInfo =  new JU.Lst ();
this.slabInfo.addLast (slabObject);
}, "~A");
Clazz_defineStatics (c$,
"STATE_UNINITIALIZED", 0,
"STATE_INITIALIZED", 1,
"STATE_DATA_READ", 2,
"STATE_DATA_COLORED", 3,
"NO_ANISOTROPY", 32,
"IS_SILENT", 64,
"IS_SOLVENTTYPE", 128,
"HAS_MAXGRID", 256,
"IS_POINTMAPPABLE", 512,
"IS_SLABBABLE", 1024,
"SURFACE_NONE", 0,
"SURFACE_SPHERE", 65,
"SURFACE_ELLIPSOID2", 66,
"SURFACE_ELLIPSOID3", 67,
"SURFACE_LOBE", 68,
"SURFACE_LCAOCARTOON", 69,
"SURFACE_LONEPAIR", 70,
"SURFACE_RADICAL", 71,
"SURFACE_FUNCTIONXY", 8,
"SURFACE_FUNCTIONXYZ", 9,
"SURFACE_GEODESIC", 74,
"SURFACE_SOLVENT", 1195,
"SURFACE_SASURFACE", 1196,
"SURFACE_MOLECULARORBITAL", 1837,
"SURFACE_ATOMICORBITAL", 1294,
"SURFACE_MEP", 1328,
"SURFACE_MLP", 1329,
"SURFACE_MOLECULAR", 1203,
"SURFACE_NCI", 1844,
"SURFACE_INTERSECT", 1333,
"SURFACE_NOMAP", 1205,
"SURFACE_PROPERTY", 1206,
"ANGSTROMS_PER_BOHR", 0.5291772,
"defaultEdgeFractionBase", 35,
"defaultEdgeFractionRange", 90,
"defaultColorFractionBase", 35,
"defaultColorFractionRange", 90,
"defaultMappedDataMin", 0,
"defaultMappedDataMax", 1.0,
"defaultCutoff", 0.02,
"defaultOrbitalCutoff", 0.04,
"defaultLobeCutoff", 0.14,
"defaultOrbitalCutoffOld", 0.14,
"defaultQMOrbitalCutoff", 0.050,
"defaultQMElectronDensityCutoff", 0.010,
"defaultContourCount", 11,
"nContourMax", 100,
"defaultColorNegative", 0xFFFF0000,
"defaultColorPositive", 0xFF0000FF,
"defaultColorNegativeLCAO", 0xFF800080,
"defaultColorPositiveLCAO", 0xFFFFA500,
"defaultSolventRadius", 1.2,
"defaultMepCutoff", 0.1,
"defaultMepMin", -0.1,
"defaultMepMax", 0.1,
"MEP_MAX_GRID", 40,
"QM_TYPE_UNKNOWN", 0,
"QM_TYPE_GAUSSIAN", 1,
"QM_TYPE_SLATER", 2,
"QM_TYPE_NCI_PRO", 3,
"QM_TYPE_NCI_SCF", 4,
"QM_TYPE_VOLUME_DATA", 5,
"MO_MAX_GRID", 80);
});
Clazz_declarePackage ("J.jvxl.readers");
Clazz_load (["J.jvxl.api.VertexDataServer", "JU.P3"], "J.jvxl.readers.SurfaceReader", ["java.lang.Float", "JU.AU", "$.BS", "J.jvxl.calc.MarchingCubes", "$.MarchingSquares", "J.jvxl.data.JvxlCoder", "$.MeshData", "JU.BoxInfo", "$.C", "$.ColorEncoder", "$.Escape", "$.Logger"], function () {
c$ = Clazz_decorateAsClass (function () {
this.sg = null;
this.meshDataServer = null;
this.params = null;
this.meshData = null;
this.jvxlData = null;
this.volumeData = null;
this.edgeData = null;
this.haveSurfaceAtoms = false;
this.allowSigma = false;
this.isProgressive = false;
this.isXLowToHigh = false;
this.assocCutoff = 0.3;
this.isQuiet = false;
this.vertexDataOnly = false;
this.hasColorData = false;
this.dataMin = 3.4028235E38;
this.dataMax = -3.4028235E38;
this.dataMean = 0;
this.xyzMin = null;
this.xyzMax = null;
this.center = null;
this.anisotropy = null;
this.isAnisotropic = false;
this.eccentricityMatrix = null;
this.eccentricityMatrixInverse = null;
this.isEccentric = false;
this.eccentricityScale = 0;
this.eccentricityRatio = 0;
this.edgeCount = 0;
this.volumetricOrigin = null;
this.volumetricVectors = null;
this.voxelCounts = null;
this.voxelData = null;
this.nBytes = 0;
this.nDataPoints = 0;
this.nPointsX = 0;
this.nPointsY = 0;
this.nPointsZ = 0;
this.isJvxl = false;
this.edgeFractionBase = 0;
this.edgeFractionRange = 0;
this.colorFractionBase = 0;
this.colorFractionRange = 0;
this.jvxlFileHeaderBuffer = null;
this.fractionData = null;
this.jvxlEdgeDataRead = "";
this.jvxlColorDataRead = "";
this.jvxlVoxelBitSet = null;
this.jvxlDataIsColorMapped = false;
this.jvxlDataIsPrecisionColor = false;
this.jvxlDataIs2dContour = false;
this.jvxlDataIsColorDensity = false;
this.jvxlCutoff = 0;
this.jvxlNSurfaceInts = 0;
this.cJvxlEdgeNaN = '\0';
this.contourVertexCount = 0;
this.marchingSquares = null;
this.marchingCubes = null;
this.yzPlanes = null;
this.yzCount = 0;
this.qpc = null;
this.ptTemp = null;
this.minMax = null;
this.haveSetAnisotropy = false;
Clazz_instantialize (this, arguments);
}, J.jvxl.readers, "SurfaceReader", null, J.jvxl.api.VertexDataServer);
Clazz_prepareFields (c$, function () {
this.ptTemp =  new JU.P3 ();
});
Clazz_makeConstructor (c$, 
function () {
});
Clazz_defineMethod (c$, "initSR", 
function (sg) {
this.sg = sg;
this.params = sg.getParams ();
this.marchingSquares = sg.getMarchingSquares ();
this.assocCutoff = this.params.assocCutoff;
this.isXLowToHigh = this.params.isXLowToHigh;
this.center = this.params.center;
this.anisotropy = this.params.anisotropy;
this.isAnisotropic = this.params.isAnisotropic;
this.eccentricityMatrix = this.params.eccentricityMatrix;
this.eccentricityMatrixInverse = this.params.eccentricityMatrixInverse;
this.isEccentric = this.params.isEccentric;
this.eccentricityScale = this.params.eccentricityScale;
this.eccentricityRatio = this.params.eccentricityRatio;
this.meshData = sg.getMeshData ();
this.jvxlData = sg.getJvxlData ();
this.setVolumeDataV (sg.getVolumeData ());
this.meshDataServer = sg.getMeshDataServer ();
this.cJvxlEdgeNaN = String.fromCharCode (125);
}, "J.jvxl.readers.SurfaceGenerator");
Clazz_defineMethod (c$, "setOutputChannel", 
function (out) {
}, "JU.OC");
Clazz_defineMethod (c$, "newVoxelDataCube", 
function () {
this.volumeData.setVoxelDataAsArray (this.voxelData =  Clazz_newFloatArray (this.nPointsX, this.nPointsY, this.nPointsZ, 0));
});
Clazz_defineMethod (c$, "setVolumeDataV", 
function (v) {
this.nBytes = 0;
this.volumetricOrigin = v.volumetricOrigin;
this.volumetricVectors = v.volumetricVectors;
this.voxelCounts = v.voxelCounts;
this.voxelData = v.getVoxelData ();
this.volumeData = v;
}, "J.jvxl.data.VolumeData");
Clazz_defineMethod (c$, "jvxlUpdateInfo", 
function () {
this.jvxlData.jvxlUpdateInfo (this.params.title, this.nBytes);
});
Clazz_defineMethod (c$, "readAndSetVolumeParameters", 
function (isMapData) {
if (!this.readVolumeParameters (isMapData)) return false;
if (this.vertexDataOnly) return true;
return (this.volumeData.setUnitVectors ());
}, "~B");
Clazz_defineMethod (c$, "createIsosurface", 
function (justForPlane) {
this.resetIsosurface ();
if (this.params.showTiming) JU.Logger.startTimer ("isosurface creation");
this.jvxlData.cutoff = NaN;
if (!this.readAndSetVolumeParameters (justForPlane)) return false;
if (!justForPlane && !Float.isNaN (this.params.sigma) && !this.allowSigma) {
if (this.params.sigma > 0) JU.Logger.error ("Reader does not support SIGMA option -- using cutoff 1.6");
this.params.cutoff = 1.6;
}if (this.params.sigma < 0) this.params.sigma = -this.params.sigma;
this.nPointsX = this.voxelCounts[0];
this.nPointsY = this.voxelCounts[1];
this.nPointsZ = this.voxelCounts[2];
this.jvxlData.isSlabbable = ((this.params.dataType & 1024) != 0);
this.jvxlData.insideOut = this.params.insideOut;
this.jvxlData.dataXYReversed = this.params.dataXYReversed;
this.jvxlData.isBicolorMap = this.params.isBicolorMap;
this.jvxlData.nPointsX = this.nPointsX;
this.jvxlData.nPointsY = this.nPointsY;
this.jvxlData.nPointsZ = this.nPointsZ;
this.jvxlData.jvxlVolumeDataXml = this.volumeData.xmlData;
this.jvxlData.voxelVolume = this.volumeData.voxelVolume;
if (justForPlane) {
this.volumeData.setMappingPlane (this.params.thePlane);
if (this.meshDataServer != null) this.meshDataServer.fillMeshData (this.meshData, 1, null);
this.params.setMapRanges (this, false);
this.generateSurfaceData ();
this.volumeData.setMappingPlane (null);
} else {
if (!this.readVolumeData (false)) return false;
this.generateSurfaceData ();
}if (this.jvxlFileHeaderBuffer != null) {
var s = this.jvxlFileHeaderBuffer.toString ();
var i = s.indexOf ('\n', s.indexOf ('\n', s.indexOf ('\n') + 1) + 1) + 1;
this.jvxlData.jvxlFileTitle = s.substring (0, i);
}if (this.params.contactPair == null) this.setBBoxAll ();
if (!this.params.isSilent) JU.Logger.info ("boundbox corners " + JU.Escape.eP (this.xyzMin) + " " + JU.Escape.eP (this.xyzMax));
this.jvxlData.boundingBox = [this.xyzMin, this.xyzMax];
this.jvxlData.dataMin = this.dataMin;
this.jvxlData.dataMax = this.dataMax;
this.jvxlData.cutoff = (this.isJvxl ? this.jvxlCutoff : this.params.cutoff);
this.jvxlData.isCutoffAbsolute = this.params.isCutoffAbsolute;
this.jvxlData.pointsPerAngstrom = 1 / this.volumeData.volumetricVectorLengths[0];
this.jvxlData.jvxlColorData = "";
this.jvxlData.jvxlPlane = this.params.thePlane;
this.jvxlData.jvxlEdgeData = this.edgeData;
this.jvxlData.isBicolorMap = this.params.isBicolorMap;
this.jvxlData.isContoured = this.params.isContoured;
this.jvxlData.colorDensity = this.params.colorDensity;
this.jvxlData.pointSize = this.params.pointSize;
if (this.jvxlData.vContours != null) this.params.nContours = this.jvxlData.vContours.length;
this.jvxlData.nContours = (this.params.contourFromZero ? this.params.nContours : -1 - this.params.nContours);
this.jvxlData.nEdges = this.edgeCount;
this.jvxlData.edgeFractionBase = this.edgeFractionBase;
this.jvxlData.edgeFractionRange = this.edgeFractionRange;
this.jvxlData.colorFractionBase = this.colorFractionBase;
this.jvxlData.colorFractionRange = this.colorFractionRange;
this.jvxlData.jvxlDataIs2dContour = this.jvxlDataIs2dContour;
this.jvxlData.jvxlDataIsColorMapped = this.jvxlDataIsColorMapped;
this.jvxlData.jvxlDataIsColorDensity = this.jvxlDataIsColorDensity;
this.jvxlData.isXLowToHigh = this.isXLowToHigh;
this.jvxlData.vertexDataOnly = this.vertexDataOnly;
this.jvxlData.saveVertexCount = 0;
if (this.jvxlDataIsColorMapped || this.jvxlData.nVertexColors > 0) {
if (this.meshDataServer != null) {
this.meshDataServer.fillMeshData (this.meshData, 1, null);
this.meshDataServer.fillMeshData (this.meshData, 2, null);
}this.jvxlData.jvxlColorData = this.readColorData ();
this.updateSurfaceData ();
if (this.meshDataServer != null) this.meshDataServer.notifySurfaceMappingCompleted ();
}if (this.params.showTiming) JU.Logger.checkTimer ("isosurface creation", false);
return true;
}, "~B");
Clazz_defineMethod (c$, "resetIsosurface", 
function () {
this.meshData =  new J.jvxl.data.MeshData ();
this.xyzMin = this.xyzMax = null;
this.jvxlData.isBicolorMap = this.params.isBicolorMap;
if (this.meshDataServer != null) this.meshDataServer.fillMeshData (null, 0, null);
this.contourVertexCount = 0;
if (this.params.cutoff == 3.4028235E38) this.params.cutoff = 0.02;
this.jvxlData.jvxlSurfaceData = "";
this.jvxlData.jvxlEdgeData = "";
this.jvxlData.jvxlColorData = "";
this.edgeCount = 0;
this.edgeFractionBase = 35;
this.edgeFractionRange = 90;
this.colorFractionBase = 35;
this.colorFractionRange = 90;
this.params.mappedDataMin = 3.4028235E38;
});
Clazz_defineMethod (c$, "discardTempData", 
function (discardAll) {
this.discardTempDataSR (discardAll);
}, "~B");
Clazz_defineMethod (c$, "discardTempDataSR", 
function (discardAll) {
if (!discardAll) return;
this.voxelData = null;
this.sg.setMarchingSquares (this.marchingSquares = null);
this.marchingCubes = null;
}, "~B");
Clazz_defineMethod (c$, "initializeVolumetricData", 
function () {
this.nPointsX = this.voxelCounts[0];
this.nPointsY = this.voxelCounts[1];
this.nPointsZ = this.voxelCounts[2];
this.setVolumeDataV (this.volumeData);
});
Clazz_defineMethod (c$, "gotoAndReadVoxelData", 
function (isMapData) {
this.initializeVolumetricData ();
if (this.nPointsX > 0 && this.nPointsY > 0 && this.nPointsZ > 0) try {
this.gotoData (this.params.fileIndex - 1, this.nPointsX * this.nPointsY * this.nPointsZ);
this.readSurfaceData (isMapData);
} catch (e) {
if (Clazz_exceptionOf (e, Exception)) {
JU.Logger.error (e.toString ());
return false;
} else {
throw e;
}
}
return true;
}, "~B");
Clazz_defineMethod (c$, "gotoData", 
function (n, nPoints) {
}, "~N,~N");
Clazz_defineMethod (c$, "readColorData", 
function () {
if (this.jvxlData.vertexColors == null) return "";
var vertexCount = this.jvxlData.vertexCount;
var colixes = this.meshData.vcs;
var vertexValues = this.meshData.vvs;
if (colixes == null || colixes.length < vertexCount) this.meshData.vcs = colixes =  Clazz_newShortArray (vertexCount, 0);
if (vertexValues == null || vertexValues.length < vertexCount) this.meshData.vvs = vertexValues =  Clazz_newFloatArray (vertexCount, 0);
for (var i = 0; i < vertexCount; i++) colixes[i] = JU.C.getColix (this.jvxlData.vertexColors[i]);

return "-";
});
Clazz_overrideMethod (c$, "getPlane", 
function (x) {
return this.getPlane2 (x);
}, "~N");
Clazz_defineMethod (c$, "getPlane2", 
function (x) {
if (this.yzCount == 0) this.initPlanes ();
if (this.qpc != null) this.qpc.getPlane (x, this.yzPlanes[x % 2]);
return this.yzPlanes[x % 2];
}, "~N");
Clazz_defineMethod (c$, "initPlanes", 
function () {
this.yzCount = this.nPointsY * this.nPointsZ;
if (!this.isQuiet) JU.Logger.info ("reading data progressively -- yzCount = " + this.yzCount);
this.yzPlanes = JU.AU.newFloat2 (2);
this.yzPlanes[0] =  Clazz_newFloatArray (this.yzCount, 0);
this.yzPlanes[1] =  Clazz_newFloatArray (this.yzCount, 0);
});
Clazz_overrideMethod (c$, "getValue", 
function (x, y, z, ptyz) {
return this.getValue2 (x, y, z, ptyz);
}, "~N,~N,~N,~N");
Clazz_defineMethod (c$, "getValue2", 
function (x, y, z, ptyz) {
return (this.yzPlanes == null ? this.voxelData[x][y][z] : this.yzPlanes[x % 2][ptyz]);
}, "~N,~N,~N,~N");
Clazz_defineMethod (c$, "generateSurfaceData", 
 function () {
this.edgeData = "";
if (this.vertexDataOnly) {
try {
this.readSurfaceData (false);
} catch (e) {
if (Clazz_exceptionOf (e, Exception)) {
System.out.println (e.toString ());
JU.Logger.error ("Exception in SurfaceReader::readSurfaceData: " + e.toString ());
} else {
throw e;
}
}
return;
}this.contourVertexCount = 0;
var contourType = -1;
this.marchingSquares = null;
if (this.params.thePlane != null || this.params.isContoured) {
this.marchingSquares =  new J.jvxl.calc.MarchingSquares (this, this.volumeData, this.params.thePlane, this.params.contoursDiscrete, this.params.nContours, this.params.thisContour, this.params.contourFromZero);
contourType = this.marchingSquares.getContourType ();
this.marchingSquares.setMinMax (this.params.valueMappedToRed, this.params.valueMappedToBlue);
}this.params.contourType = contourType;
this.params.isXLowToHigh = this.isXLowToHigh;
this.marchingCubes =  new J.jvxl.calc.MarchingCubes (this, this.volumeData, this.params, this.jvxlVoxelBitSet);
var data = this.marchingCubes.getEdgeData ();
if (this.params.thePlane == null) this.edgeData = data;
this.jvxlData.setSurfaceInfoFromBitSetPts (this.marchingCubes.getBsVoxels (), this.params.thePlane, this.params.mapLattice);
this.jvxlData.jvxlExcluded = this.params.bsExcluded;
if (this.isJvxl) this.edgeData = this.jvxlEdgeDataRead;
this.postProcessVertices ();
});
Clazz_defineMethod (c$, "postProcessVertices", 
function () {
});
Clazz_overrideMethod (c$, "getSurfacePointIndexAndFraction", 
function (cutoff, isCutoffAbsolute, x, y, z, offset, vA, vB, valueA, valueB, pointA, edgeVector, isContourType, fReturn) {
var thisValue = this.getSurfacePointAndFraction (cutoff, isCutoffAbsolute, valueA, valueB, pointA, edgeVector, x, y, z, vA, vB, fReturn, this.ptTemp);
if (this.marchingSquares != null && this.params.isContoured) return this.marchingSquares.addContourVertex (this.ptTemp, cutoff);
var assocVertex = (this.assocCutoff > 0 ? (fReturn[0] < this.assocCutoff ? vA : fReturn[0] > 1 - this.assocCutoff ? vB : -1) : -1);
if (assocVertex >= 0) assocVertex = this.marchingCubes.getLinearOffset (x, y, z, assocVertex);
var n = this.addVertexCopy (this.ptTemp, thisValue, assocVertex, true);
if (n >= 0 && this.params.iAddGridPoints) {
this.marchingCubes.calcVertexPoint (x, y, z, vB, this.ptTemp);
this.addVertexCopy (valueA < valueB ? pointA : this.ptTemp, Math.min (valueA, valueB), -3, true);
this.addVertexCopy (valueA < valueB ? this.ptTemp : pointA, Math.max (valueA, valueB), -3, true);
}return n;
}, "~N,~B,~N,~N,~N,JU.P3i,~N,~N,~N,~N,JU.T3,JU.V3,~B,~A");
Clazz_defineMethod (c$, "getSurfacePointAndFraction", 
function (cutoff, isCutoffAbsolute, valueA, valueB, pointA, edgeVector, x, y, z, vA, vB, fReturn, ptReturn) {
return this.getSPF (cutoff, isCutoffAbsolute, valueA, valueB, pointA, edgeVector, x, y, z, vA, vB, fReturn, ptReturn);
}, "~N,~B,~N,~N,JU.T3,JU.V3,~N,~N,~N,~N,~N,~A,JU.T3");
Clazz_defineMethod (c$, "getSPF", 
function (cutoff, isCutoffAbsolute, valueA, valueB, pointA, edgeVector, x, y, z, vA, vB, fReturn, ptReturn) {
var diff = valueB - valueA;
var fraction = (cutoff - valueA) / diff;
if (isCutoffAbsolute && (fraction < 0 || fraction > 1)) fraction = (-cutoff - valueA) / diff;
if (fraction < 0 || fraction > 1) {
fraction = NaN;
}fReturn[0] = fraction;
ptReturn.scaleAdd2 (fraction, edgeVector, pointA);
return valueA + fraction * diff;
}, "~N,~B,~N,~N,JU.T3,JU.V3,~N,~N,~N,~N,~N,~A,JU.T3");
Clazz_defineMethod (c$, "addVertexCopy", 
function (vertexXYZ, value, assocVertex, asCopy) {
return this.addVC (vertexXYZ, value, assocVertex, asCopy);
}, "JU.T3,~N,~N,~B");
Clazz_defineMethod (c$, "addVC", 
function (vertexXYZ, value, assocVertex, asCopy) {
return (Float.isNaN (value) && assocVertex != -3 ? -1 : this.meshDataServer == null ? this.meshData.addVertexCopy (vertexXYZ, value, assocVertex, asCopy) : this.meshDataServer.addVertexCopy (vertexXYZ, value, assocVertex, asCopy));
}, "JU.T3,~N,~N,~B");
Clazz_defineMethod (c$, "addTriangleCheck", 
function (iA, iB, iC, check, check2, isAbsolute, color) {
if (this.marchingSquares != null && this.params.isContoured) {
if (color == 0) return this.marchingSquares.addTriangle (iA, iB, iC, check, check2);
color = 0;
}return (this.meshDataServer != null ? this.meshDataServer.addTriangleCheck (iA, iB, iC, check, check2, isAbsolute, color) : isAbsolute && !J.jvxl.data.MeshData.checkCutoff (iA, iB, iC, this.meshData.vvs) ? -1 : this.meshData.addTriangleCheck (iA, iB, iC, check, check2, color));
}, "~N,~N,~N,~N,~N,~B,~N");
Clazz_defineMethod (c$, "colorIsosurface", 
function () {
if (this.params.isSquared && this.volumeData != null) this.volumeData.filterData (true, NaN);
if (this.meshDataServer != null) {
this.meshDataServer.fillMeshData (this.meshData, 1, null);
}this.jvxlData.saveVertexCount = 0;
if (this.params.isContoured && this.marchingSquares != null) {
this.initializeMapping ();
this.params.setMapRanges (this, false);
this.marchingSquares.setMinMax (this.params.valueMappedToRed, this.params.valueMappedToBlue);
this.jvxlData.saveVertexCount = this.marchingSquares.contourVertexCount;
this.contourVertexCount = this.marchingSquares.generateContourData (this.jvxlDataIs2dContour, (this.params.isSquared ? 1e-8 : 1e-4));
this.jvxlData.contourValuesUsed = this.marchingSquares.getContourValues ();
this.minMax = this.marchingSquares.getMinMax ();
if (this.meshDataServer != null) this.meshDataServer.notifySurfaceGenerationCompleted ();
this.finalizeMapping ();
}this.applyColorScale ();
this.jvxlData.nContours = (this.params.contourFromZero ? this.params.nContours : -1 - this.params.nContours);
this.jvxlData.jvxlFileMessage = "mapped: min = " + this.params.valueMappedToRed + "; max = " + this.params.valueMappedToBlue;
});
Clazz_defineMethod (c$, "applyColorScale", 
function () {
this.colorFractionBase = this.jvxlData.colorFractionBase = 35;
this.colorFractionRange = this.jvxlData.colorFractionRange = 90;
if (this.params.colorPhase == 0) this.params.colorPhase = 1;
if (this.meshDataServer == null) {
this.meshData.vcs =  Clazz_newShortArray (this.meshData.vc, 0);
} else {
this.meshDataServer.fillMeshData (this.meshData, 1, null);
if (this.params.contactPair == null) this.meshDataServer.fillMeshData (this.meshData, 2, null);
}var saveColorData = (this.params.colorDensity || this.params.isBicolorMap || this.params.colorBySign || !this.params.colorByPhase);
if (this.params.contactPair != null) saveColorData = false;
this.jvxlData.isJvxlPrecisionColor = true;
this.jvxlData.vertexCount = (this.contourVertexCount > 0 ? this.contourVertexCount : this.meshData.vc);
this.jvxlData.minColorIndex = -1;
this.jvxlData.maxColorIndex = 0;
this.jvxlData.contourValues = this.params.contoursDiscrete;
this.jvxlData.isColorReversed = this.params.isColorReversed;
if (!this.params.colorDensity) if (this.params.isBicolorMap && !this.params.isContoured || this.params.colorBySign) {
this.jvxlData.minColorIndex = JU.C.getColixTranslucent3 (JU.C.getColix (this.params.isColorReversed ? this.params.colorPos : this.params.colorNeg), this.jvxlData.translucency != 0, this.jvxlData.translucency);
this.jvxlData.maxColorIndex = JU.C.getColixTranslucent3 (JU.C.getColix (this.params.isColorReversed ? this.params.colorNeg : this.params.colorPos), this.jvxlData.translucency != 0, this.jvxlData.translucency);
}this.jvxlData.isTruncated = (this.jvxlData.minColorIndex >= 0 && !this.params.isContoured);
var useMeshDataValues = this.jvxlDataIs2dContour || this.hasColorData || this.vertexDataOnly || this.params.colorDensity || this.params.isBicolorMap && !this.params.isContoured;
if (!useMeshDataValues) {
if (this.haveSurfaceAtoms && this.meshData.vertexSource == null) this.meshData.vertexSource =  Clazz_newIntArray (this.meshData.vc, 0);
var min = 3.4028235E38;
var max = -3.4028235E38;
var value;
this.initializeMapping ();
for (var i = this.meshData.vc; --i >= this.meshData.mergeVertexCount0; ) {
if (this.params.colorBySets) {
value = this.meshData.vertexSets[i];
} else if (this.params.colorByPhase) {
value = this.getPhase (this.meshData.vs[i]);
} else {
var needSource = this.haveSurfaceAtoms;
value = this.volumeData.lookupInterpolatedVoxelValue (this.meshData.vs[i], needSource);
if (needSource) this.meshData.vertexSource[i] = this.getSurfaceAtomIndex ();
}if (value < min) min = value;
if (value > max && value != 3.4028235E38) max = value;
this.meshData.vvs[i] = value;
}
if (this.params.rangeSelected && this.minMax == null) this.minMax = [min, max];
this.finalizeMapping ();
}this.params.setMapRanges (this, true);
this.jvxlData.mappedDataMin = this.params.mappedDataMin;
this.jvxlData.mappedDataMax = this.params.mappedDataMax;
this.jvxlData.valueMappedToRed = this.params.valueMappedToRed;
this.jvxlData.valueMappedToBlue = this.params.valueMappedToBlue;
if (this.params.contactPair == null && this.jvxlData.vertexColors == null) this.colorData ();
J.jvxl.data.JvxlCoder.jvxlCreateColorData (this.jvxlData, (saveColorData ? this.meshData.vvs : null));
if (this.haveSurfaceAtoms && this.meshDataServer != null) this.meshDataServer.fillMeshData (this.meshData, 4, null);
if (this.meshDataServer != null && this.params.colorBySets) this.meshDataServer.fillMeshData (this.meshData, 3, null);
});
Clazz_defineMethod (c$, "colorData", 
 function () {
var vertexValues = this.meshData.vvs;
var vertexColixes = this.meshData.vcs;
this.meshData.pcs = null;
var valueBlue = this.jvxlData.valueMappedToBlue;
var valueRed = this.jvxlData.valueMappedToRed;
var minColorIndex = this.jvxlData.minColorIndex;
var maxColorIndex = this.jvxlData.maxColorIndex;
if (this.params.colorEncoder == null) this.params.colorEncoder =  new JU.ColorEncoder (null);
this.params.colorEncoder.setRange (this.params.valueMappedToRed, this.params.valueMappedToBlue, this.params.isColorReversed);
for (var i = this.meshData.vc; --i >= 0; ) {
var value = vertexValues[i];
if (minColorIndex >= 0) {
if (value <= 0) vertexColixes[i] = minColorIndex;
 else if (value > 0) vertexColixes[i] = maxColorIndex;
} else {
if (value <= valueRed) value = valueRed;
if (value >= valueBlue) value = valueBlue;
vertexColixes[i] = this.params.colorEncoder.getColorIndex (value);
}}
if ((this.params.nContours > 0 || this.jvxlData.contourValues != null) && this.jvxlData.contourColixes == null) {
var n = (this.jvxlData.contourValues == null ? this.params.nContours : this.jvxlData.contourValues.length);
var colors = this.jvxlData.contourColixes =  Clazz_newShortArray (n, 0);
var values = this.jvxlData.contourValues;
if (values == null) values = this.jvxlData.contourValuesUsed;
if (this.jvxlData.contourValuesUsed == null) this.jvxlData.contourValuesUsed = (values == null ?  Clazz_newFloatArray (n, 0) : values);
var dv = (valueBlue - valueRed) / (n + 1);
this.params.colorEncoder.setRange (this.params.valueMappedToRed, this.params.valueMappedToBlue, this.params.isColorReversed);
for (var i = 0; i < n; i++) {
var v = (values == null ? valueRed + (i + 1) * dv : values[i]);
this.jvxlData.contourValuesUsed[i] = v;
colors[i] = JU.C.getColixTranslucent (this.params.colorEncoder.getArgb (v));
}
this.jvxlData.contourColors = JU.C.getHexCodes (colors);
}});
c$.getColorPhaseIndex = Clazz_defineMethod (c$, "getColorPhaseIndex", 
function (color) {
var colorPhase = -1;
for (var i = 0; i < J.jvxl.readers.SurfaceReader.colorPhases.length; i++) if (color.equalsIgnoreCase (J.jvxl.readers.SurfaceReader.colorPhases[i])) {
colorPhase = i;
break;
}
return colorPhase;
}, "~S");
Clazz_defineMethod (c$, "getPhase", 
 function (pt) {
switch (this.params.colorPhase) {
case 0:
case -1:
case 1:
return (pt.x > 0 ? 1 : -1);
case 2:
return (pt.y > 0 ? 1 : -1);
case 3:
return (pt.z > 0 ? 1 : -1);
case 4:
return (pt.x * pt.y > 0 ? 1 : -1);
case 5:
return (pt.y * pt.z > 0 ? 1 : -1);
case 6:
return (pt.x * pt.z > 0 ? 1 : -1);
case 7:
return (pt.x * pt.x - pt.y * pt.y > 0 ? 1 : -1);
case 8:
return (pt.z * pt.z * 2 - pt.x * pt.x - pt.y * pt.y > 0 ? 1 : -1);
}
return 1;
}, "JU.T3");
Clazz_defineMethod (c$, "getMinMaxMappedValues", 
function (haveData) {
if (this.minMax != null && this.minMax[0] != 3.4028235E38) return this.minMax;
if (this.params.colorBySets) return (this.minMax = [0, Math.max (this.meshData.nSets - 1, 0)]);
var min = 3.4028235E38;
var max = -3.4028235E38;
if (this.params.usePropertyForColorRange && this.params.theProperty != null) {
for (var i = 0; i < this.params.theProperty.length; i++) {
if (this.params.rangeSelected && !this.params.bsSelected.get (i)) continue;
var p = this.params.theProperty[i];
if (Float.isNaN (p)) continue;
if (p < min) min = p;
if (p > max) max = p;
}
return (this.minMax = [min, max]);
}var vertexCount = (this.contourVertexCount > 0 ? this.contourVertexCount : this.meshData.vc);
var vertexes = this.meshData.vs;
var useVertexValue = (haveData || this.jvxlDataIs2dContour || this.vertexDataOnly || this.params.colorDensity);
for (var i = this.meshData.mergeVertexCount0; i < vertexCount; i++) {
var v;
if (useVertexValue) v = this.meshData.vvs[i];
 else v = this.volumeData.lookupInterpolatedVoxelValue (vertexes[i], false);
if (v < min) min = v;
if (v > max && v != 3.4028235E38) max = v;
}
return (this.minMax = [min, max]);
}, "~B");
Clazz_defineMethod (c$, "updateTriangles", 
function () {
if (this.meshDataServer == null) {
this.meshData.invalidatePolygons ();
} else {
this.meshDataServer.invalidateTriangles ();
}});
Clazz_defineMethod (c$, "updateSurfaceData", 
function () {
this.meshData.setVertexSets (true);
this.updateTriangles ();
if (this.params.bsExcluded[1] == null) this.params.bsExcluded[1] =  new JU.BS ();
this.meshData.updateInvalidatedVertices (this.params.bsExcluded[1]);
});
Clazz_defineMethod (c$, "selectPocket", 
function (doExclude) {
}, "~B");
Clazz_defineMethod (c$, "excludeMinimumSet", 
function () {
if (this.meshDataServer != null) this.meshDataServer.fillMeshData (this.meshData, 1, null);
this.meshData.getSurfaceSet ();
var bs;
for (var i = this.meshData.nSets; --i >= 0; ) if ((bs = this.meshData.surfaceSet[i]) != null && bs.cardinality () < this.params.minSet) this.meshData.invalidateSurfaceSet (i);

this.updateSurfaceData ();
if (this.meshDataServer != null) this.meshDataServer.fillMeshData (this.meshData, 3, null);
});
Clazz_defineMethod (c$, "excludeMaximumSet", 
function () {
if (this.meshDataServer != null) this.meshDataServer.fillMeshData (this.meshData, 1, null);
this.meshData.getSurfaceSet ();
var bs;
for (var i = this.meshData.nSets; --i >= 0; ) if ((bs = this.meshData.surfaceSet[i]) != null && bs.cardinality () > this.params.maxSet) this.meshData.invalidateSurfaceSet (i);

this.updateSurfaceData ();
if (this.meshDataServer != null) this.meshDataServer.fillMeshData (this.meshData, 3, null);
});
Clazz_defineMethod (c$, "slabIsosurface", 
function (slabInfo) {
if (this.meshDataServer != null) this.meshDataServer.fillMeshData (this.meshData, 1, null);
this.meshData.slabPolygonsList (slabInfo, true);
if (this.meshDataServer != null) this.meshDataServer.fillMeshData (this.meshData, 4, null);
}, "JU.Lst");
Clazz_defineMethod (c$, "setVertexAnisotropy", 
function (pt) {
pt.x *= this.anisotropy[0];
pt.y *= this.anisotropy[1];
pt.z *= this.anisotropy[2];
pt.add (this.center);
}, "JU.T3");
Clazz_defineMethod (c$, "setVectorAnisotropy", 
function (v) {
this.haveSetAnisotropy = true;
v.x *= this.anisotropy[0];
v.y *= this.anisotropy[1];
v.z *= this.anisotropy[2];
}, "JU.T3");
Clazz_defineMethod (c$, "setVolumetricAnisotropy", 
function () {
if (this.haveSetAnisotropy) return;
this.setVolumetricOriginAnisotropy ();
this.setVectorAnisotropy (this.volumetricVectors[0]);
this.setVectorAnisotropy (this.volumetricVectors[1]);
this.setVectorAnisotropy (this.volumetricVectors[2]);
});
Clazz_defineMethod (c$, "setVolumetricOriginAnisotropy", 
function () {
this.volumetricOrigin.setT (this.center);
});
Clazz_defineMethod (c$, "setBBoxAll", 
 function () {
if (this.meshDataServer != null) this.meshDataServer.fillMeshData (this.meshData, 1, null);
this.xyzMin = null;
for (var i = 0; i < this.meshData.vc; i++) {
var p = this.meshData.vs[i];
if (!Float.isNaN (p.x)) this.setBBox (p, 0);
}
});
Clazz_defineMethod (c$, "setBBox", 
function (pt, margin) {
if (this.xyzMin == null) {
this.xyzMin = JU.P3.new3 (3.4028235E38, 3.4028235E38, 3.4028235E38);
this.xyzMax = JU.P3.new3 (-3.4028235E38, -3.4028235E38, -3.4028235E38);
}JU.BoxInfo.addPoint (pt, this.xyzMin, this.xyzMax, margin);
}, "JU.T3,~N");
Clazz_defineMethod (c$, "getValueAtPoint", 
function (pt, getSource) {
return 0;
}, "JU.T3,~B");
Clazz_defineMethod (c$, "initializeMapping", 
function () {
});
Clazz_defineMethod (c$, "finalizeMapping", 
function () {
});
Clazz_defineMethod (c$, "getSurfaceAtomIndex", 
function () {
return -1;
});
Clazz_defineMethod (c$, "getSpanningVectors", 
function () {
return (this.volumeData == null ? null : this.volumeData.getSpanningVectors ());
});
Clazz_defineStatics (c$,
"ANGSTROMS_PER_BOHR", 0.5291772,
"defaultMappedDataMin", 0,
"defaultMappedDataMax", 1.0,
"defaultCutoff", 0.02,
"colorPhases", ["_orb", "x", "y", "z", "xy", "yz", "xz", "x2-y2", "z2"]);
});
Clazz_declarePackage ("J.jvxl.calc");
Clazz_load (["JU.TriangleData", "JU.BS", "$.P3", "$.SB", "$.V3"], "J.jvxl.calc.MarchingCubes", ["java.lang.Float", "J.jvxl.data.JvxlCoder"], function () {
c$ = Clazz_decorateAsClass (function () {
this.surfaceReader = null;
this.volumeData = null;
this.contourType = 0;
this.isContoured = false;
this.cutoff = 0;
this.isCutoffAbsolute = false;
this.isSquared = false;
this.isXLowToHigh = false;
this.cubeCountX = 0;
this.cubeCountY = 0;
this.cubeCountZ = 0;
this.nY = 0;
this.nZ = 0;
this.yzCount = 0;
this.colorDensity = false;
this.integrateSquared = true;
this.bsVoxels = null;
this.bsExcludedVertices = null;
this.bsExcludedTriangles = null;
this.bsExcludedPlanes = null;
this.edgeData = null;
this.excludePartialCubes = true;
this.mode = 0;
this.vertexValues = null;
this.edgeCount = 0;
this.voxelVertexVectors = null;
this.edgeVectors = null;
this.edgePointIndexes = null;
this.isoPointIndexPlanes = null;
this.yzPlanes = null;
this.mappingPlane = null;
this.allInside = false;
this.$isInside = false;
this.offset = null;
this.voxelData = null;
this.nTriangles = 0;
this.bsValues = null;
this.pt0 = null;
this.pointA = null;
this.edgeVertexPointers = null;
this.edgeVertexPlanes = null;
this.fReturn = null;
this.linearOffsets = null;
Clazz_instantialize (this, arguments);
}, J.jvxl.calc, "MarchingCubes", JU.TriangleData);
Clazz_prepareFields (c$, function () {
this.edgeData =  new JU.SB ();
this.vertexValues =  Clazz_newFloatArray (8, 0);
this.voxelVertexVectors =  new Array (8);
this.edgeVectors =  new Array (12);
{
for (var i = 12; --i >= 0; ) this.edgeVectors[i] =  new JU.V3 ();

}this.edgePointIndexes =  Clazz_newIntArray (12, 0);
this.bsValues =  new JU.BS ();
this.pt0 =  new JU.P3 ();
this.pointA =  new JU.P3 ();
this.fReturn =  Clazz_newFloatArray (1, 0);
this.linearOffsets =  Clazz_newIntArray (8, 0);
});
Clazz_defineMethod (c$, "getBsVoxels", 
function () {
return this.bsVoxels;
});
Clazz_makeConstructor (c$, 
function () {
Clazz_superConstructor (this, J.jvxl.calc.MarchingCubes, []);
});
Clazz_makeConstructor (c$, 
function (surfaceReader, volumeData, params, bsVoxels) {
Clazz_superConstructor (this, J.jvxl.calc.MarchingCubes, []);
this.excludePartialCubes = true;
this.surfaceReader = surfaceReader;
this.bsVoxels = bsVoxels;
var bsExcluded = params.bsExcluded;
this.bsExcludedVertices = (bsExcluded[0] == null ? bsExcluded[0] =  new JU.BS () : bsExcluded[0]);
this.bsExcludedPlanes = (bsExcluded[2] == null ? bsExcluded[2] =  new JU.BS () : bsExcluded[2]);
this.bsExcludedTriangles = (bsExcluded[3] == null ? bsExcluded[3] =  new JU.BS () : bsExcluded[3]);
this.mode = (volumeData.getVoxelData () != null || volumeData.mappingPlane != null ? 1 : bsVoxels != null ? 2 : 3);
this.setParameters (volumeData, params);
}, "J.jvxl.api.VertexDataServer,J.jvxl.data.VolumeData,J.jvxl.readers.Parameters,JU.BS");
Clazz_defineMethod (c$, "setParameters", 
function (volumeData, params) {
this.volumeData = volumeData;
this.colorDensity = params.colorDensity;
this.isContoured = params.thePlane == null && params.isContoured && !this.colorDensity;
this.cutoff = params.cutoff;
this.isCutoffAbsolute = params.isCutoffAbsolute;
this.contourType = params.contourType;
this.isSquared = params.isSquared;
this.isXLowToHigh = params.isXLowToHigh;
this.cubeCountX = volumeData.voxelCounts[0] - 1;
this.cubeCountY = volumeData.voxelCounts[1] - 1;
this.cubeCountZ = volumeData.voxelCounts[2] - 1;
volumeData.getYzCount ();
if (params.mapLattice != null) {
this.cubeCountX *= Math.abs (params.mapLattice.x);
this.cubeCountY *= Math.abs (params.mapLattice.y);
this.cubeCountZ *= Math.abs (params.mapLattice.z);
}this.nY = this.cubeCountY + 1;
this.nZ = this.cubeCountZ + 1;
this.yzCount = this.nY * this.nZ;
if (this.bsVoxels == null) this.bsVoxels =  new JU.BS ();
this.edgeVertexPointers = (this.isXLowToHigh ? J.jvxl.calc.MarchingCubes.edgeVertexPointersLowToHigh : J.jvxl.calc.MarchingCubes.edgeVertexPointersHighToLow);
this.edgeVertexPlanes = (this.isXLowToHigh ? J.jvxl.calc.MarchingCubes.edgeVertexPlanesLowToHigh : J.jvxl.calc.MarchingCubes.edgeVertexPlanesHighToLow);
this.isoPointIndexPlanes =  Clazz_newIntArray (2, this.yzCount, 3, 0);
this.yzPlanes = (this.mode == 3 ?  Clazz_newFloatArray (2, this.yzCount, 0) : null);
this.setLinearOffsets ();
this.calcVoxelVertexVectors ();
}, "J.jvxl.data.VolumeData,J.jvxl.readers.Parameters");
Clazz_defineMethod (c$, "calcVoxelVertexVectors", 
function () {
for (var i = 8; --i >= 0; ) this.volumeData.transform (J.jvxl.calc.MarchingCubes.cubeVertexVectors[i], this.voxelVertexVectors[i] =  new JU.V3 ());

for (var i = 12; --i >= 0; ) this.edgeVectors[i].sub2 (this.voxelVertexVectors[JU.TriangleData.edgeVertexes[i + i + 1]], this.voxelVertexVectors[JU.TriangleData.edgeVertexes[i + i]]);

});
Clazz_defineMethod (c$, "resetIndexPlane", 
function (plane) {
for (var i = 0; i < this.yzCount; i++) for (var j = 0; j < 3; j++) plane[i][j] = -2147483648;


return plane;
}, "~A");
Clazz_defineMethod (c$, "getEdgeData", 
function () {
if (this.cubeCountX < 0 || this.cubeCountY < 0 || this.cubeCountZ < 0) return "";
this.mappingPlane = this.volumeData.mappingPlane;
this.edgeCount = 0;
var x0;
var x1;
var xStep;
var ptStep;
var pt;
var ptX;
if (this.isXLowToHigh) {
x0 = 0;
x1 = this.cubeCountX + (this.colorDensity ? 1 : 0);
if (this.colorDensity) {
x1 = this.cubeCountX + 1;
ptX = this.yzCount - 1;
} else {
x1 = this.cubeCountX;
ptX = (this.yzCount - 1) - this.nZ - 1;
}xStep = 1;
ptStep = this.yzCount;
} else {
if (this.colorDensity) {
x0 = this.cubeCountX;
ptX = (this.cubeCountX + 1) * this.yzCount - 1;
} else {
x0 = this.cubeCountX - 1;
ptX = (this.cubeCountX * this.yzCount - 1) - this.nZ - 1;
}x1 = -1;
xStep = -1;
ptStep = -this.yzCount;
}pt = ptX;
this.resetIndexPlane (this.isoPointIndexPlanes[1]);
this.voxelData = null;
var y1 = this.cubeCountY + (this.colorDensity ? 1 : 0);
var z1 = this.cubeCountZ + (this.colorDensity ? 1 : 0);
switch (this.mode) {
case 3:
this.getPlane (x0, false);
break;
case 1:
this.voxelData = this.volumeData.getVoxelData ();
break;
}
this.allInside = (this.colorDensity && (this.cutoff == 0 || this.mode == 2 && this.bsVoxels.cardinality () == 0));
var colorDensityAll = (this.colorDensity && this.cutoff == 0);
var v = 0;
for (var x = x0; x != x1; x += xStep, ptX += ptStep, pt = ptX) {
if (this.mode == 3) {
if (x + xStep <= x1) this.getPlane (x + xStep, true);
}if (this.bsExcludedPlanes.get (x) && this.bsExcludedPlanes.get (x + xStep)) continue;
if (this.colorDensity) {
for (var y = y1; --y >= 0; ) for (var z = z1; --z >= 0; pt--) {
v = this.getValue (x, y, z, pt, 0);
if (colorDensityAll || this.$isInside) {
this.addVertex (x, y, z, pt, v);
}}

continue;
}var indexPlane = this.isoPointIndexPlanes[0];
this.isoPointIndexPlanes[0] = this.isoPointIndexPlanes[1];
this.isoPointIndexPlanes[1] = this.resetIndexPlane (indexPlane);
var noValues = true;
for (var y = y1; --y >= 0; pt--) {
for (var z = z1; --z >= 0; pt--) {
var insideMask = 0;
for (var i = 8; --i >= 0; ) {
v = this.getValue (x, y, z, pt, i);
if (this.$isInside) insideMask |= JU.TriangleData.Pwr2[i];
}
if (noValues && !Float.isNaN (v)) noValues = false;
if (insideMask == 0) {
continue;
}if (insideMask == 0xFF) {
continue;
}if (this.processOneCubical (insideMask, x, y, z, pt) && !this.isContoured && !this.colorDensity) {
this.processTriangles (insideMask);
}}
}
if (noValues) {
this.bsExcludedPlanes.set (x);
}}
return this.edgeData.toString ();
});
Clazz_defineMethod (c$, "getValue", 
 function (x, y, z, pt, i) {
var v;
this.offset = JU.TriangleData.cubeVertexOffsets[i];
var pti = pt + this.linearOffsets[i];
switch (this.mode) {
case 3:
v = this.vertexValues[i] = this.getValueArray (x + this.offset.x, y + this.offset.y, z + this.offset.z, pti, this.yzPlanes[J.jvxl.calc.MarchingCubes.yzPlanePts[i]]);
this.$isInside = (this.allInside || this.bsVoxels.get (pti));
break;
case 2:
this.$isInside = (this.allInside || this.bsVoxels.get (pti));
v = this.vertexValues[i] = (this.bsExcludedVertices.get (pti) ? NaN : this.$isInside ? 1 : 0);
break;
default:
case 1:
if (this.mappingPlane == null) {
v = this.vertexValues[i] = this.voxelData[x + this.offset.x][y + this.offset.y][z + this.offset.z];
} else {
this.volumeData.voxelPtToXYZ (x + this.offset.x, y + this.offset.y, z + this.offset.z, this.pt0);
v = this.vertexValues[i] = this.volumeData.distanceToMappingPlane (this.pt0);
}if (this.isSquared) this.vertexValues[i] *= this.vertexValues[i];
this.$isInside = (this.allInside ? true : J.jvxl.calc.MarchingCubes.isInside (this.vertexValues[i], this.cutoff, this.isCutoffAbsolute));
if (this.$isInside) this.bsVoxels.set (pti);
}
return v;
}, "~N,~N,~N,~N,~N");
Clazz_defineMethod (c$, "getPlane", 
 function (i, andSwap) {
if (i < 0 || i > this.cubeCountX) return;
this.surfaceReader.getPlane (i);
if (andSwap) {
var plane = this.yzPlanes[0];
this.yzPlanes[0] = this.yzPlanes[1];
this.yzPlanes[1] = plane;
}}, "~N,~B");
Clazz_defineMethod (c$, "processTriangles", 
function (insideMask) {
var triangles = JU.TriangleData.triangleTable2[insideMask];
for (var i = triangles.length; (i -= 4) >= 0; ) this.addTriangle (triangles[i], triangles[i + 1], triangles[i + 2], triangles[i + 3]);

}, "~N");
Clazz_defineMethod (c$, "addVertex", 
function (x, y, z, pti, value) {
this.volumeData.voxelPtToXYZ (x, y, z, this.pt0);
if (this.surfaceReader.addVertexCopy (this.pt0, value, -4, true) < 0) this.bsExcludedVertices.set (pti);
}, "~N,~N,~N,~N,~N");
Clazz_defineMethod (c$, "addTriangle", 
function (ia, ib, ic, edgeType) {
if (!this.bsExcludedTriangles.get (this.nTriangles) && this.surfaceReader.addTriangleCheck (this.edgePointIndexes[ia], this.edgePointIndexes[ib], this.edgePointIndexes[ic], edgeType, 0, this.isCutoffAbsolute, 0) < 0) {
this.bsExcludedTriangles.set (this.nTriangles);
}this.nTriangles++;
}, "~N,~N,~N,~N");
Clazz_defineMethod (c$, "getValueArray", 
function (x, y, z, pt, tempValues) {
var ptyz = pt % this.yzCount;
this.bsValues.set (pt);
var value = this.surfaceReader.getValue (x, y, z, ptyz);
if (this.isSquared) value *= value;
tempValues[ptyz] = value;
if (J.jvxl.calc.MarchingCubes.isInside (value, this.cutoff, this.isCutoffAbsolute)) this.bsVoxels.set (pt);
return value;
}, "~N,~N,~N,~N,~A");
c$.isInside = Clazz_defineMethod (c$, "isInside", 
function (voxelValue, max, isAbsolute) {
return ((max > 0 && (isAbsolute ? Math.abs (voxelValue) : voxelValue) >= max) || (max <= 0 && voxelValue <= max));
}, "~N,~N,~B");
Clazz_defineMethod (c$, "processOneCubical", 
function (insideMask, x, y, z, pt) {
var edgeMask = J.jvxl.calc.MarchingCubes.insideMaskTable[insideMask];
var isNaN = false;
for (var iEdge = 12; --iEdge >= 0; ) {
var xEdge = JU.TriangleData.Pwr2[iEdge];
if ((edgeMask & xEdge) == 0) continue;
var iPlane = this.edgeVertexPlanes[iEdge];
var iPt = (pt + this.linearOffsets[this.edgeVertexPointers[iEdge]]) % this.yzCount;
var iType = J.jvxl.calc.MarchingCubes.edgeTypeTable[iEdge];
var index = this.edgePointIndexes[iEdge] = this.isoPointIndexPlanes[iPlane][iPt][iType];
if (index != -2147483648) {
if (index == -1) isNaN = this.excludePartialCubes;
continue;
}var vertexA = JU.TriangleData.edgeVertexes[iEdge << 1];
var vertexB = JU.TriangleData.edgeVertexes[(iEdge << 1) + 1];
var valueA = this.vertexValues[vertexA];
var valueB = this.vertexValues[vertexB];
this.calcVertexPoint (x, y, z, vertexA, this.pointA);
this.edgeCount++;
var i = this.edgePointIndexes[iEdge] = this.isoPointIndexPlanes[iPlane][iPt][iType] = this.surfaceReader.getSurfacePointIndexAndFraction (this.cutoff, this.isCutoffAbsolute, x, y, z, JU.TriangleData.cubeVertexOffsets[vertexA], vertexA, vertexB, valueA, valueB, this.pointA, this.edgeVectors[iEdge], iType == this.contourType, this.fReturn);
this.addEdgeData (i < 0 ? NaN : this.fReturn[0]);
if (Float.isNaN (this.fReturn[0]) || i < 0) isNaN = this.excludePartialCubes;
}
return !isNaN;
}, "~N,~N,~N,~N,~N");
Clazz_defineMethod (c$, "addEdgeData", 
function (f) {
var ch = J.jvxl.data.JvxlCoder.jvxlFractionAsCharacter (f);
this.edgeData.appendC (ch);
}, "~N");
Clazz_defineMethod (c$, "calcVertexPoint", 
function (x, y, z, vertex, pt) {
this.volumeData.voxelPtToXYZ (x, y, z, this.pt0);
pt.add2 (this.pt0, this.voxelVertexVectors[vertex]);
}, "~N,~N,~N,~N,JU.P3");
Clazz_defineMethod (c$, "setLinearOffsets", 
function () {
this.linearOffsets[0] = 0;
this.linearOffsets[1] = this.yzCount;
this.linearOffsets[2] = this.yzCount + 1;
this.linearOffsets[3] = 1;
this.linearOffsets[4] = this.nZ;
this.linearOffsets[5] = this.yzCount + this.nZ;
this.linearOffsets[6] = this.yzCount + this.nZ + 1;
this.linearOffsets[7] = this.nZ + 1;
});
Clazz_defineMethod (c$, "getLinearOffset", 
function (x, y, z, offset) {
return x * this.yzCount + y * this.nZ + z + this.linearOffsets[offset];
}, "~N,~N,~N,~N");
Clazz_defineStatics (c$,
"MODE_CUBE", 1,
"MODE_JVXL", 2,
"MODE_PLANES", 3);
Clazz_defineStatics (c$,
"yzPlanePts", [0, 1, 1, 0, 0, 1, 1, 0],
"edgeVertexPointersLowToHigh", [1, 1, 2, 0, 5, 5, 6, 4, 0, 1, 2, 3],
"edgeVertexPointersHighToLow", [0, 1, 3, 0, 4, 5, 7, 4, 0, 1, 2, 3],
"edgeVertexPlanesLowToHigh", [1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0],
"edgeVertexPlanesHighToLow", [1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1]);
c$.cubeVertexVectors = c$.prototype.cubeVertexVectors = [JU.V3.new3 (0, 0, 0), JU.V3.new3 (1, 0, 0), JU.V3.new3 (1, 0, 1), JU.V3.new3 (0, 0, 1), JU.V3.new3 (0, 1, 0), JU.V3.new3 (1, 1, 0), JU.V3.new3 (1, 1, 1), JU.V3.new3 (0, 1, 1)];
Clazz_defineStatics (c$,
"edgeTypeTable", [0, 2, 0, 2, 0, 2, 0, 2, 1, 1, 1, 1],
"insideMaskTable", [0x0000, 0x0109, 0x0203, 0x030A, 0x0406, 0x050F, 0x0605, 0x070C, 0x080C, 0x0905, 0x0A0F, 0x0B06, 0x0C0A, 0x0D03, 0x0E09, 0x0F00, 0x0190, 0x0099, 0x0393, 0x029A, 0x0596, 0x049F, 0x0795, 0x069C, 0x099C, 0x0895, 0x0B9F, 0x0A96, 0x0D9A, 0x0C93, 0x0F99, 0x0E90, 0x0230, 0x0339, 0x0033, 0x013A, 0x0636, 0x073F, 0x0435, 0x053C, 0x0A3C, 0x0B35, 0x083F, 0x0936, 0x0E3A, 0x0F33, 0x0C39, 0x0D30, 0x03A0, 0x02A9, 0x01A3, 0x00AA, 0x07A6, 0x06AF, 0x05A5, 0x04AC, 0x0BAC, 0x0AA5, 0x09AF, 0x08A6, 0x0FAA, 0x0EA3, 0x0DA9, 0x0CA0, 0x0460, 0x0569, 0x0663, 0x076A, 0x0066, 0x016F, 0x0265, 0x036C, 0x0C6C, 0x0D65, 0x0E6F, 0x0F66, 0x086A, 0x0963, 0x0A69, 0x0B60, 0x05F0, 0x04F9, 0x07F3, 0x06FA, 0x01F6, 0x00FF, 0x03F5, 0x02FC, 0x0DFC, 0x0CF5, 0x0FFF, 0x0EF6, 0x09FA, 0x08F3, 0x0BF9, 0x0AF0, 0x0650, 0x0759, 0x0453, 0x055A, 0x0256, 0x035F, 0x0055, 0x015C, 0x0E5C, 0x0F55, 0x0C5F, 0x0D56, 0x0A5A, 0x0B53, 0x0859, 0x0950, 0x07C0, 0x06C9, 0x05C3, 0x04CA, 0x03C6, 0x02CF, 0x01C5, 0x00CC, 0x0FCC, 0x0EC5, 0x0DCF, 0x0CC6, 0x0BCA, 0x0AC3, 0x09C9, 0x08C0, 0x08C0, 0x09C9, 0x0AC3, 0x0BCA, 0x0CC6, 0x0DCF, 0x0EC5, 0x0FCC, 0x00CC, 0x01C5, 0x02CF, 0x03C6, 0x04CA, 0x05C3, 0x06C9, 0x07C0, 0x0950, 0x0859, 0x0B53, 0x0A5A, 0x0D56, 0x0C5F, 0x0F55, 0x0E5C, 0x015C, 0x0055, 0x035F, 0x0256, 0x055A, 0x0453, 0x0759, 0x0650, 0x0AF0, 0x0BF9, 0x08F3, 0x09FA, 0x0EF6, 0x0FFF, 0x0CF5, 0x0DFC, 0x02FC, 0x03F5, 0x00FF, 0x01F6, 0x06FA, 0x07F3, 0x04F9, 0x05F0, 0x0B60, 0x0A69, 0x0963, 0x086A, 0x0F66, 0x0E6F, 0x0D65, 0x0C6C, 0x036C, 0x0265, 0x016F, 0x0066, 0x076A, 0x0663, 0x0569, 0x0460, 0x0CA0, 0x0DA9, 0x0EA3, 0x0FAA, 0x08A6, 0x09AF, 0x0AA5, 0x0BAC, 0x04AC, 0x05A5, 0x06AF, 0x07A6, 0x00AA, 0x01A3, 0x02A9, 0x03A0, 0x0D30, 0x0C39, 0x0F33, 0x0E3A, 0x0936, 0x083F, 0x0B35, 0x0A3C, 0x053C, 0x0435, 0x073F, 0x0636, 0x013A, 0x0033, 0x0339, 0x0230, 0x0E90, 0x0F99, 0x0C93, 0x0D9A, 0x0A96, 0x0B9F, 0x0895, 0x099C, 0x069C, 0x0795, 0x049F, 0x0596, 0x029A, 0x0393, 0x0099, 0x0190, 0x0F00, 0x0E09, 0x0D03, 0x0C0A, 0x0B06, 0x0A0F, 0x0905, 0x080C, 0x070C, 0x0605, 0x050F, 0x0406, 0x030A, 0x0203, 0x0109, 0x0000]);
});
Clazz_declarePackage ("J.jvxl.calc");
Clazz_load (["JU.P3", "java.util.Hashtable"], "J.jvxl.calc.MarchingSquares", ["java.lang.Float", "JU.AU", "JU.Logger"], function () {
c$ = Clazz_decorateAsClass (function () {
this.surfaceReader = null;
this.volumeData = null;
this.nContourSegments = 0;
this.contourType = 0;
this.thisContour = 0;
this.valueMin = 0;
this.valueMax = 0;
this.pointA = null;
this.pointB = null;
this.contourFromZero = true;
this.contoursDiscrete = null;
this.contourVertexCount = 0;
this.contourVertexes = null;
if (!Clazz_isClassDefined ("J.jvxl.calc.MarchingSquares.ContourVertex")) {
J.jvxl.calc.MarchingSquares.$MarchingSquares$ContourVertex$ ();
}
this.contourPlaneMinimumValue = 0;
this.contourPlaneMaximumValue = 0;
this.contourValuesUsed = null;
this.ptTemp = null;
this.triangleCount = 0;
this.triangles = null;
this.htPts = null;
if (!Clazz_isClassDefined ("J.jvxl.calc.MarchingSquares.Triangle")) {
J.jvxl.calc.MarchingSquares.$MarchingSquares$Triangle$ ();
}
Clazz_instantialize (this, arguments);
}, J.jvxl.calc, "MarchingSquares");
Clazz_prepareFields (c$, function () {
this.pointA =  new JU.P3 ();
this.pointB =  new JU.P3 ();
this.contourVertexes =  new Array (1000);
this.ptTemp =  new JU.P3 ();
this.triangles =  new Array (1000);
this.htPts =  new java.util.Hashtable ();
});
Clazz_makeConstructor (c$, 
function (surfaceReader, volumeData, thePlane, contoursDiscrete, nContours, thisContour, contourFromZero) {
this.surfaceReader = surfaceReader;
this.volumeData = volumeData;
this.thisContour = thisContour;
this.contoursDiscrete = contoursDiscrete;
this.contourFromZero = contourFromZero;
if (contoursDiscrete == null) {
var i = 0;
this.nContourSegments = (nContours == 0 ? 9 : nContours) + i;
if (this.nContourSegments > 100) this.nContourSegments = 100;
} else {
nContours = contoursDiscrete.length;
this.nContourSegments = nContours;
this.contourFromZero = false;
}}, "J.jvxl.api.VertexDataServer,J.jvxl.data.VolumeData,JU.P4,~A,~N,~N,~B");
Clazz_defineMethod (c$, "getContourType", 
function () {
return this.contourType;
});
Clazz_defineMethod (c$, "setMinMax", 
function (valueMin, valueMax) {
this.valueMin = valueMin;
this.valueMax = valueMax;
}, "~N,~N");
Clazz_defineMethod (c$, "addContourVertex", 
function (vertexXYZ, value) {
if (this.contourVertexCount == this.contourVertexes.length) this.contourVertexes = JU.AU.doubleLength (this.contourVertexes);
var vPt = this.surfaceReader.addVertexCopy (vertexXYZ, value, -2, true);
this.contourVertexes[this.contourVertexCount++] = Clazz_innerTypeInstance (J.jvxl.calc.MarchingSquares.ContourVertex, this, null, vertexXYZ);
return vPt;
}, "JU.P3,~N");
Clazz_defineMethod (c$, "setContourData", 
function (i, value) {
this.contourVertexes[i].setValue (value);
}, "~N,~N");
Clazz_defineMethod (c$, "getContourValues", 
function () {
return this.contourValuesUsed;
});
Clazz_defineMethod (c$, "calcContourPoint", 
function (cutoff, valueA, valueB, pt) {
return this.volumeData.calculateFractionalPoint (cutoff, this.pointA, this.pointB, valueA, valueB, pt);
}, "~N,~N,~N,JU.P3");
Clazz_defineMethod (c$, "addTriangle", 
function (iA, iB, iC, check, check2) {
if (this.triangleCount == this.triangles.length) this.triangles = JU.AU.doubleLength (this.triangles);
this.triangles[this.triangleCount++] = Clazz_innerTypeInstance (J.jvxl.calc.MarchingSquares.Triangle, this, null, iA, iB, iC, check, check2);
return 0;
}, "~N,~N,~N,~N,~N");
Clazz_defineMethod (c$, "generateContourData", 
function (haveData, zeroOffset) {
JU.Logger.info ("generateContours: " + this.nContourSegments + " segments");
this.getVertexValues (haveData);
this.createContours (this.valueMin, this.valueMax, zeroOffset);
this.addAllTriangles ();
return this.contourVertexCount;
}, "~B,~N");
Clazz_defineMethod (c$, "getVertexValues", 
 function (haveData) {
this.contourPlaneMinimumValue = 3.4028235E38;
this.contourPlaneMaximumValue = -3.4028235E38;
for (var i = 0; i < this.contourVertexCount; i++) {
var c = this.contourVertexes[i];
var value;
if (haveData) {
value = c.value;
} else {
value = this.volumeData.lookupInterpolatedVoxelValue (c, false);
c.setValue (value);
}if (value < this.contourPlaneMinimumValue) this.contourPlaneMinimumValue = value;
if (value > this.contourPlaneMaximumValue) this.contourPlaneMaximumValue = value;
}
}, "~B");
Clazz_defineMethod (c$, "createContours", 
 function (min, max, zeroOffset) {
var diff = max - min;
this.contourValuesUsed =  Clazz_newFloatArray (this.nContourSegments, 0);
for (var i = this.triangleCount; --i >= 0; ) this.triangles[i].check = 0;

var minCutoff = -3.4028235E38;
var cutoff = minCutoff;
for (var i = 0; i < this.nContourSegments; i++) {
cutoff = (this.contoursDiscrete != null ? this.contoursDiscrete[i] : this.contourFromZero ? min + (i * 1 / this.nContourSegments) * diff : i == 0 ? -3.4028235E38 : i == this.nContourSegments - 1 ? 3.4028235E38 : min + ((i - 1) * 1 / (this.nContourSegments - 1)) * diff);
if (this.contoursDiscrete == null && Math.abs (cutoff) < zeroOffset) cutoff = (cutoff < 0 ? -zeroOffset : zeroOffset);
this.contourValuesUsed[i] = cutoff;
JU.Logger.info ("#contour " + (i + 1) + " " + cutoff);
this.htPts.clear ();
for (var ii = this.triangleCount; --ii >= 0; ) if (this.triangles[ii].isValid) this.checkContour (this.triangles[ii], i, cutoff);

if (this.thisContour > 0) {
if (i + 1 == this.thisContour) minCutoff = cutoff;
} else {
}}
if (this.contoursDiscrete != null) {
minCutoff = this.contoursDiscrete[0];
}this.valueMin = this.contourValuesUsed[0];
this.valueMax = (this.contourValuesUsed.length == 0 ? this.valueMin : this.contourValuesUsed[this.contourValuesUsed.length - 1]);
return true;
}, "~N,~N,~N");
Clazz_defineMethod (c$, "intercept", 
 function (t, i, value) {
var iA = t.pts[i];
var iB = t.pts[(i + 1) % 3];
if (iA == 2147483647 || iB == 2147483647) return -1;
var key = (iA < iB ? iA + "_" + iB : iB + "_" + iA);
if (this.htPts.containsKey (key)) return this.htPts.get (key).intValue ();
var valueA = this.contourVertexes[iA].value;
var valueB = this.contourVertexes[iB].value;
var iPt = -1;
if (valueA != valueB) {
var f = (value - valueA) / (valueB - valueA);
if (f >= 0 && f <= 1) {
this.pointA.setT (this.contourVertexes[iA]);
this.pointB.setT (this.contourVertexes[iB]);
value = this.calcContourPoint (value, valueA, valueB, this.ptTemp);
if (!Float.isNaN (value)) {
iPt = this.addContourVertex (this.ptTemp, value);
if (iPt < 0) return -1;
this.contourVertexes[iPt].setValue (value);
} else {
}}}this.htPts.put (key, Integer.$valueOf (iPt));
return iPt;
}, "J.jvxl.calc.MarchingSquares.Triangle,~N,~N");
Clazz_defineMethod (c$, "checkContour", 
 function (t, i, value) {
var ipt0 = this.intercept (t, 0, value);
var ipt1 = this.intercept (t, 1, value);
var ipt2 = this.intercept (t, 2, value);
var pts = t.pts;
var mode = 0;
if (ipt0 >= 0) {
mode += 1;
}if (ipt1 >= 0) {
mode += 2;
}if (ipt2 >= 0) {
mode += 4;
}switch (mode) {
case 3:
this.addTriangle (pts[0], ipt0, ipt1, 2 | (t.check & 1), i);
this.addTriangle (ipt0, pts[1], ipt1, 4 | (t.check & 3), i);
this.addTriangle (pts[0], ipt1, pts[2], (t.check & 6), i);
break;
case 5:
this.addTriangle (pts[0], ipt0, ipt2, 2 | (t.check & 5), i);
this.addTriangle (ipt0, pts[1], ipt2, 4 | (t.check & 1), i);
this.addTriangle (ipt2, pts[1], pts[2], (t.check & 6), i);
break;
case 6:
this.addTriangle (pts[0], pts[1], ipt2, (t.check & 5), i);
this.addTriangle (ipt2, pts[1], ipt1, 4 | (t.check & 2), i);
this.addTriangle (ipt2, ipt1, pts[2], 1 | (t.check & 6), i);
break;
default:
return;
}
t.isValid = false;
}, "J.jvxl.calc.MarchingSquares.Triangle,~N,~N");
Clazz_defineMethod (c$, "getMinMax", 
function () {
return [this.valueMin, this.valueMax];
});
Clazz_defineMethod (c$, "addAllTriangles", 
 function () {
for (var i = 0; i < this.triangleCount; i++) if (this.triangles[i].isValid) {
var t = this.triangles[i];
this.surfaceReader.addTriangleCheck (t.pts[0], t.pts[1], t.pts[2], t.check, t.contourIndex, false, -1);
}
});
c$.$MarchingSquares$ContourVertex$ = function () {
Clazz_pu$h(self.c$);
c$ = Clazz_decorateAsClass (function () {
Clazz_prepareCallback (this, arguments);
this.value = 0;
Clazz_instantialize (this, arguments);
}, J.jvxl.calc.MarchingSquares, "ContourVertex", JU.P3);
Clazz_makeConstructor (c$, 
function (a) {
Clazz_superConstructor (this, J.jvxl.calc.MarchingSquares.ContourVertex, []);
this.setT (a);
}, "JU.P3");
Clazz_defineMethod (c$, "setValue", 
function (a) {
this.value = a;
}, "~N");
c$ = Clazz_p0p ();
};
c$.$MarchingSquares$Triangle$ = function () {
Clazz_pu$h(self.c$);
c$ = Clazz_decorateAsClass (function () {
Clazz_prepareCallback (this, arguments);
this.pts = null;
this.check = 0;
this.isValid = true;
this.contourIndex = 0;
Clazz_instantialize (this, arguments);
}, J.jvxl.calc.MarchingSquares, "Triangle");
Clazz_makeConstructor (c$, 
function (a, b, c, d, e) {
this.pts = [a, b, c];
this.check = d;
this.contourIndex = e;
}, "~N,~N,~N,~N,~N");
c$ = Clazz_p0p ();
};
Clazz_defineStatics (c$,
"CONTOUR_POINT", -1,
"VERTEX_POINT", -2,
"EDGE_POINT", -3,
"nContourMax", 100,
"defaultContourCount", 9);
});
Clazz_declarePackage ("J.shapesurface");
Clazz_load (["J.shape.Mesh", "J.jvxl.data.JvxlData"], "J.shapesurface.IsosurfaceMesh", ["java.lang.Character", "$.Float", "java.util.Hashtable", "JU.AU", "$.BS", "$.CU", "$.Lst", "$.M4", "$.P3", "$.P4", "$.PT", "$.SB", "$.V3", "J.api.Interface", "J.jvxl.data.JvxlCoder", "JS.T", "JU.BoxInfo", "$.C", "$.ColorEncoder", "$.Escape", "$.Logger", "$.Measure", "JV.Viewer"], function () {
c$ = Clazz_decorateAsClass (function () {
this.jvxlData = null;
this.vertexIncrement = 1;
this.firstRealVertex = -1;
this.dataType = 0;
this.hasGridPoints = false;
this.calculatedArea = null;
this.calculatedVolume = null;
this.info = null;
this.assocGridPointMap = null;
this.assocGridPointNormals = null;
this.mergeAssociatedNormalCount = 0;
this.centers = null;
this.contourValues = null;
this.contourColixes = null;
this.colorEncoder = null;
this.bsVdw = null;
this.colorPhased = false;
Clazz_instantialize (this, arguments);
}, J.shapesurface, "IsosurfaceMesh", J.shape.Mesh);
Clazz_prepareFields (c$, function () {
this.jvxlData =  new J.jvxl.data.JvxlData ();
});
Clazz_makeConstructor (c$, 
function (thisID, colix, index) {
Clazz_superConstructor (this, J.shapesurface.IsosurfaceMesh, []);
this.mesh1 (thisID, colix, index);
this.checkByteCount = 2;
this.jvxlData.version = JV.Viewer.getJmolVersion ();
}, "~S,~N,~N");
Clazz_defineMethod (c$, "clearType", 
function (meshType, iAddGridPoints) {
this.clear (meshType);
this.jvxlData.clear ();
this.assocGridPointMap = null;
this.assocGridPointNormals = null;
this.bsVdw = null;
this.calculatedVolume = null;
this.calculatedArea = null;
this.centers = null;
this.colorEncoder = null;
this.colorPhased = false;
this.firstRealVertex = -1;
this.hasGridPoints = iAddGridPoints;
this.isColorSolid = true;
this.mergeAssociatedNormalCount = 0;
this.nSets = 0;
this.pcs = null;
this.showPoints = iAddGridPoints;
this.surfaceSet = null;
this.vcs = null;
this.vertexColorMap = null;
this.vertexIncrement = 1;
this.vertexSets = null;
this.vvs = null;
}, "~S,~B");
Clazz_defineMethod (c$, "allocVertexColixes", 
function () {
if (this.vcs == null) {
this.vcs =  Clazz_newShortArray (this.vc, 0);
for (var i = this.vc; --i >= 0; ) this.vcs[i] = this.colix;

}this.isColorSolid = false;
});
Clazz_defineMethod (c$, "addVertexCopy", 
function (vertex, value, assocVertex, associateNormals, asCopy) {
var vPt = this.addVCVal (vertex, value, asCopy);
switch (assocVertex) {
case -1:
if (this.firstRealVertex < 0) this.firstRealVertex = vPt;
break;
case -2:
this.hasGridPoints = true;
break;
case -3:
this.vertexIncrement = 3;
break;
default:
if (this.firstRealVertex < 0) this.firstRealVertex = vPt;
if (associateNormals) {
if (this.assocGridPointMap == null) this.assocGridPointMap =  new java.util.Hashtable ();
this.assocGridPointMap.put (Integer.$valueOf (vPt), Integer.$valueOf (assocVertex + this.mergeAssociatedNormalCount));
}}
return vPt;
}, "JU.T3,~N,~N,~B,~B");
Clazz_overrideMethod (c$, "setTranslucent", 
function (isTranslucent, iLevel) {
this.colix = JU.C.getColixTranslucent3 (this.colix, isTranslucent, iLevel);
if (this.vcs != null) for (var i = this.vc; --i >= 0; ) this.vcs[i] = JU.C.getColixTranslucent3 (this.vcs[i], isTranslucent, iLevel);

}, "~B,~N");
Clazz_defineMethod (c$, "setMerged", 
function (TF) {
this.isMerged = TF;
this.mergePolygonCount0 = (TF ? this.pc : 0);
this.mergeVertexCount0 = (TF ? this.vc : 0);
if (TF) {
this.mergeAssociatedNormalCount += this.jvxlData.nPointsX * this.jvxlData.nPointsY * this.jvxlData.nPointsZ;
this.assocGridPointNormals = null;
}}, "~B");
Clazz_overrideMethod (c$, "sumVertexNormals", 
function (vertices, vectorSums) {
this.sumVertexNormals2 (vertices, vectorSums);
if (this.assocGridPointMap != null && vectorSums.length > 0 && !this.isMerged) {
if (this.assocGridPointNormals == null) this.assocGridPointNormals =  new java.util.Hashtable ();
for (var entry, $entry = this.assocGridPointMap.entrySet ().iterator (); $entry.hasNext () && ((entry = $entry.next ()) || true);) {
var gridPoint = entry.getValue ();
if (!this.assocGridPointNormals.containsKey (gridPoint)) this.assocGridPointNormals.put (gridPoint, JU.V3.new3 (0, 0, 0));
this.assocGridPointNormals.get (gridPoint).add (vectorSums[entry.getKey ().intValue ()]);
}
for (var entry, $entry = this.assocGridPointMap.entrySet ().iterator (); $entry.hasNext () && ((entry = $entry.next ()) || true);) vectorSums[entry.getKey ().intValue ()] = this.assocGridPointNormals.get (entry.getValue ());

}}, "~A,~A");
Clazz_defineMethod (c$, "getCenters", 
function () {
if (this.centers != null) return this.centers;
this.centers =  new Array (this.pc);
for (var i = 0; i < this.pc; i++) {
var p = this.pis[i];
if (p == null) continue;
var pt = this.centers[i] = JU.P3.newP (this.vs[p[0]]);
pt.add (this.vs[p[1]]);
pt.add (this.vs[p[2]]);
pt.scale (0.33333334);
}
return this.centers;
});
Clazz_defineMethod (c$, "getFacePlane", 
function (i, vNorm) {
var plane =  new JU.P4 ();
JU.Measure.getPlaneThroughPoints (this.vs[this.pis[i][0]], this.vs[this.pis[i][1]], this.vs[this.pis[i][2]], vNorm, this.vAB, this.vAC, plane);
return plane;
}, "~N,JU.V3");
Clazz_defineMethod (c$, "getContours", 
function () {
var n = this.jvxlData.nContours;
if (n == 0 || this.pis == null) return null;
this.havePlanarContours = (this.jvxlData.jvxlPlane != null);
if (this.havePlanarContours) return null;
if (n < 0) n = -1 - n;
var vContours = this.jvxlData.vContours;
if (vContours != null) {
for (var i = 0; i < n; i++) {
if (vContours[i].size () > 6) return this.jvxlData.vContours;
J.jvxl.data.JvxlCoder.set3dContourVector (vContours[i], this.pis, this.vs);
}
return this.jvxlData.vContours;
}vContours =  new Array (n);
for (var i = 0; i < n; i++) {
vContours[i] =  new JU.Lst ();
}
if (this.jvxlData.contourValuesUsed == null) {
var dv = (this.jvxlData.valueMappedToBlue - this.jvxlData.valueMappedToRed) / (n + 1);
for (var i = 0; i < n; i++) {
var value = this.jvxlData.valueMappedToRed + (i + 1) * dv;
this.get3dContour (vContours[i], value, this.jvxlData.contourColixes[i]);
}
JU.Logger.info (n + " contour lines; separation = " + dv);
} else {
for (var i = 0; i < n; i++) {
var value = this.jvxlData.contourValuesUsed[i];
this.get3dContour (vContours[i], value, this.jvxlData.contourColixes[i]);
}
}this.jvxlData.contourColixes =  Clazz_newShortArray (n, 0);
this.jvxlData.contourValues =  Clazz_newFloatArray (n, 0);
for (var i = 0; i < n; i++) {
this.jvxlData.contourValues[i] = (vContours[i].get (2)).floatValue ();
this.jvxlData.contourColixes[i] = (vContours[i].get (3))[0];
}
return this.jvxlData.vContours = vContours;
});
Clazz_defineMethod (c$, "get3dContour", 
 function (v, value, colix) {
var bsContour = JU.BS.newN (this.pc);
var fData =  new JU.SB ();
var color = JU.C.getArgb (colix);
J.shapesurface.IsosurfaceMesh.setContourVector (v, this.pc, bsContour, value, colix, color, fData);
for (var i = 0; i < this.pc; i++) if (this.setABC (i)) J.shapesurface.IsosurfaceMesh.addContourPoints (v, bsContour, i, fData, this.vs, this.vvs, this.iA, this.iB, this.iC, value);

}, "JU.Lst,~N,~N");
c$.setContourVector = Clazz_defineMethod (c$, "setContourVector", 
function (v, nPolygons, bsContour, value, colix, color, fData) {
v.add (0, Integer.$valueOf (nPolygons));
v.add (1, bsContour);
v.add (2, Float.$valueOf (value));
v.add (3, [colix]);
v.add (4, [color]);
v.add (5, fData);
}, "JU.Lst,~N,JU.BS,~N,~N,~N,JU.SB");
c$.addContourPoints = Clazz_defineMethod (c$, "addContourPoints", 
function (v, bsContour, i, fData, vertices, vertexValues, iA, iB, iC, value) {
var pt1 = null;
var pt2 = null;
var type = 0;
var f1 = J.shapesurface.IsosurfaceMesh.checkPt (vertexValues, iA, iB, value);
if (!Float.isNaN (f1)) {
pt1 = J.shapesurface.IsosurfaceMesh.getContourPoint (vertices, iA, iB, f1);
type |= 1;
}var f2 = (f1 == 1 ? NaN : J.shapesurface.IsosurfaceMesh.checkPt (vertexValues, iB, iC, value));
if (!Float.isNaN (f2)) {
pt2 = J.shapesurface.IsosurfaceMesh.getContourPoint (vertices, iB, iC, f2);
if (type == 0) {
pt1 = pt2;
f1 = f2;
}type |= 2;
}switch (type) {
case 0:
return;
case 1:
if (f1 == 0) return;
case 2:
f2 = (f2 == 1 ? NaN : J.shapesurface.IsosurfaceMesh.checkPt (vertexValues, iC, iA, value));
if (!Float.isNaN (f2)) {
pt2 = J.shapesurface.IsosurfaceMesh.getContourPoint (vertices, iC, iA, f2);
type |= 4;
}break;
}
switch (type) {
case 3:
case 5:
case 6:
break;
default:
return;
}
bsContour.set (i);
J.jvxl.data.JvxlCoder.appendContourTriangleIntersection (type, f1, f2, fData);
v.addLast (pt1);
v.addLast (pt2);
}, "JU.Lst,JU.BS,~N,JU.SB,~A,~A,~N,~N,~N,~N");
c$.checkPt = Clazz_defineMethod (c$, "checkPt", 
 function (vertexValues, i, j, v) {
var v1;
var v2;
return (v == (v1 = vertexValues[i]) ? 0 : v == (v2 = vertexValues[j]) ? 1 : (v1 < v) == (v < v2) ? (v - v1) / (v2 - v1) : NaN);
}, "~A,~N,~N,~N");
c$.getContourPoint = Clazz_defineMethod (c$, "getContourPoint", 
 function (vertices, i, j, f) {
var pt =  new JU.P3 ();
pt.sub2 (vertices[j], vertices[i]);
pt.scaleAdd2 (f, pt, vertices[i]);
return pt;
}, "~A,~N,~N,~N");
Clazz_defineMethod (c$, "setDiscreteColixes", 
function (values, colixes) {
if (values != null) this.jvxlData.contourValues = values;
if (values == null || values.length == 0) values = this.jvxlData.contourValues = this.jvxlData.contourValuesUsed;
if (colixes == null && this.jvxlData.contourColixes != null) {
colixes = this.jvxlData.contourColixes;
} else {
this.jvxlData.contourColixes = colixes;
this.jvxlData.contourColors = JU.C.getHexCodes (colixes);
}if (this.vs == null || this.vvs == null || values == null) return;
var n = values.length;
var vMax = values[n - 1];
this.colorCommand = null;
var haveColixes = (colixes != null && colixes.length > 0);
this.isColorSolid = (haveColixes && this.jvxlData.jvxlPlane != null);
if (this.jvxlData.vContours != null) {
if (haveColixes) for (var i = 0; i < this.jvxlData.vContours.length; i++) {
var colix = colixes[i % colixes.length];
(this.jvxlData.vContours[i].get (3))[0] = colix;
(this.jvxlData.vContours[i].get (4))[0] = JU.C.getArgb (colix);
}
return;
}var defaultColix = 0;
this.pcs =  Clazz_newShortArray (this.pc, 0);
for (var i = 0; i < this.pc; i++) {
var p = this.pis[i];
if (p == null) continue;
this.pcs[i] = defaultColix;
var v = (this.vvs[p[0]] + this.vvs[p[1]] + this.vvs[p[2]]) / 3;
for (var j = n; --j >= 0; ) {
if (v >= values[j] && v < vMax) {
this.pcs[i] = (haveColixes ? colixes[j % colixes.length] : 0);
break;
}}
}
}, "~A,~A");
Clazz_defineMethod (c$, "getContourList", 
function (vwr) {
var ht =  new java.util.Hashtable ();
ht.put ("values", (this.jvxlData.contourValuesUsed == null ? this.jvxlData.contourValues : this.jvxlData.contourValuesUsed));
var colors =  new JU.Lst ();
if (this.jvxlData.contourColixes != null) {
for (var i = 0; i < this.jvxlData.contourColixes.length; i++) {
colors.addLast (JU.CU.colorPtFromInt (JU.C.getArgb (this.jvxlData.contourColixes[i]), null));
}
ht.put ("colors", colors);
}return ht;
}, "JV.Viewer");
Clazz_defineMethod (c$, "deleteContours", 
function () {
this.jvxlData.contourValuesUsed = null;
this.jvxlData.contourValues = null;
this.jvxlData.contourColixes = null;
this.jvxlData.vContours = null;
});
Clazz_defineMethod (c$, "setVertexColorMap", 
function () {
this.vertexColorMap =  new java.util.Hashtable ();
var lastColix = -999;
var bs = null;
for (var i = this.vc; --i >= 0; ) {
var c = this.vcs[i];
if (c != lastColix) {
var color = JU.C.getHexCode (lastColix = c);
bs = this.vertexColorMap.get (color);
if (bs == null) this.vertexColorMap.put (color, bs =  new JU.BS ());
}bs.set (i);
}
});
Clazz_defineMethod (c$, "setVertexColixesForAtoms", 
function (vwr, colixes, atomMap, bs) {
this.jvxlData.vertexDataOnly = true;
this.jvxlData.vertexColors =  Clazz_newIntArray (this.vc, 0);
this.jvxlData.nVertexColors = this.vc;
var atoms = vwr.ms.at;
for (var i = this.mergeVertexCount0; i < this.vc; i++) {
var iAtom = this.vertexSource[i];
if (iAtom < 0 || !bs.get (iAtom)) continue;
this.jvxlData.vertexColors[i] = vwr.getColorArgbOrGray (this.vcs[i] = JU.C.copyColixTranslucency (this.colix, atoms[iAtom].getColix ()));
var colix = (colixes == null ? 0 : colixes[atomMap[iAtom]]);
if (colix == 0) colix = atoms[iAtom].getColix ();
this.vcs[i] = JU.C.copyColixTranslucency (this.colix, colix);
}
}, "JV.Viewer,~A,~A,JU.BS");
Clazz_defineMethod (c$, "colorVertices", 
function (colix, bs, isAtoms) {
if (this.vertexSource == null) return;
colix = JU.C.copyColixTranslucency (this.colix, colix);
var bsVertices = (isAtoms ?  new JU.BS () : bs);
this.checkAllocColixes ();
if (isAtoms) for (var i = 0; i < this.vc; i++) {
var pt = this.vertexSource[i];
if (pt >= 0 && bs.get (pt)) {
this.vcs[i] = colix;
if (bsVertices != null) bsVertices.set (i);
}}
 else for (var i = 0; i < this.vc; i++) if (bsVertices.get (i)) this.vcs[i] = colix;

if (!isAtoms) {
return;
}var color = JU.C.getHexCode (colix);
if (this.vertexColorMap == null) this.vertexColorMap =  new java.util.Hashtable ();
J.shapesurface.IsosurfaceMesh.addColorToMap (this.vertexColorMap, color, bs);
}, "~N,JU.BS,~B");
Clazz_defineMethod (c$, "checkAllocColixes", 
function () {
if (this.vcs == null || this.vertexColorMap == null && this.isColorSolid) this.allocVertexColixes ();
this.isColorSolid = false;
});
c$.addColorToMap = Clazz_defineMethod (c$, "addColorToMap", 
 function (colorMap, color, bs) {
var bsMap = null;
for (var entry, $entry = colorMap.entrySet ().iterator (); $entry.hasNext () && ((entry = $entry.next ()) || true);) if (entry.getKey () === color) {
bsMap = entry.getValue ();
bsMap.or (bs);
} else {
entry.getValue ().andNot (bs);
}
if (bsMap == null) colorMap.put (color, bs);
}, "java.util.Map,~S,JU.BS");
Clazz_defineMethod (c$, "setJvxlColorMap", 
function (isAll) {
this.jvxlData.diameter = this.diameter;
this.jvxlData.color = JU.C.getHexCode (this.colix);
this.jvxlData.meshColor = (this.meshColix == 0 ? null : JU.C.getHexCode (this.meshColix));
this.jvxlData.translucency = JU.C.getColixTranslucencyFractional (this.colix);
this.jvxlData.rendering = this.getRendering ().substring (1);
this.jvxlData.colorScheme = (this.colorEncoder == null ? null : this.colorEncoder.getColorScheme ());
if (this.jvxlData.vertexColors == null) this.jvxlData.nVertexColors = (this.vertexColorMap == null ? 0 : this.vertexColorMap.size ());
if (this.vertexColorMap == null || this.vertexSource == null || !isAll) return;
if (this.jvxlData.vertexColorMap == null) this.jvxlData.vertexColorMap =  new java.util.Hashtable ();
for (var entry, $entry = this.vertexColorMap.entrySet ().iterator (); $entry.hasNext () && ((entry = $entry.next ()) || true);) {
var bsMap = entry.getValue ();
if (bsMap.isEmpty ()) continue;
var color = entry.getKey ();
var bs =  new JU.BS ();
for (var i = 0; i < this.vc; i++) if (bsMap.get (this.vertexSource[i])) bs.set (i);

J.shapesurface.IsosurfaceMesh.addColorToMap (this.jvxlData.vertexColorMap, color, bs);
}
this.jvxlData.nVertexColors = this.jvxlData.vertexColorMap.size ();
if (this.jvxlData.vertexColorMap.size () == 0) this.jvxlData.vertexColorMap = null;
}, "~B");
Clazz_defineMethod (c$, "setColorCommand", 
function () {
if (this.colorEncoder == null) return;
this.colorCommand = this.colorEncoder.getColorScheme ();
if (this.colorCommand.equals ("inherit")) {
this.colorCommand = "#inherit;";
return;
}if (this.colorCommand == null) return;
this.colorCommand = "color $" + (Character.isLetter (this.thisID.charAt (0)) && this.thisID.indexOf (" ") < 0 ? this.thisID : "\"" + this.thisID + "\"") + " \"" + this.colorCommand + "\" range " + (this.jvxlData.isColorReversed ? this.jvxlData.valueMappedToBlue + " " + this.jvxlData.valueMappedToRed : this.jvxlData.valueMappedToRed + " " + this.jvxlData.valueMappedToBlue);
});
Clazz_defineMethod (c$, "setColorsFromJvxlData", 
function (colorRgb) {
this.diameter = this.jvxlData.diameter;
if (colorRgb == -1) {
} else if (colorRgb != -2147483648 && colorRgb != 2147483647) {
this.colix = JU.C.getColix (colorRgb);
} else if (this.jvxlData.color != null) {
this.colix = JU.C.getColixS (this.jvxlData.color);
}if (this.colix == 0) this.colix = 5;
this.colix = JU.C.getColixTranslucent3 (this.colix, this.jvxlData.translucency != 0, this.jvxlData.translucency);
if (this.jvxlData.meshColor != null) this.meshColix = JU.C.getColixS (this.jvxlData.meshColor);
this.setJvxlDataRendering ();
this.isColorSolid = !this.jvxlData.isBicolorMap && this.jvxlData.vertexColors == null && this.jvxlData.vertexColorMap == null;
if (this.colorEncoder != null) {
if (this.jvxlData.vertexColorMap == null) {
if (this.jvxlData.colorScheme != null) {
var colorScheme = this.jvxlData.colorScheme;
var isTranslucent = colorScheme.startsWith ("translucent ");
if (isTranslucent) colorScheme = colorScheme.substring (12);
this.colorEncoder.setColorScheme (colorScheme, isTranslucent);
this.remapColors (null, null, NaN);
}} else {
if (this.jvxlData.baseColor != null) {
for (var i = this.vc; --i >= 0; ) this.vcs[i] = this.colix;

}for (var entry, $entry = this.jvxlData.vertexColorMap.entrySet ().iterator (); $entry.hasNext () && ((entry = $entry.next ()) || true);) {
var bsMap = entry.getValue ();
var colix = JU.C.copyColixTranslucency (this.colix, JU.C.getColixS (entry.getKey ()));
for (var i = bsMap.nextSetBit (0); i >= 0; i = bsMap.nextSetBit (i + 1)) this.vcs[i] = colix;

}
}}}, "~N");
Clazz_defineMethod (c$, "setJvxlDataRendering", 
function () {
if (this.jvxlData.rendering != null) {
var tokens = JU.PT.getTokens (this.jvxlData.rendering);
for (var i = 0; i < tokens.length; i++) this.setTokenProperty (JS.T.getTokFromName (tokens[i]), true);

}});
Clazz_defineMethod (c$, "remapColors", 
function (vwr, ce, translucentLevel) {
if (ce == null) ce = this.colorEncoder;
if (ce == null) ce = this.colorEncoder =  new JU.ColorEncoder (null);
this.colorEncoder = ce;
this.setColorCommand ();
if (Float.isNaN (translucentLevel)) {
translucentLevel = JU.C.getColixTranslucencyLevel (this.colix);
} else {
this.colix = JU.C.getColixTranslucent3 (this.colix, true, translucentLevel);
}var min = ce.lo;
var max = ce.hi;
var inherit = (this.vertexSource != null && ce.currentPalette == 15);
this.vertexColorMap = null;
this.pcs = null;
this.jvxlData.baseColor = null;
this.jvxlData.vertexCount = this.vc;
if (this.vvs == null || this.jvxlData.vertexCount == 0) return;
if (this.vcs == null || this.vcs.length != this.vc) this.allocVertexColixes ();
if (inherit) {
this.jvxlData.vertexDataOnly = true;
this.jvxlData.vertexColors =  Clazz_newIntArray (this.vc, 0);
this.jvxlData.nVertexColors = this.vc;
var atoms = vwr.ms.at;
for (var i = this.mergeVertexCount0; i < this.vc; i++) {
var pt = this.vertexSource[i];
if (pt >= 0 && pt < atoms.length) this.jvxlData.vertexColors[i] = vwr.getColorArgbOrGray (this.vcs[i] = JU.C.copyColixTranslucency (this.colix, atoms[pt].getColix ()));
}
return;
}this.jvxlData.vertexColors = null;
this.jvxlData.vertexColorMap = null;
if (this.jvxlData.isBicolorMap) {
for (var i = this.mergeVertexCount0; i < this.vc; i++) this.vcs[i] = JU.C.copyColixTranslucency (this.colix, this.vvs[i] < 0 ? this.jvxlData.minColorIndex : this.jvxlData.maxColorIndex);

return;
}this.jvxlData.isColorReversed = ce.isReversed;
if (max != 3.4028235E38) {
this.jvxlData.valueMappedToRed = min;
this.jvxlData.valueMappedToBlue = max;
}ce.setRange (this.jvxlData.valueMappedToRed, this.jvxlData.valueMappedToBlue, this.jvxlData.isColorReversed);
var isTranslucent = JU.C.isColixTranslucent (this.colix);
if (ce.isTranslucent) {
if (!isTranslucent) this.colix = JU.C.getColixTranslucent3 (this.colix, true, 0.5);
isTranslucent = false;
}for (var i = this.vc; --i >= this.mergeVertexCount0; ) this.vcs[i] = ce.getColorIndex (this.vvs[i]);

this.setTranslucent (isTranslucent, translucentLevel);
this.colorEncoder = ce;
var contours = this.getContours ();
if (contours != null) {
for (var i = contours.length; --i >= 0; ) {
var value = (contours[i].get (2)).floatValue ();
var colix = (contours[i].get (3));
colix[0] = ce.getColorIndex (value);
var color = (contours[i].get (4));
color[0] = JU.C.getArgb (colix[0]);
}
}if (this.contourValues != null) {
this.contourColixes =  Clazz_newShortArray (this.contourValues.length, 0);
for (var i = 0; i < this.contourValues.length; i++) this.contourColixes[i] = ce.getColorIndex (this.contourValues[i]);

this.setDiscreteColixes (null, null);
}this.jvxlData.isJvxlPrecisionColor = true;
J.jvxl.data.JvxlCoder.jvxlCreateColorData (this.jvxlData, this.vvs);
this.setColorCommand ();
this.isColorSolid = false;
}, "JV.Viewer,JU.ColorEncoder,~N");
Clazz_defineMethod (c$, "reinitializeLightingAndColor", 
function (vwr) {
this.initialize (this.lighting, null, null);
if (this.colorEncoder != null || this.jvxlData.isBicolorMap) {
this.vcs = null;
this.remapColors (vwr, null, NaN);
}}, "JV.Viewer");
Clazz_overrideMethod (c$, "getBoundingBox", 
function () {
return this.jvxlData.boundingBox;
});
Clazz_defineMethod (c$, "resetBoundingBox", 
 function () {
var bi =  new JU.BoxInfo ();
if (this.pc == 0) for (var i = this.vc; --i >= 0; ) {
bi.addBoundBoxPoint (this.vs[i]);
}
 else {
var bsDone =  new JU.BS ();
for (var i = this.pc; --i >= 0; ) {
if (!this.setABC (i)) continue;
if (!bsDone.get (this.iA)) {
bi.addBoundBoxPoint (this.vs[this.iA]);
bsDone.set (this.iA);
}if (!bsDone.get (this.iB)) {
bi.addBoundBoxPoint (this.vs[this.iB]);
bsDone.set (this.iB);
}if (!bsDone.get (this.iC)) {
bi.addBoundBoxPoint (this.vs[this.iC]);
bsDone.set (this.iC);
}}
}this.jvxlData.boundingBox = bi.getBoundBoxPoints (false);
});
Clazz_defineMethod (c$, "merge", 
function (m) {
var nV = this.vc + (m == null ? 0 : m.vc);
if (this.pis == null) this.pis =  Clazz_newIntArray (0, 0, 0);
if (m != null && m.pis == null) m.pis =  Clazz_newIntArray (0, 0, 0);
var nP = (this.bsSlabDisplay == null || this.pc == 0 ? this.pc : this.bsSlabDisplay.cardinality ()) + (m == null || m.pc == 0 ? 0 : m.bsSlabDisplay == null ? m.pc : m.bsSlabDisplay.cardinality ());
if (this.vs == null) this.vs =  new Array (0);
this.vs = JU.AU.ensureLength (this.vs, nV);
this.vvs = JU.AU.ensureLengthA (this.vvs, nV);
var haveSources = (this.vertexSource != null && (m == null || m.vertexSource != null));
this.vertexSource = JU.AU.ensureLengthI (this.vertexSource, nV);
var newPolygons = JU.AU.newInt2 (nP);
var ipt = J.shapesurface.IsosurfaceMesh.mergePolygons (this, 0, 0, newPolygons);
if (m != null) {
ipt = J.shapesurface.IsosurfaceMesh.mergePolygons (m, ipt, this.vc, newPolygons);
for (var i = 0; i < m.vc; i++, this.vc++) {
this.vs[this.vc] = m.vs[i];
this.vvs[this.vc] = m.vvs[i];
if (haveSources) this.vertexSource[this.vc] = m.vertexSource[i];
}
}this.pc = this.polygonCount0 = nP;
this.vc = this.vertexCount0 = nV;
if (nP > 0) this.resetSlab ();
this.pis = newPolygons;
}, "J.jvxl.data.MeshData");
c$.mergePolygons = Clazz_defineMethod (c$, "mergePolygons", 
 function (m, ipt, vertexCount, newPolygons) {
var p;
for (var i = 0; i < m.pc; i++) {
if ((p = m.pis[i]) == null || m.bsSlabDisplay != null && !m.bsSlabDisplay.get (i)) continue;
newPolygons[ipt++] = m.pis[i];
if (vertexCount > 0) for (var j = 0; j < 3; j++) p[j] += vertexCount;

}
return ipt;
}, "JU.MeshSurface,~N,~N,~A");
Clazz_overrideMethod (c$, "getUnitCell", 
function () {
return (this.spanningVectors == null ? null : J.api.Interface.getSymmetry ().getUnitCell (this.spanningVectors, true, null));
});
Clazz_overrideMethod (c$, "slabBrillouin", 
function (unitCellPoints) {
var vectors = (unitCellPoints == null ? this.spanningVectors : unitCellPoints);
if (vectors == null) return;
var pts =  new Array (27);
pts[0] = JU.P3.newP (vectors[0]);
var pt = 0;
for (var i = -1; i <= 1; i++) for (var j = -1; j <= 1; j++) for (var k = -1; k <= 1; k++) if (i != 0 || j != 0 || k != 0) {
pts[++pt] = JU.P3.newP (pts[0]);
pts[pt].scaleAdd2 (i, vectors[1], pts[pt]);
pts[pt].scaleAdd2 (j, vectors[2], pts[pt]);
pts[pt].scaleAdd2 (k, vectors[3], pts[pt]);
}


System.out.println ("draw line1 {0 0 0} color red" + JU.Escape.eP (this.spanningVectors[1]));
System.out.println ("draw line2 {0 0 0} color green" + JU.Escape.eP (this.spanningVectors[2]));
System.out.println ("draw line3 {0 0 0} color blue" + JU.Escape.eP (this.spanningVectors[3]));
var ptTemp =  new JU.P3 ();
var planeGammaK =  new JU.P4 ();
var vGammaToKPoint =  new JU.V3 ();
var vTemp =  new JU.V3 ();
var bsMoved =  new JU.BS ();
var mapEdge =  new java.util.Hashtable ();
this.bsSlabGhost =  new JU.BS ();
for (var i = 1; i < 27; i++) {
vGammaToKPoint.setT (pts[i]);
JU.Measure.getBisectingPlane (pts[0], vGammaToKPoint, ptTemp, vTemp, planeGammaK);
this.getIntersection (1, planeGammaK, null, null, null, null, null, false, false, 135266319, true);
bsMoved.clearAll ();
mapEdge.clear ();
for (var j = this.bsSlabGhost.nextSetBit (0); j >= 0; j = this.bsSlabGhost.nextSetBit (j + 1)) {
if (!this.setABC (j)) continue;
var p = JU.AU.arrayCopyRangeI (this.pis[j], 0, -1);
for (var k = 0; k < 3; k++) {
var pk = p[k];
p[k] = this.addIntersectionVertex (this.vs[pk], this.vvs[pk], this.vertexSource == null ? 0 : this.vertexSource[pk], this.vertexSets == null ? 0 : this.vertexSets[pk], mapEdge, 0, pk);
if (pk != p[k] && bsMoved.get (pk)) bsMoved.set (p[k]);
}
this.addPolygonC (p, 0, this.bsSlabDisplay);
for (var k = 0; k < 3; k++) if (!bsMoved.get (p[k])) {
bsMoved.set (p[k]);
this.vs[p[k]].sub (vGammaToKPoint);
}
}
if (this.bsSlabGhost.nextSetBit (0) >= 0) {
this.bsSlabGhost.clearAll ();
i = 0;
}}
this.bsSlabGhost = null;
this.resetBoundingBox ();
}, "~A");
Clazz_overrideMethod (c$, "getMinDistance2ForVertexGrouping", 
function () {
if (this.jvxlData.boundingBox != null && this.jvxlData.boundingBox[0] != null) {
var d2 = this.jvxlData.boundingBox[1].distanceSquared (this.jvxlData.boundingBox[0]);
if (d2 < 5) return 1e-10;
}return 1e-8;
});
Clazz_overrideMethod (c$, "getVisibleVertexBitSet", 
function () {
var bs = this.getVisibleVBS ();
if (this.jvxlData.thisSet >= 0) for (var i = 0; i < this.vc; i++) if (this.vertexSets[i] != this.jvxlData.thisSet) bs.clear (i);

return bs;
});
Clazz_defineMethod (c$, "updateCoordinates", 
function (m, bs) {
var doUpdate = (bs == null);
if (!doUpdate) for (var i = 0; i < this.connections.length; i++) if (this.connections[i] >= 0 && bs.get (this.connections[i])) {
doUpdate = true;
break;
}
if (!doUpdate) return;
if (this.mat4 == null) this.mat4 = JU.M4.newM4 (null);
this.mat4.mul2 (m, this.mat4);
this.recalcAltVertices = true;
}, "JU.M4,JU.BS");
});
Clazz_declarePackage ("J.jvxl.readers");
Clazz_load (["J.jvxl.readers.SurfaceReader"], "J.jvxl.readers.VolumeDataReader", ["JU.AU", "$.SB", "J.jvxl.data.JvxlCoder", "JU.Logger"], function () {
c$ = Clazz_decorateAsClass (function () {
this.dataType = 0;
this.precalculateVoxelData = false;
this.allowMapData = false;
this.point = null;
this.ptsPerAngstrom = 0;
this.maxGrid = 0;
this.atomDataServer = null;
this.useOriginStepsPoints = false;
Clazz_instantialize (this, arguments);
}, J.jvxl.readers, "VolumeDataReader", J.jvxl.readers.SurfaceReader);
Clazz_makeConstructor (c$, 
function () {
Clazz_superConstructor (this, J.jvxl.readers.VolumeDataReader, []);
});
Clazz_defineMethod (c$, "initVDR", 
function (sg) {
this.initSR (sg);
this.useOriginStepsPoints = (this.params.origin != null && this.params.points != null && this.params.steps != null);
this.dataType = this.params.dataType;
this.precalculateVoxelData = true;
this.allowMapData = true;
}, "J.jvxl.readers.SurfaceGenerator");
Clazz_defineMethod (c$, "setup", 
function (isMapData) {
this.jvxlFileHeaderBuffer =  new JU.SB ().append ("volume data read from file\n\n");
J.jvxl.data.JvxlCoder.jvxlCreateHeaderWithoutTitleOrAtoms (this.volumeData, this.jvxlFileHeaderBuffer);
}, "~B");
Clazz_overrideMethod (c$, "readVolumeParameters", 
function (isMapData) {
this.setup (isMapData);
this.initializeVolumetricData ();
return true;
}, "~B");
Clazz_overrideMethod (c$, "readVolumeData", 
function (isMapData) {
try {
this.readSurfaceData (isMapData);
} catch (e) {
if (Clazz_exceptionOf (e, Exception)) {
System.out.println (e.toString ());
return false;
} else {
throw e;
}
}
return true;
}, "~B");
Clazz_defineMethod (c$, "readVoxelDataIndividually", 
function (isMapData) {
if (isMapData && !this.allowMapData) return;
if (!isMapData || this.volumeData.sr != null) {
this.volumeData.setVoxelDataAsArray (this.voxelData = null);
return;
}this.newVoxelDataCube ();
for (var x = 0; x < this.nPointsX; ++x) {
var plane = JU.AU.newFloat2 (this.nPointsY);
this.voxelData[x] = plane;
var ptyz = 0;
for (var y = 0; y < this.nPointsY; ++y) {
var strip = plane[y] =  Clazz_newFloatArray (this.nPointsZ, 0);
for (var z = 0; z < this.nPointsZ; ++z, ++ptyz) {
strip[z] = this.getValue (x, y, z, ptyz);
}
}
}
}, "~B");
Clazz_defineMethod (c$, "setVolumeData", 
function () {
});
Clazz_defineMethod (c$, "setVolumeDataParams", 
function () {
if (this.params.volumeData != null) {
this.setVolumeDataV (this.params.volumeData);
return true;
}if (!this.useOriginStepsPoints) {
return false;
}this.volumetricOrigin.setT (this.params.origin);
this.volumetricVectors[0].set (this.params.steps.x, 0, 0);
this.volumetricVectors[1].set (0, this.params.steps.y, 0);
this.volumetricVectors[2].set (0, 0, this.params.steps.z);
this.voxelCounts[0] = Clazz_floatToInt (this.params.points.x);
this.voxelCounts[1] = Clazz_floatToInt (this.params.points.y);
this.voxelCounts[2] = Clazz_floatToInt (this.params.points.z);
if (this.voxelCounts[0] < 1 || this.voxelCounts[1] < 1 || this.voxelCounts[2] < 1) return false;
this.showGridInfo ();
return true;
});
Clazz_defineMethod (c$, "showGridInfo", 
function () {
JU.Logger.info ("grid origin  = " + this.params.origin);
JU.Logger.info ("grid steps   = " + this.params.steps);
JU.Logger.info ("grid points  = " + this.params.points);
this.ptTemp.x = this.params.steps.x * this.params.points.x;
this.ptTemp.y = this.params.steps.y * this.params.points.y;
this.ptTemp.z = this.params.steps.z * this.params.points.z;
JU.Logger.info ("grid lengths = " + this.ptTemp);
this.ptTemp.add (this.params.origin);
JU.Logger.info ("grid max xyz = " + this.ptTemp);
});
Clazz_defineMethod (c$, "setVoxelRange", 
function (index, min, max, ptsPerAngstrom, gridMax, minPointsPerAngstrom) {
var nGrid;
var d;
if (min >= max) {
min = -10;
max = 10;
}var range = max - min;
var resolution = this.params.resolution;
if (resolution != 3.4028235E38) ptsPerAngstrom = resolution;
nGrid = Clazz_doubleToInt (Math.floor (range * ptsPerAngstrom)) + 1;
if (nGrid > gridMax) {
if ((this.dataType & 256) > 0) {
if (resolution == 3.4028235E38) {
if (!this.isQuiet) JU.Logger.info ("Maximum number of voxels for index=" + index + " exceeded (" + nGrid + ") -- set to " + gridMax);
nGrid = gridMax;
} else {
if (!this.isQuiet) JU.Logger.info ("Warning -- high number of grid points: " + nGrid);
}} else if (resolution == 3.4028235E38) {
nGrid = gridMax;
}}ptsPerAngstrom = (nGrid - 1) / range;
if (ptsPerAngstrom < minPointsPerAngstrom) {
ptsPerAngstrom = minPointsPerAngstrom;
nGrid = Clazz_doubleToInt (Math.floor (ptsPerAngstrom * range + 1));
ptsPerAngstrom = (nGrid - 1) / range;
}d = this.volumeData.volumetricVectorLengths[index] = 1 / ptsPerAngstrom;
this.voxelCounts[index] = nGrid;
if (!this.isQuiet) JU.Logger.info ("isosurface resolution for axis " + (index + 1) + " set to " + ptsPerAngstrom + " points/Angstrom; " + this.voxelCounts[index] + " voxels");
switch (index) {
case 0:
this.volumetricVectors[0].set (d, 0, 0);
this.volumetricOrigin.x = min;
break;
case 1:
this.volumetricVectors[1].set (0, d, 0);
this.volumetricOrigin.y = min;
break;
case 2:
this.volumetricVectors[2].set (0, 0, d);
this.volumetricOrigin.z = min;
if (this.isEccentric) this.eccentricityMatrix.rotate (this.volumetricOrigin);
if (this.center != null && this.center.x != 3.4028235E38) this.volumetricOrigin.add (this.center);
}
if (this.isEccentric) this.eccentricityMatrix.rotate (this.volumetricVectors[index]);
return this.voxelCounts[index];
}, "~N,~N,~N,~N,~N,~N");
Clazz_overrideMethod (c$, "readSurfaceData", 
function (isMapData) {
this.readSurfaceDataVDR (isMapData);
}, "~B");
Clazz_defineMethod (c$, "readSurfaceDataVDR", 
function (isMapData) {
if (this.isProgressive && !isMapData) {
this.nDataPoints = this.volumeData.setVoxelCounts (this.nPointsX, this.nPointsY, this.nPointsZ);
this.voxelData = null;
return;
}if (this.precalculateVoxelData) this.generateCube ();
 else this.readVoxelDataIndividually (isMapData);
}, "~B");
Clazz_defineMethod (c$, "generateCube", 
function () {
JU.Logger.info ("data type: user volumeData");
JU.Logger.info ("voxel grid origin:" + this.volumetricOrigin);
for (var i = 0; i < 3; ++i) JU.Logger.info ("voxel grid vector:" + this.volumetricVectors[i]);

JU.Logger.info ("Read " + this.nPointsX + " x " + this.nPointsY + " x " + this.nPointsZ + " data points");
});
Clazz_overrideMethod (c$, "closeReader", 
function () {
});
});
Clazz_declarePackage ("J.jvxl.readers");
Clazz_load (["J.jvxl.readers.VolumeDataReader", "JU.BS", "$.P3", "$.P3i", "J.atomdata.AtomData"], "J.jvxl.readers.AtomDataReader", ["java.lang.Float", "java.util.Date", "JU.AU", "$.SB", "$.V3", "J.atomdata.RadiusData", "J.c.VDW", "J.jvxl.data.JvxlCoder", "JU.BSUtil", "$.Logger", "$.Txt"], function () {
c$ = Clazz_decorateAsClass (function () {
this.maxDistance = 0;
this.contactPair = null;
this.fileName = null;
this.fileDotModel = null;
this.modelIndex = 0;
this.atomData = null;
this.atomXyz = null;
this.atomRadius = null;
this.atomProp = null;
this.atomNo = null;
this.atomIndex = null;
this.myIndex = null;
this.ac = 0;
this.myAtomCount = 0;
this.nearbyAtomCount = 0;
this.firstNearbyAtom = 0;
this.bsMySelected = null;
this.bsMyIgnored = null;
this.bsNearby = null;
this.doAddHydrogens = false;
this.havePlane = false;
this.doUseIterator = false;
this.theProperty = 0;
this.haveOneProperty = false;
this.minPtsPerAng = 0;
this.sr = 0;
this.rs = null;
this.rs2 = null;
this.maxRS = 0;
this.thisPlane = null;
this.thisAtomSet = null;
this.thisX = 0;
this.margin = 0;
this.vl0 = 0;
this.vl1 = 0;
this.vl2 = 0;
this.bsSurfaceVoxels = null;
this.validSpheres = null;
this.noFaceSpheres = null;
this.voxelSource = null;
this.ptY0 = null;
this.ptZ0 = null;
this.pt0 = null;
this.pt1 = null;
this.ptV = null;
Clazz_instantialize (this, arguments);
}, J.jvxl.readers, "AtomDataReader", J.jvxl.readers.VolumeDataReader);
Clazz_prepareFields (c$, function () {
this.atomData =  new J.atomdata.AtomData ();
this.bsMySelected =  new JU.BS ();
this.bsMyIgnored =  new JU.BS ();
this.ptY0 =  new JU.P3 ();
this.ptZ0 =  new JU.P3 ();
this.pt0 =  new JU.P3i ();
this.pt1 =  new JU.P3i ();
this.ptV =  new JU.P3 ();
});
Clazz_makeConstructor (c$, 
function () {
Clazz_superConstructor (this, J.jvxl.readers.AtomDataReader, []);
});
Clazz_defineMethod (c$, "initADR", 
function (sg) {
this.initVDR (sg);
this.precalculateVoxelData = true;
this.atomDataServer = sg.getAtomDataServer ();
}, "J.jvxl.readers.SurfaceGenerator");
Clazz_overrideMethod (c$, "setup", 
function (isMapData) {
this.setup2 ();
}, "~B");
Clazz_defineMethod (c$, "setup2", 
function () {
this.contactPair = this.params.contactPair;
this.doAddHydrogens = (this.sg.getAtomDataServer () != null && this.params.addHydrogens);
this.modelIndex = this.params.modelIndex;
if (this.params.bsIgnore != null) this.bsMyIgnored = this.params.bsIgnore;
if (this.params.volumeData != null) {
this.setVolumeDataV (this.params.volumeData);
this.setBBox (this.volumeData.volumetricOrigin, 0);
this.ptV.setT (this.volumeData.volumetricOrigin);
for (var i = 0; i < 3; i++) this.ptV.scaleAdd2 (this.volumeData.voxelCounts[i] - 1, this.volumeData.volumetricVectors[i], this.ptV);

this.setBBox (this.ptV, 0);
}this.havePlane = (this.params.thePlane != null);
if (this.havePlane) this.volumeData.setPlaneParameters (this.params.thePlane);
});
Clazz_defineMethod (c$, "markPlaneVoxels", 
function (p, r) {
for (var i = 0, pt = this.thisX * this.yzCount, pt1 = pt + this.yzCount; pt < pt1; pt++, i++) {
this.volumeData.getPoint (pt, this.ptV);
this.thisPlane[i] = this.ptV.distance (p) - r;
}
}, "JU.P3,~N");
Clazz_defineMethod (c$, "setVolumeForPlane", 
function () {
if (this.useOriginStepsPoints) {
this.xyzMin = JU.P3.newP (this.params.origin);
this.xyzMax = JU.P3.newP (this.params.origin);
this.xyzMax.add3 ((this.params.points.x - 1) * this.params.steps.x, (this.params.points.y - 1) * this.params.steps.y, (this.params.points.z - 1) * this.params.steps.z);
} else if (this.params.boundingBox == null) {
this.getAtoms (this.params.bsSelected, false, true, false, false, false, false, this.params.mep_marginAngstroms);
if (this.xyzMin == null) {
this.xyzMin = JU.P3.new3 (-10, -10, -10);
this.xyzMax = JU.P3.new3 (10, 10, 10);
}} else {
this.xyzMin = JU.P3.newP (this.params.boundingBox[0]);
this.xyzMax = JU.P3.newP (this.params.boundingBox[1]);
}this.setRanges (this.params.plane_ptsPerAngstrom, this.params.plane_gridMax, 0);
});
Clazz_defineMethod (c$, "getAtoms", 
function (bsSelected, doAddHydrogens, getRadii, getMolecules, getAllModels, addNearbyAtoms, getAtomMinMax, marginAtoms) {
if (addNearbyAtoms) getRadii = true;
if (getRadii) {
if (this.params.atomRadiusData == null) this.params.atomRadiusData =  new J.atomdata.RadiusData (null, 1, J.atomdata.RadiusData.EnumType.FACTOR, J.c.VDW.AUTO);
this.atomData.radiusData = this.params.atomRadiusData;
this.atomData.radiusData.valueExtended = this.params.solventExtendedAtomRadius;
if (doAddHydrogens) this.atomData.radiusData.vdwType = J.c.VDW.NOJMOL;
}this.atomData.modelIndex = this.modelIndex;
this.atomData.bsSelected = bsSelected;
this.atomData.bsIgnored = this.bsMyIgnored;
this.sg.fillAtomData (this.atomData, 1 | (getAllModels ? 16 : 0) | (getMolecules ? 4 : 0) | (getRadii ? 2 : 0));
if (this.doUseIterator) this.atomData.bsSelected = null;
this.ac = this.atomData.ac;
this.modelIndex = this.atomData.firstModelIndex;
var needRadius = false;
for (var i = 0; i < this.ac; i++) {
if ((bsSelected == null || bsSelected.get (i)) && (!this.bsMyIgnored.get (i))) {
if (this.havePlane && Math.abs (this.volumeData.distancePointToPlane (this.atomData.atomXyz[i])) > 2 * (this.atomData.atomRadius[i] = this.getWorkingRadius (i, marginAtoms))) continue;
this.bsMySelected.set (i);
needRadius = !this.havePlane;
}if (getRadii && (addNearbyAtoms || needRadius)) this.atomData.atomRadius[i] = this.getWorkingRadius (i, marginAtoms);
}
var rH = (getRadii && doAddHydrogens ? this.getWorkingRadius (-1, marginAtoms) : 0);
this.myAtomCount = JU.BSUtil.cardinalityOf (this.bsMySelected);
var atomSet = JU.BSUtil.copy (this.bsMySelected);
var nH = 0;
this.atomProp = null;
this.theProperty = 3.4028235E38;
this.haveOneProperty = false;
var props = this.params.theProperty;
if (this.myAtomCount > 0) {
var hAtoms = null;
if (doAddHydrogens) {
this.atomData.bsSelected = atomSet;
this.atomDataServer.fillAtomData (this.atomData, 8);
hAtoms =  new Array (nH = this.atomData.hydrogenAtomCount);
for (var i = 0; i < this.atomData.hAtoms.length; i++) if (this.atomData.hAtoms[i] != null) for (var j = this.atomData.hAtoms[i].length; --j >= 0; ) hAtoms[--nH] = this.atomData.hAtoms[i][j];


nH = hAtoms.length;
JU.Logger.info (nH + " attached hydrogens added");
}var n = nH + this.myAtomCount;
if (getRadii) this.atomRadius =  Clazz_newFloatArray (n, 0);
this.atomXyz =  new Array (n);
if (this.params.theProperty != null) this.atomProp =  Clazz_newFloatArray (n, 0);
this.atomNo =  Clazz_newIntArray (n, 0);
this.atomIndex =  Clazz_newIntArray (n, 0);
this.myIndex =  Clazz_newIntArray (this.ac, 0);
for (var i = 0; i < nH; i++) {
if (getRadii) this.atomRadius[i] = rH;
this.atomXyz[i] = hAtoms[i];
this.atomNo[i] = -1;
if (this.atomProp != null) this.addAtomProp (i, NaN);
}
this.myAtomCount = nH;
for (var i = atomSet.nextSetBit (0); i >= 0; i = atomSet.nextSetBit (i + 1)) {
if (this.atomProp != null) this.addAtomProp (this.myAtomCount, (props != null && i < props.length ? props[i] : NaN));
this.atomXyz[this.myAtomCount] = this.atomData.atomXyz[i];
this.atomNo[this.myAtomCount] = this.atomData.atomicNumber[i];
this.atomIndex[this.myAtomCount] = i;
this.myIndex[i] = this.myAtomCount;
if (getRadii) this.atomRadius[this.myAtomCount] = this.atomData.atomRadius[i];
this.myAtomCount++;
}
}this.firstNearbyAtom = this.myAtomCount;
JU.Logger.info (this.myAtomCount + " atoms will be used in the surface calculation");
if (this.myAtomCount == 0) {
this.setBBox (JU.P3.new3 (10, 10, 10), 0);
this.setBBox (JU.P3.new3 (-10, -10, -10), 0);
}for (var i = 0; i < this.myAtomCount; i++) this.setBBox (this.atomXyz[i], getRadii ? this.atomRadius[i] + 0.5 : 0);

if (!Float.isNaN (this.params.scale)) {
var v = JU.V3.newVsub (this.xyzMax, this.xyzMin);
v.scale (0.5);
this.xyzMin.add (v);
v.scale (this.params.scale);
this.xyzMax.add2 (this.xyzMin, v);
this.xyzMin.sub (v);
}if (!addNearbyAtoms || this.myAtomCount == 0) return;
var pt =  new JU.P3 ();
this.bsNearby =  new JU.BS ();
for (var i = 0; i < this.ac; i++) {
if (atomSet.get (i) || this.bsMyIgnored.get (i)) continue;
var rA = this.atomData.atomRadius[i];
if (this.params.thePlane != null && Math.abs (this.volumeData.distancePointToPlane (this.atomData.atomXyz[i])) > 2 * rA) continue;
if (this.params.theProperty != null) rA += this.maxDistance;
pt = this.atomData.atomXyz[i];
if (pt.x + rA > this.xyzMin.x && pt.x - rA < this.xyzMax.x && pt.y + rA > this.xyzMin.y && pt.y - rA < this.xyzMax.y && pt.z + rA > this.xyzMin.z && pt.z - rA < this.xyzMax.z) {
this.bsNearby.set (i);
this.nearbyAtomCount++;
}}
var nAtoms = this.myAtomCount;
if (this.nearbyAtomCount != 0) {
nAtoms += this.nearbyAtomCount;
this.atomRadius = JU.AU.arrayCopyF (this.atomRadius, nAtoms);
this.atomXyz = JU.AU.arrayCopyObject (this.atomXyz, nAtoms);
if (this.atomIndex != null) this.atomIndex = JU.AU.arrayCopyI (this.atomIndex, nAtoms);
if (props != null) this.atomProp = JU.AU.arrayCopyF (this.atomProp, nAtoms);
for (var i = this.bsNearby.nextSetBit (0); i >= 0; i = this.bsNearby.nextSetBit (i + 1)) {
if (props != null) this.addAtomProp (this.myAtomCount, props[i]);
this.myIndex[i] = this.myAtomCount;
this.atomIndex[this.myAtomCount] = i;
this.atomXyz[this.myAtomCount] = this.atomData.atomXyz[i];
this.atomRadius[this.myAtomCount++] = this.atomData.atomRadius[i];
}
}if (getRadii) this.setRadii ();
this.haveOneProperty = (!Float.isNaN (this.theProperty));
}, "JU.BS,~B,~B,~B,~B,~B,~B,~N");
Clazz_defineMethod (c$, "setRadii", 
function () {
if (this.rs != null) return;
this.maxRS = 0;
this.rs =  Clazz_newFloatArray (this.myAtomCount, 0);
this.rs2 =  Clazz_newFloatArray (this.myAtomCount, 0);
for (var i = 0; i < this.myAtomCount; i++) {
var r = this.rs[i] = this.atomRadius[i] + this.sr;
if (r > this.maxRS) this.maxRS = r;
this.rs2[i] = this.rs[i] * this.rs[i];
}
});
Clazz_defineMethod (c$, "addAtomProp", 
 function (i, f) {
this.atomProp[i] = f;
if (!Float.isNaN (this.theProperty)) if (f != this.theProperty) this.theProperty = (this.theProperty == 3.4028235E38 ? f : NaN);
}, "~N,~N");
Clazz_defineMethod (c$, "getWorkingRadius", 
 function (i, marginAtoms) {
var r = (i < 0 ? this.atomData.hAtomRadius : this.atomData.atomRadius[i]);
return (Float.isNaN (marginAtoms) ? Math.max (r, 0.1) : r + marginAtoms);
}, "~N,~N");
Clazz_defineMethod (c$, "setHeader", 
function (calcType, line2) {
this.jvxlFileHeaderBuffer =  new JU.SB ();
if (this.atomData.programInfo != null) this.jvxlFileHeaderBuffer.append ("#created by ").append (this.atomData.programInfo).append (" on ").append ("" +  new java.util.Date ()).append ("\n");
this.jvxlFileHeaderBuffer.append (calcType).append ("\n").append (line2).append ("\n");
}, "~S,~S");
Clazz_defineMethod (c$, "setRanges", 
function (ptsPerAngstrom, maxGrid, minPtsPerAng) {
if (this.xyzMin == null) return;
this.ptsPerAngstrom = ptsPerAngstrom;
this.maxGrid = maxGrid;
this.minPtsPerAng = minPtsPerAng;
this.setVolumeData ();
J.jvxl.data.JvxlCoder.jvxlCreateHeader (this.volumeData, this.jvxlFileHeaderBuffer);
}, "~N,~N,~N");
Clazz_overrideMethod (c$, "setVolumeData", 
function () {
this.setVolumeDataADR ();
});
Clazz_defineMethod (c$, "setVolumeDataADR", 
function () {
if (!this.setVolumeDataParams ()) {
this.setVoxelRange (0, this.xyzMin.x, this.xyzMax.x, this.ptsPerAngstrom, this.maxGrid, this.minPtsPerAng);
this.setVoxelRange (1, this.xyzMin.y, this.xyzMax.y, this.ptsPerAngstrom, this.maxGrid, this.minPtsPerAng);
this.setVoxelRange (2, this.xyzMin.z, this.xyzMax.z, this.ptsPerAngstrom, this.maxGrid, this.minPtsPerAng);
}});
Clazz_defineMethod (c$, "fixTitleLine", 
function (iLine) {
if (this.params.title == null) return false;
var line = this.params.title[iLine];
if (line.indexOf ("%F") > 0) line = this.params.title[iLine] = JU.Txt.formatStringS (line, "F", this.atomData.fileName);
if (line.indexOf ("%M") > 0) this.params.title[iLine] = JU.Txt.formatStringS (line, "M", this.atomData.modelName);
return true;
}, "~N");
Clazz_defineMethod (c$, "setVertexSource", 
function () {
if (this.meshDataServer != null) this.meshDataServer.fillMeshData (this.meshData, 1, null);
if (this.params.vertexSource != null) {
this.params.vertexSource = JU.AU.arrayCopyI (this.params.vertexSource, this.meshData.vc);
for (var i = 0; i < this.meshData.vc; i++) this.params.vertexSource[i] = Math.abs (this.params.vertexSource[i]) - 1;

}});
Clazz_defineMethod (c$, "resetPlane", 
function (value) {
for (var i = 0; i < this.yzCount; i++) this.thisPlane[i] = value;

}, "~N");
Clazz_defineMethod (c$, "resetVoxelData", 
function (value) {
for (var x = 0; x < this.nPointsX; ++x) for (var y = 0; y < this.nPointsY; ++y) for (var z = 0; z < this.nPointsZ; ++z) this.voxelData[x][y][z] = value;



}, "~N");
Clazz_defineMethod (c$, "getVoxel", 
 function (i, j, k, ipt) {
return (this.isProgressive ? this.thisPlane[ipt % this.yzCount] : this.voxelData[i][j][k]);
}, "~N,~N,~N,~N");
Clazz_defineMethod (c$, "unsetVoxelData", 
function () {
this.unsetVoxelData2 ();
});
Clazz_defineMethod (c$, "unsetVoxelData2", 
function () {
if (this.isProgressive) for (var i = 0; i < this.yzCount; i++) {
if (this.thisPlane[i] == 3.4028235E38) this.thisPlane[i] = NaN;
}
 else for (var x = 0; x < this.nPointsX; ++x) for (var y = 0; y < this.nPointsY; ++y) for (var z = 0; z < this.nPointsZ; ++z) if (this.voxelData[x][y][z] == 3.4028235E38) this.voxelData[x][y][z] = NaN;



});
Clazz_defineMethod (c$, "setGridLimitsForAtom", 
function (ptA, rA, pt0, pt1) {
rA += this.margin;
this.volumeData.xyzToVoxelPt (ptA.x, ptA.y, ptA.z, pt0);
var x = Clazz_doubleToInt (Math.floor (rA / this.volumeData.volumetricVectorLengths[0]));
var y = Clazz_doubleToInt (Math.floor (rA / this.volumeData.volumetricVectorLengths[1]));
var z = Clazz_doubleToInt (Math.floor (rA / this.volumeData.volumetricVectorLengths[2]));
pt1.set (pt0.x + x, pt0.y + y, pt0.z + z);
pt0.set (pt0.x - x, pt0.y - y, pt0.z - z);
pt0.x = Math.max (pt0.x - 1, 0);
pt0.y = Math.max (pt0.y - 1, 0);
pt0.z = Math.max (pt0.z - 1, 0);
pt1.x = Math.min (pt1.x + 1, this.nPointsX);
pt1.y = Math.min (pt1.y + 1, this.nPointsY);
pt1.z = Math.min (pt1.z + 1, this.nPointsZ);
}, "JU.P3,~N,JU.P3i,JU.P3i");
Clazz_defineMethod (c$, "getAtomMinMax", 
function (bs, bsAtomMinMax) {
for (var i = 0; i < this.nPointsX; i++) bsAtomMinMax[i] =  new JU.BS ();

for (var iAtom = this.myAtomCount; --iAtom >= 0; ) {
if (bs != null && !bs.get (iAtom)) continue;
this.setGridLimitsForAtom (this.atomXyz[iAtom], this.atomRadius[iAtom], this.pt0, this.pt1);
for (var i = this.pt0.x; i < this.pt1.x; i++) bsAtomMinMax[i].set (iAtom);

}
}, "JU.BS,~A");
Clazz_defineMethod (c$, "markSphereVoxels", 
function (r0, distance) {
var isWithin = (distance != 3.4028235E38 && this.point != null);
var v0 = this.volumetricVectors[0];
var v1 = this.volumetricVectors[1];
var v2 = this.volumetricVectors[2];
for (var iAtom = this.thisAtomSet.nextSetBit (0); iAtom >= 0; iAtom = this.thisAtomSet.nextSetBit (iAtom + 1)) {
if (!this.havePlane && this.validSpheres != null && !this.validSpheres.get (iAtom)) continue;
var isSurface = (this.noFaceSpheres != null && this.noFaceSpheres.get (iAtom));
var isNearby = (iAtom >= this.firstNearbyAtom);
var ptA = this.atomXyz[iAtom];
var rA = this.atomRadius[iAtom];
if (isWithin && ptA.distance (this.point) > distance + rA + 0.5) continue;
var rA0 = rA + r0;
this.setGridLimitsForAtom (ptA, rA0, this.pt0, this.pt1);
if (this.isProgressive) {
this.pt0.x = this.thisX;
this.pt1.x = this.thisX + 1;
}this.volumeData.voxelPtToXYZ (this.pt0.x, this.pt0.y, this.pt0.z, this.ptV);
for (var i = this.pt0.x; i < this.pt1.x; i++, this.ptV.add2 (v0, this.ptY0)) {
this.ptY0.setT (this.ptV);
for (var j = this.pt0.y; j < this.pt1.y; j++, this.ptV.add2 (v1, this.ptZ0)) {
this.ptZ0.setT (this.ptV);
for (var k = this.pt0.z; k < this.pt1.z; k++, this.ptV.add (v2)) {
var value = this.ptV.distance (ptA) - rA;
var ipt = this.volumeData.getPointIndex (i, j, k);
if ((r0 == 0 || value <= rA0) && value < this.getVoxel (i, j, k, ipt)) {
if (isNearby || isWithin && this.ptV.distance (this.point) > distance) value = NaN;
this.setVoxel (i, j, k, ipt, value);
if (!Float.isNaN (value)) {
if (this.voxelSource != null) this.voxelSource[ipt] = iAtom + 1;
if (value < 0 && isSurface) this.bsSurfaceVoxels.set (ipt);
}}}
}
}
}
}, "~N,~N");
Clazz_defineMethod (c$, "setVoxel", 
function (i, j, k, ipt, value) {
if (this.isProgressive) this.thisPlane[ipt % this.yzCount] = value;
 else this.voxelData[i][j][k] = value;
}, "~N,~N,~N,~N,~N");
});
Clazz_declarePackage ("J.jvxl.readers");
Clazz_load (["JU.P3", "J.jvxl.readers.AtomDataReader", "JU.P4", "$.V3"], "J.jvxl.readers.IsoSolventReader", ["java.lang.Float", "java.util.Hashtable", "JU.BS", "$.Lst", "J.jvxl.data.MeshData", "JU.BSUtil", "$.Logger", "$.Measure", "$.MeshSurface", "$.TempArray"], function () {
c$ = Clazz_decorateAsClass (function () {
this.cavityRadius = 0;
this.envelopeRadius = 0;
this.dots = null;
this.doCalculateTroughs = false;
this.isCavity = false;
this.isPocket = false;
this.iter = null;
this.bsSurfacePoints = null;
this.bsSurfaceDone = null;
this.bsLocale = null;
this.htEdges = null;
this.vEdges = null;
this.vFaces = null;
this.ptS1 = null;
this.ptS2 = null;
this.vTemp3 = null;
this.vTemp = null;
this.plane = null;
this.ptTemp2 = null;
this.vTemp2 = null;
this.p = null;
this.bsAtomMinMax = null;
this.isSurfacePoint = false;
this.iAtomSurface = 0;
this.nTest = 0;
if (!Clazz_isClassDefined ("J.jvxl.readers.IsoSolventReader.Edge")) {
J.jvxl.readers.IsoSolventReader.$IsoSolventReader$Edge$ ();
}
if (!Clazz_isClassDefined ("J.jvxl.readers.IsoSolventReader.Face")) {
J.jvxl.readers.IsoSolventReader.$IsoSolventReader$Face$ ();
}
this.rAS = 0;
this.rBS = 0;
this.rAS2 = 0;
this.rBS2 = 0;
this.dAB = 0;
this.dAB2 = 0;
this.ecosASB2 = 0;
Clazz_instantialize (this, arguments);
}, J.jvxl.readers, "IsoSolventReader", J.jvxl.readers.AtomDataReader);
Clazz_prepareFields (c$, function () {
this.ptS1 =  new JU.P3 ();
this.ptS2 =  new JU.P3 ();
this.vTemp3 =  new JU.V3 ();
this.vTemp =  new JU.V3 ();
this.plane =  new JU.P4 ();
this.ptTemp2 =  new JU.P3 ();
this.vTemp2 =  new JU.V3 ();
this.p =  new JU.P3 ();
});
Clazz_makeConstructor (c$, 
function () {
Clazz_superConstructor (this, J.jvxl.readers.IsoSolventReader, []);
});
Clazz_overrideMethod (c$, "init", 
function (sg) {
this.initADR (sg);
}, "J.jvxl.readers.SurfaceGenerator");
Clazz_overrideMethod (c$, "readVolumeParameters", 
function (isMapData) {
this.setup (isMapData);
this.initializeVolumetricData ();
this.volumeData.setUnitVectors ();
this.vl0 = this.volumeData.volumetricVectorLengths[0];
this.vl1 = this.volumeData.volumetricVectorLengths[1];
this.vl2 = this.volumeData.volumetricVectorLengths[2];
if (this.isProgressive) {
this.volumeData.getYzCount ();
this.bsAtomMinMax =  new Array (this.nPointsX);
this.getAtomMinMax (null, this.bsAtomMinMax);
this.voxelSource =  Clazz_newIntArray (this.volumeData.nPoints, 0);
}return true;
}, "~B");
Clazz_overrideMethod (c$, "setup", 
function (isMapData) {
this.setup2 ();
if (this.contactPair == null) {
this.cavityRadius = this.params.cavityRadius;
this.envelopeRadius = this.params.envelopeRadius;
this.sr = this.params.solventRadius;
this.point = this.params.point;
this.isCavity = (this.params.isCavity && this.meshDataServer != null);
this.isPocket = (this.params.pocket != null && this.meshDataServer != null);
this.doCalculateTroughs = (!isMapData && this.atomDataServer != null && !this.isCavity && this.sr > 0 && (this.dataType == 1195 || this.dataType == 1203));
this.doUseIterator = this.doCalculateTroughs;
this.getAtoms (this.params.bsSelected, this.doAddHydrogens, true, false, false, true, false, NaN);
if (this.isCavity || this.isPocket) this.dots = this.meshDataServer.calculateGeodesicSurface (this.bsMySelected, this.envelopeRadius);
this.setHeader ("solvent/molecular surface", this.params.calculationType);
if (this.havePlane || !isMapData) {
var minPtsPerAng = 0;
this.setRanges (this.params.solvent_ptsPerAngstrom, this.params.solvent_gridMax, minPtsPerAng);
this.volumeData.getYzCount ();
this.margin = this.volumeData.maxGrid * 2.0;
}if (this.bsNearby != null) this.bsMySelected.or (this.bsNearby);
} else if (!isMapData) {
this.setVolumeData ();
}if (!this.doCalculateTroughs) {
if (isMapData) {
this.precalculateVoxelData = false;
this.volumeData.sr = this;
} else if (!this.isCavity) {
this.isProgressive = this.isXLowToHigh = true;
}}if (this.thisAtomSet == null) this.thisAtomSet = JU.BSUtil.setAll (this.myAtomCount);
}, "~B");
Clazz_overrideMethod (c$, "generateCube", 
function () {
if (this.isCavity && this.params.theProperty != null) return;
if (this.isCavity && this.dataType != 1205 && this.dataType != 1206) {
this.params.vertexSource = null;
this.newVoxelDataCube ();
this.resetVoxelData (3.4028235E38);
this.markSphereVoxels (this.cavityRadius, this.params.distance);
this.generateSolventCavity ();
this.resetVoxelData (3.4028235E38);
this.markSphereVoxels (0, NaN);
} else {
this.voxelSource =  Clazz_newIntArray (this.volumeData.nPoints, 0);
this.generateSolventCube ();
}this.unsetVoxelData ();
var info = this.params.slabInfo;
if (info != null) for (var i = 0; i < info.size (); i++) if ((info.get (i)[2]).booleanValue () && Clazz_instanceOf (info.get (i)[0], JU.P4)) {
this.volumeData.capData (info.get (i)[0], this.params.cutoff);
info.remove (i--);
}
});
Clazz_overrideMethod (c$, "getSurfacePointAndFraction", 
function (cutoff, isCutoffAbsolute, valueA, valueB, pointA, edgeVector, x, y, z, vA0, vB0, fReturn, ptReturn) {
var vA = this.marchingCubes.getLinearOffset (x, y, z, vA0);
var vB = this.marchingCubes.getLinearOffset (x, y, z, vB0);
this.isSurfacePoint = (this.bsSurfaceVoxels != null && (this.bsSurfaceVoxels.get (vA) || this.bsSurfaceVoxels.get (vB)));
if (J.jvxl.readers.IsoSolventReader.testLinear || this.voxelSource == null || this.voxelSource[vA] == 0 || this.voxelSource[vA] != this.voxelSource[vB]) return this.getSPF (cutoff, isCutoffAbsolute, valueA, valueB, pointA, edgeVector, x, y, z, vA, vB, fReturn, ptReturn);
var iAtom = Math.abs (valueA < valueB ? this.voxelSource[vA] : this.voxelSource[vB]);
this.iAtomSurface = this.atomIndex[iAtom - 1];
var fraction = fReturn[0] = JU.MeshSurface.getSphericalInterpolationFraction ((this.voxelSource[vA] < 0 ? this.sr : this.atomRadius[this.voxelSource[vA] - 1]), valueA, valueB, edgeVector.length ());
ptReturn.scaleAdd2 (fraction, edgeVector, pointA);
var diff = valueB - valueA;
return valueA + fraction * diff;
}, "~N,~B,~N,~N,JU.T3,JU.V3,~N,~N,~N,~N,~N,~A,JU.T3");
Clazz_overrideMethod (c$, "addVertexCopy", 
function (vertexXYZ, value, assocVertex, asCopy) {
var i = this.addVC (vertexXYZ, value, assocVertex, asCopy);
if (i < 0) return i;
if (this.isSurfacePoint) this.bsSurfacePoints.set (i);
if (this.params.vertexSource != null) this.params.vertexSource[i] = this.iAtomSurface;
return i;
}, "JU.T3,~N,~N,~B");
Clazz_overrideMethod (c$, "selectPocket", 
function (doExclude) {
if (this.meshDataServer != null) this.meshDataServer.fillMeshData (this.meshData, 1, null);
var v = this.meshData.vs;
var nVertices = this.meshData.vc;
var vv = this.meshData.vvs;
var nDots = this.dots.length;
for (var i = 0; i < nVertices; i++) {
for (var j = 0; j < nDots; j++) {
if (this.dots[j].distance (v[i]) < this.envelopeRadius) {
vv[i] = NaN;
continue;
}}
}
this.meshData.getSurfaceSet ();
var nSets = this.meshData.nSets;
var pocketSet = JU.BS.newN (nSets);
var ss;
for (var i = 0; i < nSets; i++) if ((ss = this.meshData.surfaceSet[i]) != null) for (var j = ss.nextSetBit (0); j >= 0; j = ss.nextSetBit (j + 1)) if (Float.isNaN (this.meshData.vvs[j])) {
pocketSet.set (i);
break;
}

for (var i = 0; i < nSets; i++) if (this.meshData.surfaceSet[i] != null && pocketSet.get (i) == doExclude) this.meshData.invalidateSurfaceSet (i);

this.updateSurfaceData ();
if (!doExclude) this.meshData.surfaceSet = null;
if (this.meshDataServer != null) {
this.meshDataServer.fillMeshData (this.meshData, 3, null);
this.meshData =  new J.jvxl.data.MeshData ();
}}, "~B");
Clazz_overrideMethod (c$, "postProcessVertices", 
function () {
this.setVertexSource ();
if (this.doCalculateTroughs && this.bsSurfacePoints != null) {
var bsAll =  new JU.BS ();
var bsSurfaces = this.meshData.getSurfaceSet ();
var bsSources = null;
var volumes = (this.isPocket ? null : this.meshData.calculateVolumeOrArea (-2147483648, false, false));
var minVolume = (1.5 * 3.141592653589793 * Math.pow (this.sr, 3));
var maxVolume = 0;
var maxIsNegative = false;
if (volumes != null && !this.isCavity) for (var i = 0; i < this.meshData.nSets; i++) {
var v = volumes[i];
if (Math.abs (v) > maxVolume) {
maxVolume = Math.abs (v);
maxIsNegative = (v < 0);
}}
var factor = (maxIsNegative ? -1 : 1);
for (var i = 0; i < this.meshData.nSets; i++) {
var bss = bsSurfaces[i];
if (bss.intersects (this.bsSurfacePoints)) {
if (volumes == null || volumes[i] * factor > minVolume) if (this.params.vertexSource != null) {
var bs =  new JU.BS ();
if (bsSources == null) bsSources =  new Array (bsSurfaces.length);
for (var j = bss.nextSetBit (0); j >= 0; j = bss.nextSetBit (j + 1)) {
var iatom = this.params.vertexSource[j];
if (iatom < 0) continue;
if (bsAll.get (iatom)) {
this.meshData.invalidateSurfaceSet (i);
break;
}bs.set (iatom);
}
bsAll.or (bs);
continue;
}}this.meshData.invalidateSurfaceSet (i);
}
this.updateSurfaceData ();
if (this.meshDataServer != null) {
this.meshDataServer.fillMeshData (this.meshData, 3, null);
this.meshData =  new J.jvxl.data.MeshData ();
}}if (this.params.thePlane != null && this.params.slabInfo == null) this.params.addSlabInfo (JU.TempArray.getSlabWithinRange (-100, 0));
});
Clazz_defineMethod (c$, "generateSolventCavity", 
 function () {
var bs = JU.BS.newN (this.nPointsX * this.nPointsY * this.nPointsZ);
var i = 0;
var nDots = this.dots.length;
var n = 0;
var d;
var r2 = this.envelopeRadius;
for (var x = 0; x < this.nPointsX; ++x) for (var y = 0; y < this.nPointsY; ++y) {
out : for (var z = 0; z < this.nPointsZ; ++z, ++i) if ((d = this.voxelData[x][y][z]) < 3.4028235E38 && d >= this.cavityRadius) {
this.volumeData.voxelPtToXYZ (x, y, z, this.ptV);
for (var j = 0; j < nDots; j++) {
if (this.dots[j].distance (this.ptV) < r2) continue out;
}
bs.set (i);
n++;
}
}

JU.Logger.info ("cavities include " + n + " voxel points");
this.atomRadius =  Clazz_newFloatArray (n, 0);
this.atomXyz =  new Array (n);
for (var x = 0, ipt = 0, apt = 0; x < this.nPointsX; ++x) for (var y = 0; y < this.nPointsY; ++y) for (var z = 0; z < this.nPointsZ; ++z) if (bs.get (ipt++)) {
this.volumeData.voxelPtToXYZ (x, y, z, (this.atomXyz[apt] =  new JU.P3 ()));
this.atomRadius[apt++] = this.voxelData[x][y][z];
}


this.myAtomCount = this.firstNearbyAtom = n;
this.thisAtomSet = JU.BSUtil.setAll (this.myAtomCount);
this.rs = null;
this.setRadii ();
});
Clazz_defineMethod (c$, "generateSolventCube", 
 function () {
if (this.dataType == 1205) return;
this.params.vertexSource =  Clazz_newIntArray (this.volumeData.nPoints, 0);
this.bsSurfaceDone =  new JU.BS ();
this.bsSurfaceVoxels =  new JU.BS ();
this.bsSurfacePoints =  new JU.BS ();
if (this.doCalculateTroughs) {
this.iter = this.atomDataServer.getSelectedAtomIterator (this.bsMySelected, true, false, false);
this.vEdges =  new JU.Lst ();
this.bsLocale =  new Array (this.myAtomCount);
this.htEdges =  new java.util.Hashtable ();
this.getEdges ();
JU.Logger.info (this.vEdges.size () + " edges");
this.vFaces =  new JU.Lst ();
this.getFaces ();
JU.Logger.info (this.vFaces.size () + " faces");
this.bsLocale = null;
this.htEdges = null;
this.iter.release ();
this.iter = null;
this.newVoxelDataCube ();
this.resetVoxelData (3.4028235E38);
this.markFaceVoxels (true);
this.markToroidVoxels ();
this.validSpheres.or (this.noFaceSpheres);
this.vEdges = null;
this.markFaceVoxels (false);
this.vFaces = null;
} else {
this.newVoxelDataCube ();
this.resetVoxelData (3.4028235E38);
}this.markSphereVoxels (0, this.doCalculateTroughs ? 3.4028235E38 : this.params.distance);
this.noFaceSpheres = null;
this.validSpheres = null;
});
Clazz_defineMethod (c$, "getEdges", 
 function () {
for (var iatomA = 0; iatomA < this.myAtomCount; iatomA++) this.bsLocale[iatomA] =  new JU.BS ();

for (var iatomA = 0; iatomA < this.myAtomCount; iatomA++) {
var ptA = this.atomXyz[iatomA];
var rA = this.rs[iatomA];
this.atomDataServer.setIteratorForAtom (this.iter, this.atomIndex[iatomA], rA + this.maxRS);
while (this.iter.hasNext ()) {
var iB = this.iter.next ();
var iatomB = this.myIndex[iB];
if (iatomA >= this.firstNearbyAtom && iatomB >= this.firstNearbyAtom) continue;
var ptB = this.atomXyz[iatomB];
var rB = this.rs[iatomB];
var dAB = ptA.distance (ptB);
if (dAB >= rA + rB) continue;
var edge = Clazz_innerTypeInstance (J.jvxl.readers.IsoSolventReader.Edge, this, null, this, iatomA, iatomB, dAB);
this.vEdges.addLast (edge);
this.bsLocale[iatomA].set (iatomB);
this.bsLocale[iatomB].set (iatomA);
this.htEdges.put (edge.toString (), edge);
}
}
});
Clazz_defineMethod (c$, "findEdge", 
function (i, j) {
return this.htEdges.get (i < j ? i + "_" + j : j + "_" + i);
}, "~N,~N");
Clazz_defineMethod (c$, "getFaces", 
 function () {
var bs =  new JU.BS ();
this.validSpheres =  new JU.BS ();
this.noFaceSpheres = JU.BSUtil.setAll (this.myAtomCount);
for (var i = this.vEdges.size (); --i >= 0; ) {
var edge = this.vEdges.get (i);
var ia = edge.ia;
var ib = edge.ib;
bs.clearAll ();
bs.or (this.bsLocale[ia]);
bs.and (this.bsLocale[ib]);
for (var ic = bs.nextSetBit (ib + 1); ic >= 0; ic = bs.nextSetBit (ic + 1)) {
if (this.getSolventPoints (edge, ia, ib, ic)) {
var f;
var isOK = false;
if ((f = this.validateFace (ia, ib, ic, edge, this.ptS1)) != null) {
this.vFaces.addLast (f);
isOK = true;
}if ((f = this.validateFace (ia, ib, ic, edge, this.ptS2)) != null) {
this.vFaces.addLast (f);
isOK = true;
}if (isOK) {
this.noFaceSpheres.clear (ia);
this.noFaceSpheres.clear (ib);
this.noFaceSpheres.clear (ic);
}}}
}
});
Clazz_defineMethod (c$, "validateFace", 
 function (ia, ib, ic, edge, ptS) {
this.atomDataServer.setIteratorForPoint (this.iter, this.modelIndex, ptS, this.maxRS);
var isValid = true;
while (this.iter.hasNext ()) {
var iia = this.iter.next ();
var iatom = this.myIndex[iia];
if (iatom == ia || iatom == ib || iatom == ic) continue;
var d = this.atomData.atomXyz[iia].distance (ptS);
if (d < this.atomData.atomRadius[iia] + this.sr) {
isValid = false;
break;
}}
var bc = this.findEdge (ib, ic);
var ca = this.findEdge (ia, ic);
var f = (isValid ? Clazz_innerTypeInstance (J.jvxl.readers.IsoSolventReader.Face, this, null, ia, ib, ic, ptS) : null);
edge.addFace (f);
bc.addFace (f);
ca.addFace (f);
if (!isValid) return null;
this.validSpheres.set (ia);
this.validSpheres.set (ib);
this.validSpheres.set (ic);
return f;
}, "~N,~N,~N,J.jvxl.readers.IsoSolventReader.Edge,JU.P3");
Clazz_defineMethod (c$, "markFaceVoxels", 
 function (firstPass) {
var bsThisPass =  new JU.BS ();
var v0 = this.volumetricVectors[0];
var v1 = this.volumetricVectors[1];
var v2 = this.volumetricVectors[2];
for (var fi = this.vFaces.size (); --fi >= 0; ) {
var f = this.vFaces.get (fi);
var ptA = this.atomXyz[f.ia];
var ptB = this.atomXyz[f.ib];
var ptC = this.atomXyz[f.ic];
var ptS = f.pS;
this.setGridLimitsForAtom (ptS, this.sr, this.pt0, this.pt1);
this.volumeData.voxelPtToXYZ (this.pt0.x, this.pt0.y, this.pt0.z, this.ptV);
for (var i = this.pt0.x; i < this.pt1.x; i++, this.ptV.add2 (v0, this.ptY0)) {
this.ptY0.setT (this.ptV);
for (var j = this.pt0.y; j < this.pt1.y; j++, this.ptV.add2 (v1, this.ptZ0)) {
this.ptZ0.setT (this.ptV);
for (var k = this.pt0.z; k < this.pt1.z; k++, this.ptV.add (v2)) {
var value = this.sr - this.ptV.distance (ptS);
var v = this.voxelData[i][j][k];
var ipt = this.volumeData.getPointIndex (i, j, k);
if (firstPass && value > 0) this.bsSurfaceDone.set (ipt);
if (JU.Measure.isInTetrahedron (this.ptV, ptA, ptB, ptC, ptS, this.plane, this.vTemp, this.vTemp2, this.vTemp3, false)) {
if (!firstPass ? !this.bsSurfaceDone.get (ipt) && value < 0 && value > -this.volumeData.maxGrid * 1.8 && (value > v) == bsThisPass.get (ipt) : (value > 0 && (v < 0 || v == 3.4028235E38 || (value > v) == bsThisPass.get (ipt)))) {
bsThisPass.set (ipt);
this.setVoxel (i, j, k, ipt, value);
if (this.voxelSource != null) this.voxelSource[ipt] = -1 - f.ia;
if (value > 0) {
this.bsSurfaceVoxels.set (ipt);
}}}}
}
}
}
}, "~B");
Clazz_defineMethod (c$, "markToroidVoxels", 
 function () {
var v0 = this.volumetricVectors[0];
var v1 = this.volumetricVectors[1];
var v2 = this.volumetricVectors[2];
for (var ei = this.vEdges.size (); --ei >= 0; ) {
var edge = this.vEdges.get (ei);
if (!edge.isValid ()) continue;
var ia = edge.ia;
var ib = edge.ib;
var ptA = this.atomXyz[ia];
var ptB = this.atomXyz[ib];
this.rAS = this.rs[ia];
this.rBS = this.rs[ib];
this.rAS2 = this.rs2[ia];
this.rBS2 = this.rs2[ib];
this.dAB = edge.d;
this.dAB2 = edge.d2;
this.ecosASB2 = edge.cosASB2;
this.setGridLimitsForAtom (edge, edge.maxr, this.pt0, this.pt1);
this.volumeData.voxelPtToXYZ (this.pt0.x, this.pt0.y, this.pt0.z, this.ptV);
for (var i = this.pt0.x; i < this.pt1.x; i++, this.ptV.add2 (v0, this.ptY0)) {
this.ptY0.setT (this.ptV);
for (var j = this.pt0.y; j < this.pt1.y; j++, this.ptV.add2 (v1, this.ptZ0)) {
this.ptZ0.setT (this.ptV);
for (var k = this.pt0.z; k < this.pt1.z; k++, this.ptV.add (v2)) {
var dVS = this.checkSpecialVoxel (ptA, ptB, this.ptV);
if (Float.isNaN (dVS)) continue;
var value = this.sr - dVS;
if (value < this.voxelData[i][j][k]) {
var ipt = this.volumeData.getPointIndex (i, j, k);
this.setVoxel (i, j, k, ipt, value);
if (this.voxelSource != null) this.voxelSource[ipt] = -1 - ia;
}}
}
}
}
});
Clazz_overrideMethod (c$, "unsetVoxelData", 
function () {
if (!this.havePlane) {
this.unsetVoxelData2 ();
return;
}if (this.isProgressive) for (var i = 0; i < this.yzCount; i++) {
if (this.thisPlane[i] < 0.001) {
} else {
this.thisPlane[i] = 0.001;
}}
 else for (var x = 0; x < this.nPointsX; ++x) for (var y = 0; y < this.nPointsY; ++y) for (var z = 0; z < this.nPointsZ; ++z) if (this.voxelData[x][y][z] < 0.001) {
} else {
this.voxelData[x][y][z] = 0.001;
}


});
Clazz_defineMethod (c$, "getSolventPoints", 
 function (edge, ia, ib, ic) {
var rAS = this.rs[ia];
var v = edge.v;
var cosAngleBAS = (edge.d2 + this.rs2[ia] - this.rs2[ib]) / (2 * edge.d * rAS);
var angleBAS = Math.acos (cosAngleBAS);
this.p.scaleAdd2 (cosAngleBAS * rAS, v, this.atomXyz[ia]);
JU.Measure.getPlaneThroughPoint (this.p, v, this.plane);
var dPS = (Math.sin (angleBAS) * rAS);
var ptC = this.atomXyz[ic];
var rCS = this.rs[ic];
var dCT = JU.Measure.distanceToPlane (this.plane, ptC);
if (Math.abs (dCT) >= rCS * 0.9) return false;
this.ptTemp.scaleAdd2 (-dCT, v, ptC);
var dpT = this.p.distance (this.ptTemp);
var dsp2 = dPS * dPS;
var dST2 = this.rs2[ic] - dCT * dCT;
var cosTheta = (dsp2 + dpT * dpT - dST2) / (2 * dPS * dpT);
if (Math.abs (cosTheta) >= 0.99) return false;
var vXS = this.vTemp2;
vXS.sub2 (this.ptTemp, this.p);
vXS.normalize ();
this.ptTemp.scaleAdd2 (dPS * cosTheta, vXS, this.p);
vXS.cross (v, vXS);
vXS.normalize ();
vXS.scale ((Math.sqrt (1 - cosTheta * cosTheta) * dPS));
this.ptS1.add2 (this.ptTemp, vXS);
this.ptS2.sub2 (this.ptTemp, vXS);
return true;
}, "J.jvxl.readers.IsoSolventReader.Edge,~N,~N,~N");
Clazz_defineMethod (c$, "checkSpecialVoxel", 
 function (ptA, ptB, ptV) {
var dAV = ptA.distance (ptV);
var dAV2 = ptA.distanceSquared (ptV);
var f = this.rAS / dAV;
if (f > 1) {
this.p.set (ptA.x + (ptV.x - ptA.x) * f, ptA.y + (ptV.y - ptA.y) * f, ptA.z + (ptV.z - ptA.z) * f);
return (ptB.distanceSquared (this.p) >= this.rBS2 ? NaN : this.solventDistance (this.rAS, this.rAS2, this.rBS2, dAV, dAV2, ptB.distanceSquared (ptV)));
}var dBV = ptB.distance (ptV);
if ((f = this.rBS / dBV) > 1) {
this.p.set (ptB.x + (ptV.x - ptB.x) * f, ptB.y + (ptV.y - ptB.y) * f, ptB.z + (ptV.z - ptB.z) * f);
return (ptA.distanceSquared (this.p) >= this.rAS2 ? NaN : this.solventDistance (this.rBS, this.rBS2, this.rAS2, dBV, dBV * dBV, dAV2));
}return NaN;
}, "JU.P3,JU.P3,JU.P3");
Clazz_defineMethod (c$, "solventDistance", 
 function (rAS, rAS2, rBS2, dAV, dAV2, dBV2) {
var angleVAB = Math.acos ((dAV2 + this.dAB2 - dBV2) / (2 * dAV * this.dAB));
var angleSAB = Math.acos ((rAS2 + this.dAB2 - rBS2) / (2 * rAS * this.dAB));
var dVS2 = (rAS2 + dAV2 - 2 * rAS * dAV * Math.cos (angleSAB - angleVAB));
var dVS = Math.sqrt (dVS2);
return (this.ecosASB2 < (rAS2 + dVS2 - dAV * dAV) / (dVS * rAS) ? dVS : NaN);
}, "~N,~N,~N,~N,~N,~N");
Clazz_defineMethod (c$, "dumpLine", 
function (pt1, pt2, label, color) {
this.sg.log ("draw ID \"x" + label + (this.nTest++) + "\" " + JU.P3.newP (pt1) + " " + JU.P3.newP (pt2) + " color " + color);
}, "JU.P3,JU.T3,~S,~S");
Clazz_defineMethod (c$, "dumpLine2", 
function (pt1, pt2, label, d, color1, color2) {
var pt =  new JU.V3 ();
pt.setT (pt2);
pt.sub (pt1);
pt.normalize ();
pt.scale (d);
pt.add (pt1);
this.sg.log ("draw ID \"" + label + (this.nTest++) + "\" " + JU.P3.newP (pt1) + " " + JU.P3.newP (pt) + " color " + color1);
this.sg.log ("draw ID \"" + label + (this.nTest++) + "\"" + JU.P3.newP (pt) + " " + JU.P3.newP (pt2) + " color " + color2 + "\"" + label + "\"");
}, "JU.P3,JU.P3,~S,~N,~S,~S");
Clazz_defineMethod (c$, "dumpPoint", 
function (pt, label, color) {
this.sg.log ("draw ID \"" + label + (this.nTest++) + "\"" + JU.P3.newP (pt) + " color " + color);
}, "JU.P3,~S,~S");
Clazz_overrideMethod (c$, "getValueAtPoint", 
function (pt, getSource) {
if (this.contactPair != null) return pt.distance (this.contactPair.myAtoms[1]) - this.contactPair.radii[1];
var value = 3.4028235E38;
for (var iAtom = 0; iAtom < this.firstNearbyAtom; iAtom++) {
if (this.rs == null || this.atomXyz == null || this.atomXyz[iAtom] == null || pt == null) System.out.println ("HOH");
var r = pt.distance (this.atomXyz[iAtom]) - this.rs[iAtom];
if (r < value) value = r;
}
return (value == 3.4028235E38 ? NaN : value);
}, "JU.T3,~B");
Clazz_overrideMethod (c$, "discardTempData", 
function (discardAll) {
this.rs = null;
this.rs2 = null;
this.discardTempDataSR (discardAll);
}, "~B");
Clazz_overrideMethod (c$, "getPlane", 
function (x) {
if (this.yzCount == 0) {
this.initPlanes ();
}this.thisX = x;
this.thisPlane = this.yzPlanes[x % 2];
if (this.contactPair == null) {
this.resetPlane (3.4028235E38);
this.thisAtomSet = this.bsAtomMinMax[x];
this.markSphereVoxels (0, this.params.distance);
this.unsetVoxelData ();
} else {
this.markPlaneVoxels (this.contactPair.myAtoms[0], this.contactPair.radii[0]);
}return this.thisPlane;
}, "~N");
c$.$IsoSolventReader$Edge$ = function () {
Clazz_pu$h(self.c$);
c$ = Clazz_decorateAsClass (function () {
Clazz_prepareCallback (this, arguments);
this.ia = 0;
this.ib = 0;
this.nFaces = 0;
this.nInvalid = 0;
this.d = 0;
this.d2 = 0;
this.maxr = 0;
this.cosASB2 = 0;
this.v = null;
Clazz_instantialize (this, arguments);
}, J.jvxl.readers.IsoSolventReader, "Edge", JU.P3);
Clazz_overrideConstructor (c$, 
function (a, b, c, d) {
this.ia = Math.min (b, c);
this.ib = Math.max (b, c);
this.d = d;
this.d2 = d * d;
this.maxr = Math.sqrt (this.d2 / 4 + Math.max (a.rs2[b], a.rs2[c]));
this.ave (a.atomXyz[b], a.atomXyz[c]);
this.cosASB2 = (a.rs2[b] + a.rs2[c] - this.d2) / (a.rs[c] * a.rs[b]);
this.v = JU.V3.newVsub (a.atomXyz[c], a.atomXyz[b]);
this.v.normalize ();
}, "J.jvxl.readers.IsoSolventReader,~N,~N,~N");
Clazz_defineMethod (c$, "addFace", 
function (a) {
this.nFaces++;
if (a == null) {
this.nInvalid++;
return;
}}, "J.jvxl.readers.IsoSolventReader.Face");
Clazz_defineMethod (c$, "isValid", 
function () {
return (this.nFaces == 0 || this.nInvalid != this.nFaces);
});
Clazz_overrideMethod (c$, "toString", 
function () {
return this.ia + "_" + this.ib;
});
c$ = Clazz_p0p ();
};
c$.$IsoSolventReader$Face$ = function () {
Clazz_pu$h(self.c$);
c$ = Clazz_decorateAsClass (function () {
Clazz_prepareCallback (this, arguments);
this.ia = 0;
this.ib = 0;
this.ic = 0;
this.pS = null;
Clazz_instantialize (this, arguments);
}, J.jvxl.readers.IsoSolventReader, "Face");
Clazz_makeConstructor (c$, 
function (a, b, c, d) {
this.ia = a;
this.ib = b;
this.ic = c;
this.pS = JU.P3.newP (d);
}, "~N,~N,~N,JU.P3");
Clazz_overrideMethod (c$, "toString", 
function () {
return this.ia + "_" + this.ib + "_" + this.ic + "_" + this.pS;
});
c$ = Clazz_p0p ();
};
Clazz_defineStatics (c$,
"testLinear", false);
});
Clazz_declarePackage ("J.jvxl.readers");
Clazz_load (["J.jvxl.readers.SurfaceReader"], "J.jvxl.readers.SurfaceFileReader", ["JU.PT", "J.api.Interface"], function () {
c$ = Clazz_decorateAsClass (function () {
this.br = null;
this.binarydoc = null;
this.out = null;
this.line = null;
this.next = null;
Clazz_instantialize (this, arguments);
}, J.jvxl.readers, "SurfaceFileReader", J.jvxl.readers.SurfaceReader);
Clazz_prepareFields (c$, function () {
this.next =  Clazz_newIntArray (1, 0);
});
Clazz_makeConstructor (c$, 
function () {
Clazz_superConstructor (this, J.jvxl.readers.SurfaceFileReader, []);
});
Clazz_overrideMethod (c$, "init", 
function (sg) {
this.initSR (sg);
}, "J.jvxl.readers.SurfaceGenerator");
Clazz_defineMethod (c$, "init2", 
function (sg, br) {
this.init2SFR (sg, br);
}, "J.jvxl.readers.SurfaceGenerator,java.io.BufferedReader");
Clazz_defineMethod (c$, "init2SFR", 
function (sg, br) {
this.init (sg);
this.br = br;
}, "J.jvxl.readers.SurfaceGenerator,java.io.BufferedReader");
Clazz_defineMethod (c$, "newBinaryDocument", 
function () {
return J.api.Interface.getInterface ("JU.BinaryDocument");
});
Clazz_overrideMethod (c$, "setOutputChannel", 
function (out) {
if (this.binarydoc == null) this.out = out;
 else this.sg.setOutputChannel (this.binarydoc, out);
}, "JU.OC");
Clazz_overrideMethod (c$, "closeReader", 
function () {
this.closeReaderSFR ();
});
Clazz_defineMethod (c$, "closeReaderSFR", 
function () {
if (this.br != null) try {
this.br.close ();
} catch (e) {
if (Clazz_exceptionOf (e, java.io.IOException)) {
} else {
throw e;
}
}
if (this.out != null) this.out.closeChannel ();
if (this.binarydoc != null) this.binarydoc.close ();
});
Clazz_overrideMethod (c$, "discardTempData", 
function (discardAll) {
this.closeReader ();
this.discardTempDataSR (discardAll);
}, "~B");
Clazz_defineMethod (c$, "getTokens", 
function () {
return JU.PT.getTokensAt (this.line, 0);
});
Clazz_defineMethod (c$, "parseFloat", 
function () {
return JU.PT.parseFloatNext (this.line, this.next);
});
Clazz_defineMethod (c$, "parseFloatStr", 
function (s) {
this.next[0] = 0;
return JU.PT.parseFloatNext (s, this.next);
}, "~S");
Clazz_defineMethod (c$, "parseFloatRange", 
function (s, iStart, iEnd) {
this.next[0] = iStart;
return JU.PT.parseFloatRange (s, iEnd, this.next);
}, "~S,~N,~N");
Clazz_defineMethod (c$, "parseInt", 
function () {
return JU.PT.parseIntNext (this.line, this.next);
});
Clazz_defineMethod (c$, "parseIntStr", 
function (s) {
this.next[0] = 0;
return JU.PT.parseIntNext (s, this.next);
}, "~S");
Clazz_defineMethod (c$, "parseIntNext", 
function (s) {
return JU.PT.parseIntNext (s, this.next);
}, "~S");
Clazz_defineMethod (c$, "parseFloatArrayStr", 
function (s) {
this.next[0] = 0;
return JU.PT.parseFloatArrayNext (s, this.next, null, null, null);
}, "~S");
Clazz_defineMethod (c$, "parseFloatArray", 
function (a, strStart, strEnd) {
return JU.PT.parseFloatArrayNext (this.line, this.next, a, strStart, strEnd);
}, "~A,~S,~S");
Clazz_defineMethod (c$, "getQuotedStringNext", 
function () {
return JU.PT.getQuotedStringNext (this.line, this.next);
});
Clazz_defineMethod (c$, "skipTo", 
function (info, what) {
if (info != null) while (this.readLine ().indexOf (info) < 0) {
}
if (what != null) this.next[0] = this.line.indexOf (what) + what.length + 2;
}, "~S,~S");
Clazz_defineMethod (c$, "readLine", 
function () {
this.line = this.br.readLine ();
if (this.line != null) {
this.nBytes += this.line.length;
if (this.out != null) {
var b = this.line.getBytes ();
this.out.write (b, 0, b.length);
this.out.writeByteAsInt (0x0A);
}}return this.line;
});
});
Clazz_declarePackage ("J.jvxl.readers");
Clazz_load (["J.jvxl.readers.SurfaceFileReader"], "J.jvxl.readers.VolumeFileReader", ["java.lang.Float", "JU.AU", "$.PT", "$.SB", "J.api.Interface", "J.atomdata.AtomData", "JU.Logger"], function () {
c$ = Clazz_decorateAsClass (function () {
this.endOfData = false;
this.negativeAtomCount = false;
this.ac = 0;
this.nSurfaces = 0;
this.isAngstroms = false;
this.canDownsample = false;
this.downsampleRemainders = null;
this.preProcessPlanes = false;
this.nData = 0;
this.readerClosed = false;
this.downsampleFactor = 0;
this.nSkipX = 0;
this.nSkipY = 0;
this.nSkipZ = 0;
this.yzPlanesRaw = null;
this.iPlaneRaw = 0;
this.boundingBox = null;
this.isScaledAlready = false;
Clazz_instantialize (this, arguments);
}, J.jvxl.readers, "VolumeFileReader", J.jvxl.readers.SurfaceFileReader);
Clazz_makeConstructor (c$, 
function () {
Clazz_superConstructor (this, J.jvxl.readers.VolumeFileReader, []);
});
Clazz_overrideMethod (c$, "init2", 
function (sg, br) {
this.init2VFR (sg, br);
}, "J.jvxl.readers.SurfaceGenerator,java.io.BufferedReader");
Clazz_defineMethod (c$, "init2VFR", 
function (sg, br) {
this.init2SFR (sg, br);
this.canDownsample = this.isProgressive = this.isXLowToHigh = true;
this.jvxlData.wasCubic = true;
this.boundingBox = this.params.boundingBox;
if (this.params.qmOrbitalType == 4) {
this.hasColorData = (this.params.parameters == null || this.params.parameters[1] >= 0);
this.preProcessPlanes = true;
this.params.insideOut = !this.params.insideOut;
}}, "J.jvxl.readers.SurfaceGenerator,java.io.BufferedReader");
Clazz_defineMethod (c$, "recordData", 
function (value) {
if (Float.isNaN (value)) return value;
if (value < this.dataMin) this.dataMin = value;
if (value > this.dataMax) this.dataMax = value;
this.dataMean += value;
this.nData++;
return value;
}, "~N");
Clazz_overrideMethod (c$, "closeReader", 
function () {
if (this.readerClosed) return;
this.readerClosed = true;
this.closeReaderSFR ();
if (this.nData == 0 || this.dataMax == -3.4028235E38) return;
this.dataMean /= this.nData;
JU.Logger.info ("VolumeFileReader closing file: " + this.nData + " points read \ndata min/max/mean = " + this.dataMin + "/" + this.dataMax + "/" + this.dataMean);
});
Clazz_overrideMethod (c$, "readVolumeParameters", 
function (isMapData) {
this.endOfData = false;
this.nSurfaces = this.readVolumetricHeader ();
if (this.nSurfaces == 0) return false;
if (this.nSurfaces < this.params.fileIndex) {
JU.Logger.warn ("not enough surfaces in file -- resetting params.fileIndex to " + this.nSurfaces);
this.params.fileIndex = this.nSurfaces;
}return true;
}, "~B");
Clazz_overrideMethod (c$, "readVolumeData", 
function (isMapData) {
return this.readVolumeDataVFR (isMapData);
}, "~B");
Clazz_defineMethod (c$, "readVolumeDataVFR", 
function (isMapData) {
if (!this.gotoAndReadVoxelData (isMapData)) return false;
if (!this.vertexDataOnly) JU.Logger.info ("JVXL read: " + this.nPointsX + " x " + this.nPointsY + " x " + this.nPointsZ + " data points");
return true;
}, "~B");
Clazz_defineMethod (c$, "readVolumetricHeader", 
 function () {
try {
this.readParameters ();
if (this.ac == -2147483648) return 0;
if (!this.vertexDataOnly) JU.Logger.info ("voxel grid origin:" + this.volumetricOrigin);
var downsampleFactor = this.params.downsampleFactor;
var downsampling = (this.canDownsample && downsampleFactor > 1);
if (downsampling) {
this.downsampleRemainders =  Clazz_newIntArray (3, 0);
JU.Logger.info ("downsample factor = " + downsampleFactor);
for (var i = 0; i < 3; ++i) {
var n = this.voxelCounts[i];
this.downsampleRemainders[i] = n % downsampleFactor;
this.voxelCounts[i] /= downsampleFactor;
this.volumetricVectors[i].scale (downsampleFactor);
JU.Logger.info ("downsampling axis " + (i + 1) + " from " + n + " to " + this.voxelCounts[i]);
}
}if (!this.vertexDataOnly) for (var i = 0; i < 3; ++i) {
if (!this.isAngstroms) this.volumetricVectors[i].scale (0.5291772);
this.line = this.voxelCounts[i] + " " + this.volumetricVectors[i].x + " " + this.volumetricVectors[i].y + " " + this.volumetricVectors[i].z;
this.jvxlFileHeaderBuffer.append (this.line).appendC ('\n');
JU.Logger.info ("voxel grid count/vector:" + this.line);
}
this.scaleIsosurface (this.params.scale);
this.volumeData.setVolumetricXml ();
return this.nSurfaces;
} catch (e) {
if (Clazz_exceptionOf (e, Exception)) {
JU.Logger.error (e.toString ());
return 0;
} else {
throw e;
}
}
});
Clazz_defineMethod (c$, "skipComments", 
function (allowBlankLines) {
var sb =  new JU.SB ();
while (this.readLine () != null && (allowBlankLines && this.line.length == 0 || this.line.indexOf ("#") == 0)) sb.append (this.line).appendC ('\n');

return sb.toString ();
}, "~B");
Clazz_defineMethod (c$, "readVoxelVector", 
function (voxelVectorIndex) {
this.readLine ();
var voxelVector = this.volumetricVectors[voxelVectorIndex];
if ((this.voxelCounts[voxelVectorIndex] = this.parseIntStr (this.line)) == -2147483648) this.next[0] = this.line.indexOf (" ");
voxelVector.set (this.parseFloat (), this.parseFloat (), this.parseFloat ());
if (this.isAnisotropic) this.setVectorAnisotropy (voxelVector);
}, "~N");
Clazz_defineMethod (c$, "initializeSurfaceData", 
function () {
this.downsampleFactor = this.params.downsampleFactor;
this.nSkipX = 0;
this.nSkipY = 0;
this.nSkipZ = 0;
if (this.canDownsample && this.downsampleFactor > 0) {
this.nSkipX = this.downsampleFactor - 1;
this.nSkipY = this.downsampleRemainders[2] + (this.downsampleFactor - 1) * (this.nSkipZ = (this.nPointsZ * this.downsampleFactor + this.downsampleRemainders[2]));
this.nSkipZ = this.downsampleRemainders[1] * this.nSkipZ + (this.downsampleFactor - 1) * this.nSkipZ * (this.nPointsY * this.downsampleFactor + this.downsampleRemainders[1]);
}if (this.params.thePlane != null) {
this.params.cutoff = 0;
} else if (this.isJvxl) {
this.params.cutoff = (this.params.isBicolorMap || this.params.colorBySign ? 0.01 : 0.5);
}this.nDataPoints = 0;
this.next[0] = 0;
this.line = "";
this.jvxlNSurfaceInts = 0;
});
Clazz_overrideMethod (c$, "readSurfaceData", 
function (isMapData) {
this.readSurfaceDataVFR (isMapData);
}, "~B");
Clazz_defineMethod (c$, "readSurfaceDataVFR", 
function (isMapData) {
this.initializeSurfaceData ();
if (this.isProgressive && !isMapData || this.isJvxl) {
this.nDataPoints = this.volumeData.setVoxelCounts (this.nPointsX, this.nPointsY, this.nPointsZ);
this.voxelData = null;
if (this.isJvxl) this.jvxlVoxelBitSet = this.getVoxelBitSet (this.nDataPoints);
} else if (isMapData && this.volumeData.hasPlane ()) {
this.volumeData.setVoxelMap ();
var f = this.volumeData.getToPlaneParameter ();
for (var x = 0; x < this.nPointsX; ++x) {
for (var y = 0; y < this.nPointsY; ++y) {
for (var z = 0; z < this.nPointsZ; ++z) {
var v = this.recordData (this.getNextVoxelValue ());
if (this.volumeData.isNearPlane (x, y, z, f)) this.volumeData.setVoxelMapValue (x, y, z, v);
if (this.nSkipX != 0) this.skipVoxels (this.nSkipX);
}
if (this.nSkipY != 0) this.skipVoxels (this.nSkipY);
}
if (this.nSkipZ != 0) this.skipVoxels (this.nSkipZ);
}
} else {
this.voxelData = JU.AU.newFloat3 (this.nPointsX, -1);
for (var x = 0; x < this.nPointsX; ++x) {
var plane = JU.AU.newFloat2 (this.nPointsY);
this.voxelData[x] = plane;
for (var y = 0; y < this.nPointsY; ++y) {
var strip =  Clazz_newFloatArray (this.nPointsZ, 0);
plane[y] = strip;
for (var z = 0; z < this.nPointsZ; ++z) {
strip[z] = this.recordData (this.getNextVoxelValue ());
if (this.nSkipX != 0) this.skipVoxels (this.nSkipX);
}
if (this.nSkipY != 0) this.skipVoxels (this.nSkipY);
}
if (this.nSkipZ != 0) this.skipVoxels (this.nSkipZ);
}
}this.volumeData.setVoxelDataAsArray (this.voxelData);
}, "~B");
Clazz_defineMethod (c$, "getPlane", 
function (x) {
if (x == 0) this.initPlanes ();
if (this.preProcessPlanes) return this.getPlaneProcessed (x);
var plane = this.getPlane2 (x);
if (this.qpc == null) this.getPlane (plane, true);
return plane;
}, "~N");
Clazz_defineMethod (c$, "getPlaneProcessed", 
function (x) {
var plane;
if (this.iPlaneRaw == 0) {
this.qpc = J.api.Interface.getOption ("quantum.NciCalculation");
var atomData =  new J.atomdata.AtomData ();
atomData.modelIndex = -1;
atomData.bsSelected = this.params.bsSelected;
this.sg.fillAtomData (atomData, 1);
this.qpc.setupCalculation (this.volumeData, this.sg.getBsSelected (), null, null, null, atomData.atomXyz, -1, null, null, null, null, null, null, this.params.isSquaredLinear, null, this.params.theProperty, true, null, this.params.parameters, this.params.testFlags);
this.iPlaneRaw = 1;
this.qpc.setPlanes (this.yzPlanesRaw =  Clazz_newFloatArray (4, this.yzCount, 0));
if (this.hasColorData) {
this.getPlane (this.yzPlanesRaw[0], false);
this.getPlane (this.yzPlanesRaw[1], false);
plane = this.yzPlanes[0];
for (var i = 0; i < this.yzCount; i++) plane[i] = NaN;

return plane;
}this.iPlaneRaw = -1;
}var nan = this.qpc.getNoValue ();
var x1 = this.nPointsX - 1;
switch (this.iPlaneRaw) {
case -1:
plane = this.yzPlanes[x % 2];
x1++;
break;
case 3:
plane = this.yzPlanesRaw[0];
this.yzPlanesRaw[0] = this.yzPlanesRaw[1];
this.yzPlanesRaw[1] = this.yzPlanesRaw[2];
this.yzPlanesRaw[2] = this.yzPlanesRaw[3];
this.yzPlanesRaw[3] = plane;
plane = this.yzPlanesRaw[this.iPlaneRaw];
break;
default:
this.iPlaneRaw++;
plane = this.yzPlanesRaw[this.iPlaneRaw];
}
if (x < x1) {
this.getPlane (plane, false);
this.qpc.calcPlane (x, plane = this.yzPlanes[x % 2]);
for (var i = 0; i < this.yzCount; i++) if (plane[i] != nan) this.recordData (plane[i]);

} else {
for (var i = 0; i < this.yzCount; i++) plane[i] = NaN;

}return plane;
}, "~N");
Clazz_defineMethod (c$, "getPlane", 
 function (plane, doRecord) {
try {
for (var y = 0, ptyz = 0; y < this.nPointsY; ++y) {
for (var z = 0; z < this.nPointsZ; ++z) {
var v = this.getNextVoxelValue ();
if (doRecord) this.recordData (v);
plane[ptyz++] = v;
if (this.nSkipX != 0) this.skipVoxels (this.nSkipX);
}
if (this.nSkipY != 0) this.skipVoxels (this.nSkipY);
}
if (this.nSkipZ != 0) this.skipVoxels (this.nSkipZ);
} catch (e) {
if (Clazz_exceptionOf (e, Exception)) {
} else {
throw e;
}
}
}, "~A,~B");
Clazz_overrideMethod (c$, "getValue", 
function (x, y, z, ptyz) {
if (this.boundingBox != null) {
this.volumeData.voxelPtToXYZ (x, y, z, this.ptTemp);
if (this.ptTemp.x < this.boundingBox[0].x || this.ptTemp.x > this.boundingBox[1].x || this.ptTemp.y < this.boundingBox[0].y || this.ptTemp.y > this.boundingBox[1].y || this.ptTemp.z < this.boundingBox[0].z || this.ptTemp.z > this.boundingBox[1].z) return NaN;
}return this.getValue2 (x, y, z, ptyz);
}, "~N,~N,~N,~N");
Clazz_defineMethod (c$, "skipVoxels", 
 function (n) {
for (var i = n; --i >= 0; ) this.getNextVoxelValue ();

}, "~N");
Clazz_defineMethod (c$, "getVoxelBitSet", 
function (nPoints) {
return null;
}, "~N");
Clazz_defineMethod (c$, "getNextVoxelValue", 
function () {
var voxelValue = 0;
if (this.nSurfaces > 1 && !this.params.blockCubeData) {
for (var i = 1; i < this.params.fileIndex; i++) this.nextVoxel ();

voxelValue = this.nextVoxel ();
for (var i = this.params.fileIndex; i < this.nSurfaces; i++) this.nextVoxel ();

} else {
voxelValue = this.nextVoxel ();
}return voxelValue;
});
Clazz_defineMethod (c$, "nextVoxel", 
function () {
var voxelValue = this.parseFloat ();
if (Float.isNaN (voxelValue)) {
while (this.readLine () != null && Float.isNaN (voxelValue = this.parseFloatStr (this.line))) {
}
if (this.line == null) {
if (!this.endOfData) JU.Logger.warn ("end of file reading cube voxel data? nBytes=" + this.nBytes + " nDataPoints=" + this.nDataPoints + " (line):" + this.line);
this.endOfData = true;
this.line = "0 0 0 0 0 0 0 0 0 0";
}}return voxelValue;
});
Clazz_overrideMethod (c$, "gotoData", 
function (n, nPoints) {
if (!this.params.blockCubeData) return;
if (n > 0) JU.Logger.info ("skipping " + n + " data sets, " + nPoints + " points each");
for (var i = 0; i < n; i++) this.skipData (nPoints);

}, "~N,~N");
Clazz_defineMethod (c$, "skipData", 
function (nPoints) {
this.skipDataVFR (nPoints);
}, "~N");
Clazz_defineMethod (c$, "skipDataVFR", 
function (nPoints) {
var iV = 0;
while (iV < nPoints) iV += this.countData (this.readLine ());

}, "~N");
Clazz_defineMethod (c$, "countData", 
 function (str) {
var count = 0;
var ich = 0;
var ichMax = str.length;
var ch;
while (ich < ichMax) {
while (ich < ichMax && ((ch = str.charAt (ich)) == ' ' || ch == '\t')) ++ich;

if (ich < ichMax) ++count;
while (ich < ichMax && ((ch = str.charAt (ich)) != ' ' && ch != '\t')) ++ich;

}
return count;
}, "~S");
c$.checkAtomLine = Clazz_defineMethod (c$, "checkAtomLine", 
function (isXLowToHigh, isAngstroms, strAtomCount, atomLine, bs) {
if (atomLine.indexOf ("ANGSTROMS") >= 0) isAngstroms = true;
var ac = (strAtomCount == null ? 2147483647 : JU.PT.parseInt (strAtomCount));
switch (ac) {
case -2147483648:
ac = 0;
atomLine = " " + atomLine.substring (atomLine.indexOf (" ") + 1);
break;
case 2147483647:
ac = -2147483648;
break;
default:
var s = "" + ac;
atomLine = atomLine.substring (atomLine.indexOf (s) + s.length);
}
if (isAngstroms) {
if (atomLine.indexOf ("ANGSTROM") < 0) atomLine += " ANGSTROMS";
} else {
if (atomLine.indexOf ("BOHR") < 0) atomLine += " BOHR";
}atomLine = (ac == -2147483648 ? "" : (isXLowToHigh ? "+" : "-") + Math.abs (ac)) + atomLine + "\n";
bs.append (atomLine);
return isAngstroms;
}, "~B,~B,~S,~S,JU.SB");
Clazz_overrideMethod (c$, "getSurfacePointAndFraction", 
function (cutoff, isCutoffAbsolute, valueA, valueB, pointA, edgeVector, x, y, z, vA, vB, fReturn, ptReturn) {
return this.getSPFv (cutoff, isCutoffAbsolute, valueA, valueB, pointA, edgeVector, x, y, z, vA, vB, fReturn, ptReturn);
}, "~N,~B,~N,~N,JU.T3,JU.V3,~N,~N,~N,~N,~N,~A,JU.T3");
Clazz_defineMethod (c$, "getSPFv", 
function (cutoff, isCutoffAbsolute, valueA, valueB, pointA, edgeVector, x, y, z, vA, vB, fReturn, ptReturn) {
var zero = this.getSPF (cutoff, isCutoffAbsolute, valueA, valueB, pointA, edgeVector, x, y, z, vA, vB, fReturn, ptReturn);
if (this.qpc == null || Float.isNaN (zero) || !this.hasColorData) return zero;
vA = this.marchingCubes.getLinearOffset (x, y, z, vA);
vB = this.marchingCubes.getLinearOffset (x, y, z, vB);
return this.qpc.process (vA, vB, fReturn[0]);
}, "~N,~B,~N,~N,JU.T3,JU.V3,~N,~N,~N,~N,~N,~A,JU.T3");
Clazz_defineMethod (c$, "scaleIsosurface", 
 function (scale) {
if (this.isScaledAlready) return;
this.isScaledAlready = true;
if (this.isAnisotropic) this.setVolumetricAnisotropy ();
if (Float.isNaN (scale)) return;
JU.Logger.info ("applying scaling factor of " + scale);
this.volumetricOrigin.scaleAdd2 ((1 - scale) / 2, this.volumetricVectors[0], this.volumetricOrigin);
this.volumetricOrigin.scaleAdd2 ((1 - scale) / 2, this.volumetricVectors[1], this.volumetricOrigin);
this.volumetricOrigin.scaleAdd2 ((1 - scale) / 2, this.volumetricVectors[2], this.volumetricOrigin);
this.volumetricVectors[0].scale (scale);
this.volumetricVectors[1].scale (scale);
this.volumetricVectors[2].scale (scale);
}, "~N");
});
Clazz_declarePackage ("J.jvxl.readers");
Clazz_load (["J.jvxl.readers.VolumeFileReader"], "J.jvxl.readers.JvxlXmlReader", ["java.lang.Float", "$.NullPointerException", "java.util.Hashtable", "JU.AU", "$.BS", "$.CU", "$.Lst", "$.P3", "$.P4", "$.PT", "$.SB", "J.jvxl.data.JvxlCoder", "$.MeshData", "J.jvxl.readers.XmlReader", "J.shapesurface.IsosurfaceMesh", "JU.C", "$.ColorEncoder", "$.Escape", "$.Logger"], function () {
c$ = Clazz_decorateAsClass (function () {
this.JVXL_VERSION = "2.3";
this.surfaceDataCount = 0;
this.edgeDataCount = 0;
this.colorDataCount = 0;
this.excludedTriangleCount = 0;
this.excludedVertexCount = 0;
this.invalidatedVertexCount = 0;
this.haveContourData = false;
this.xr = null;
this.isXmlFile = true;
this.thisInside = false;
this.tempDataXml = null;
this.bsVoxelBitSet = null;
this.includeValueNaN = true;
this.valueCount = 0;
this.valueMin = NaN;
this.valueRange = NaN;
this.fractionPtr = 0;
this.colorPtr = 0;
this.strFractionTemp = "";
this.haveReadColorData = false;
this.jvxlColorEncodingRead = null;
Clazz_instantialize (this, arguments);
}, J.jvxl.readers, "JvxlXmlReader", J.jvxl.readers.VolumeFileReader);
Clazz_makeConstructor (c$, 
function () {
Clazz_superConstructor (this, J.jvxl.readers.JvxlXmlReader, []);
});
Clazz_overrideMethod (c$, "init2", 
function (sg, br) {
this.init2JXR (sg, br);
}, "J.jvxl.readers.SurfaceGenerator,java.io.BufferedReader");
Clazz_defineMethod (c$, "init2JXR", 
function (sg, br) {
this.init2VFR (sg, br);
this.jvxlData.wasJvxl = this.isJvxl = true;
this.isXLowToHigh = false;
this.xr =  new J.jvxl.readers.XmlReader (br);
}, "J.jvxl.readers.SurfaceGenerator,java.io.BufferedReader");
Clazz_overrideMethod (c$, "readVolumeData", 
function (isMapData) {
if (!this.readVolumeDataVFR (isMapData)) return false;
this.strFractionTemp = this.jvxlEdgeDataRead;
this.fractionPtr = 0;
return true;
}, "~B");
Clazz_overrideMethod (c$, "gotoAndReadVoxelData", 
function (isMapData) {
this.initializeVolumetricData ();
if (this.nPointsX < 0 || this.nPointsY < 0 || this.nPointsZ < 0) return true;
try {
this.gotoData (this.params.fileIndex - 1, this.nPointsX * this.nPointsY * this.nPointsZ);
if (this.vertexDataOnly) return true;
this.volumeData.setMappingPlane (this.params.thePlane);
this.readSurfaceData (isMapData);
this.volumeData.setMappingPlane (null);
if (this.edgeDataCount > 0) this.jvxlEdgeDataRead = this.jvxlReadFractionData ("edge", this.edgeDataCount);
this.params.bsExcluded = this.jvxlData.jvxlExcluded =  new Array (4);
this.hasColorData = (this.colorDataCount > 0);
if (this.hasColorData) this.jvxlColorDataRead = this.jvxlReadFractionData ("color", this.colorDataCount);
if (this.excludedVertexCount > 0) {
this.jvxlData.jvxlExcluded[0] = J.jvxl.data.JvxlCoder.jvxlDecodeBitSet (this.xr.getXmlData ("jvxlExcludedVertexData", null, false, false));
if (this.xr.isNext ("jvxlExcludedPlaneData")) this.jvxlData.jvxlExcluded[2] = J.jvxl.data.JvxlCoder.jvxlDecodeBitSet (this.xr.getXmlData ("jvxlExcludedPlaneData", null, false, false));
}if (this.excludedTriangleCount > 0) this.jvxlData.jvxlExcluded[3] = J.jvxl.data.JvxlCoder.jvxlDecodeBitSet (this.xr.getXmlData ("jvxlExcludedTriangleData", null, false, false));
if (this.invalidatedVertexCount > 0) this.jvxlData.jvxlExcluded[1] = J.jvxl.data.JvxlCoder.jvxlDecodeBitSet (this.xr.getXmlData ("jvxlInvalidatedVertexData", null, false, false));
if (this.haveContourData) this.jvxlDecodeContourData (this.jvxlData, this.xr.getXmlData ("jvxlContourData", null, false, false));
if (this.jvxlDataIsColorMapped && this.jvxlData.nVertexColors > 0) {
this.jvxlData.vertexColorMap =  new java.util.Hashtable ();
var vdata = this.xr.getXmlData ("jvxlVertexColorData", null, true, false);
var baseColor = J.jvxl.readers.XmlReader.getXmlAttrib (vdata, "baseColor");
this.jvxlData.baseColor = (baseColor.length > 0 ? baseColor : null);
for (var i = 0; i < this.jvxlData.nVertexColors; i++) {
var s = this.xr.getXmlData ("jvxlColorMap", vdata, true, false);
var color = J.jvxl.readers.XmlReader.getXmlAttrib (s, "color");
var bs = J.jvxl.data.JvxlCoder.jvxlDecodeBitSet (this.xr.getXmlData ("jvxlColorMap", s, false, false));
this.jvxlData.vertexColorMap.put (color, bs);
}
}} catch (e) {
if (Clazz_exceptionOf (e, Exception)) {
JU.Logger.error (e.toString ());
return false;
} else {
throw e;
}
}
return true;
}, "~B");
Clazz_overrideMethod (c$, "readParameters", 
function () {
var s = this.xr.getXmlData ("jvxlFileTitle", null, false, false);
this.jvxlFileHeaderBuffer = JU.SB.newS (s);
this.xr.toTag ("jvxlVolumeData");
var data = this.tempDataXml = this.xr.getXmlData ("jvxlVolumeData", null, true, false);
this.volumetricOrigin.setT (this.xr.getXmlPoint (data, "origin"));
this.isAngstroms = true;
this.readVector (0);
this.readVector (1);
this.readVector (2);
this.line = this.xr.toTag ("jvxlSurfaceSet");
this.nSurfaces = this.parseIntStr (J.jvxl.readers.XmlReader.getXmlAttrib (this.line, "count"));
JU.Logger.info ("jvxl file surfaces: " + this.nSurfaces);
JU.Logger.info ("using default edge fraction base and range");
JU.Logger.info ("using default color fraction base and range");
this.cJvxlEdgeNaN = String.fromCharCode (this.edgeFractionBase + this.edgeFractionRange);
});
Clazz_defineMethod (c$, "readVector", 
function (voxelVectorIndex) {
var data = this.xr.getXmlData ("jvxlVolumeVector", this.tempDataXml, true, true);
this.tempDataXml = this.tempDataXml.substring (this.tempDataXml.indexOf (data) + data.length);
var n = this.parseIntStr (J.jvxl.readers.XmlReader.getXmlAttrib (data, "count"));
if (n == -2147483648) this.vertexDataOnly = true;
this.voxelCounts[voxelVectorIndex] = (n < 0 ? 0 : n);
this.volumetricVectors[voxelVectorIndex].setT (this.xr.getXmlPoint (data, "vector"));
if (this.isAnisotropic) this.setVectorAnisotropy (this.volumetricVectors[voxelVectorIndex]);
}, "~N");
Clazz_overrideMethod (c$, "gotoData", 
function (n, nPoints) {
if (n > 0) JU.Logger.info ("skipping " + n + " data sets, " + nPoints + " points each");
this.vertexDataOnly = this.jvxlData.vertexDataOnly = (nPoints == 0);
for (var i = 0; i < n; i++) {
this.jvxlSkipData (nPoints, true);
}
this.xr.toTag ("jvxlSurface");
this.jvxlReadSurfaceInfo ();
}, "~N,~N");
Clazz_defineMethod (c$, "jvxlSkipData", 
function (nPoints, doSkipColorData) {
this.readLine ();
this.xr.skipTag ("jvxlSurface");
}, "~N,~B");
Clazz_defineMethod (c$, "jvxlReadSurfaceInfo", 
function () {
var s;
var data = this.xr.getXmlData ("jvxlSurfaceInfo", null, true, true);
this.isXLowToHigh = J.jvxl.readers.XmlReader.getXmlAttrib (data, "isXLowToHigh").equals ("true");
this.jvxlCutoff = this.parseFloatStr (J.jvxl.readers.XmlReader.getXmlAttrib (data, "cutoff"));
if (!Float.isNaN (this.jvxlCutoff)) JU.Logger.info ("JVXL read: cutoff " + this.jvxlCutoff);
var nContourData = this.parseIntStr (J.jvxl.readers.XmlReader.getXmlAttrib (data, "nContourData"));
this.haveContourData = (nContourData > 0);
this.params.isContoured = J.jvxl.readers.XmlReader.getXmlAttrib (data, "contoured").equals ("true");
if (this.params.isContoured) {
var nContoursRead = this.parseIntStr (J.jvxl.readers.XmlReader.getXmlAttrib (data, "nContours"));
if (nContoursRead <= 0) {
nContoursRead = 0;
} else {
s = J.jvxl.readers.XmlReader.getXmlAttrib (data, "contourValues");
if (s.length > 0) {
this.jvxlData.contourValues = this.params.contoursDiscrete = this.parseFloatArrayStr (s);
JU.Logger.info ("JVXL read: contourValues " + JU.Escape.eAF (this.jvxlData.contourValues));
}s = J.jvxl.readers.XmlReader.getXmlAttrib (data, "contourColors");
if (s.length > 0) {
this.jvxlData.contourColixes = this.params.contourColixes = JU.C.getColixArray (s);
this.jvxlData.contourColors = JU.C.getHexCodes (this.jvxlData.contourColixes);
JU.Logger.info ("JVXL read: contourColixes " + JU.C.getHexCodes (this.jvxlData.contourColixes));
}this.params.contourFromZero = J.jvxl.readers.XmlReader.getXmlAttrib (data, "contourFromZero").equals ("true");
}this.params.nContours = (this.haveContourData ? nContourData : nContoursRead);
}this.jvxlData.nVertexColors = this.parseIntStr (J.jvxl.readers.XmlReader.getXmlAttrib (data, "nVertexColors"));
this.params.isBicolorMap = J.jvxl.readers.XmlReader.getXmlAttrib (data, "bicolorMap").equals ("true");
if (this.params.isBicolorMap) {
s = J.jvxl.readers.XmlReader.getXmlAttrib (data, "colorPositive");
if (s.length > 0 && this.params.colorRgb == -2147483648 && this.params.colorPos == -16776961) this.params.colorPos = JU.CU.getArgbFromString (s);
s = J.jvxl.readers.XmlReader.getXmlAttrib (data, "colorNegative");
if (s.length > 0 && this.params.colorRgb == -2147483648 && this.params.colorNeg == -65536) this.params.colorNeg = JU.CU.getArgbFromString (s);
}if (this.params.isBicolorMap || this.params.colorBySign) this.jvxlCutoff = 0;
this.jvxlDataIsColorMapped = ((this.params.colorRgb == -2147483648 || this.params.colorRgb == 2147483647) && (this.params.isBicolorMap || J.jvxl.readers.XmlReader.getXmlAttrib (data, "colorMapped").equals ("true")));
this.jvxlData.isJvxlPrecisionColor = J.jvxl.readers.XmlReader.getXmlAttrib (data, "precisionColor").equals ("true");
this.jvxlData.jvxlDataIsColorDensity = this.params.colorDensity = (this.params.colorRgb == -2147483648 && J.jvxl.readers.XmlReader.getXmlAttrib (data, "colorDensity").equals ("true"));
if (this.jvxlData.jvxlDataIsColorDensity && Float.isNaN (this.params.pointSize)) {
s = J.jvxl.readers.XmlReader.getXmlAttrib (data, "pointSize");
if (s.length > 0) this.jvxlData.pointSize = this.params.pointSize = this.parseFloatStr (s);
}s = J.jvxl.readers.XmlReader.getXmlAttrib (data, "allowVolumeRender");
this.jvxlData.allowVolumeRender = this.params.allowVolumeRender = (s.length == 0 || s.equalsIgnoreCase ("true"));
s = J.jvxl.readers.XmlReader.getXmlAttrib (data, "plane");
if (s.indexOf ("{") >= 0) {
this.params.thePlane = null;
this.params.mapLattice = null;
try {
this.params.thePlane = JU.Escape.uP (s);
s = J.jvxl.readers.XmlReader.getXmlAttrib (data, "maplattice");
JU.Logger.info ("JVXL read: plane " + this.params.thePlane);
if (s.indexOf ("{") >= 0) {
this.params.mapLattice = JU.Escape.uP (s);
JU.Logger.info ("JVXL read: mapLattice " + this.params.mapLattice);
}if (this.params.scale3d == 0) this.params.scale3d = this.parseFloatStr (J.jvxl.readers.XmlReader.getXmlAttrib (data, "scale3d"));
if (Float.isNaN (this.params.scale3d)) this.params.scale3d = 0;
} catch (e) {
if (Clazz_exceptionOf (e, Exception)) {
if (this.params.thePlane == null) {
JU.Logger.error ("JVXL Error reading plane definition -- setting to 0 0 1 0  (z=0)");
this.params.thePlane = JU.P4.new4 (0, 0, 1, 0);
} else {
JU.Logger.error ("JVXL Error reading mapLattice definition -- ignored");
}} else {
throw e;
}
}
this.surfaceDataCount = 0;
this.edgeDataCount = 0;
} else {
this.params.thePlane = null;
this.surfaceDataCount = this.parseIntStr (J.jvxl.readers.XmlReader.getXmlAttrib (data, "nSurfaceInts"));
this.edgeDataCount = this.parseIntStr (J.jvxl.readers.XmlReader.getXmlAttrib (data, "nBytesUncompressedEdgeData"));
}this.excludedVertexCount = this.parseIntStr (J.jvxl.readers.XmlReader.getXmlAttrib (data, "nExcludedVertexes"));
this.excludedTriangleCount = this.parseIntStr (J.jvxl.readers.XmlReader.getXmlAttrib (data, "nExcludedTriangles"));
this.invalidatedVertexCount = this.parseIntStr (J.jvxl.readers.XmlReader.getXmlAttrib (data, "nInvalidatedVertexes"));
s = J.jvxl.readers.XmlReader.getXmlAttrib (data, "slabInfo");
if (s.length > 0) this.jvxlData.slabInfo = s;
this.colorDataCount = Math.max (0, this.parseIntStr (J.jvxl.readers.XmlReader.getXmlAttrib (data, "nBytesUncompressedColorData")));
this.jvxlDataIs2dContour = (this.params.thePlane != null && this.jvxlDataIsColorMapped);
this.jvxlData.color = J.jvxl.readers.XmlReader.getXmlAttrib (data, "color");
if (this.jvxlData.color.length == 0 || this.jvxlData.color.indexOf ("null") >= 0) this.jvxlData.color = "orange";
this.jvxlData.translucency = this.parseFloatStr (J.jvxl.readers.XmlReader.getXmlAttrib (data, "translucency"));
if (Float.isNaN (this.jvxlData.translucency)) this.jvxlData.translucency = 0;
s = J.jvxl.readers.XmlReader.getXmlAttrib (data, "meshColor");
if (s.length > 0) this.jvxlData.meshColor = s;
s = J.jvxl.readers.XmlReader.getXmlAttrib (data, "rendering");
if (s.length > 0) this.jvxlData.rendering = s;
this.jvxlData.colorScheme = J.jvxl.readers.XmlReader.getXmlAttrib (data, "colorScheme");
if (this.jvxlData.colorScheme.length == 0) this.jvxlData.colorScheme = null;
if (this.jvxlData.thisSet < 0) {
var n = this.parseIntStr (J.jvxl.readers.XmlReader.getXmlAttrib (data, "set"));
if (n > 0) this.jvxlData.thisSet = n - 1;
}this.jvxlData.slabValue = this.parseIntStr (J.jvxl.readers.XmlReader.getXmlAttrib (data, "slabValue"));
this.jvxlData.isSlabbable = (J.jvxl.readers.XmlReader.getXmlAttrib (data, "slabbable").equalsIgnoreCase ("true"));
this.jvxlData.diameter = this.parseIntStr (J.jvxl.readers.XmlReader.getXmlAttrib (data, "diameter"));
if (this.jvxlData.diameter == -2147483648) this.jvxlData.diameter = 0;
if (this.jvxlDataIs2dContour) this.params.isContoured = true;
if (this.params.colorBySign) this.params.isBicolorMap = true;
var insideOut = J.jvxl.readers.XmlReader.getXmlAttrib (data, "insideOut").equals ("true");
var dataMin = NaN;
var dataMax = NaN;
var red = NaN;
var blue = NaN;
if (this.jvxlDataIsColorMapped) {
dataMin = this.parseFloatStr (J.jvxl.readers.XmlReader.getXmlAttrib (data, "dataMinimum"));
dataMax = this.parseFloatStr (J.jvxl.readers.XmlReader.getXmlAttrib (data, "dataMaximum"));
red = this.parseFloatStr (J.jvxl.readers.XmlReader.getXmlAttrib (data, "valueMappedToRed"));
blue = this.parseFloatStr (J.jvxl.readers.XmlReader.getXmlAttrib (data, "valueMappedToBlue"));
if (Float.isNaN (dataMin)) {
dataMin = red = -1.0;
dataMax = blue = 1;
}}this.jvxlSetColorRanges (dataMin, dataMax, red, blue, insideOut);
});
Clazz_defineMethod (c$, "jvxlSetColorRanges", 
function (dataMin, dataMax, red, blue, insideOut) {
if (this.jvxlDataIsColorMapped) {
if (!Float.isNaN (dataMin) && !Float.isNaN (dataMax)) {
if (dataMax == 0 && dataMin == 0) {
dataMin = -1;
dataMax = 1;
}this.params.mappedDataMin = dataMin;
this.params.mappedDataMax = dataMax;
JU.Logger.info ("JVXL read: data_min/max " + this.params.mappedDataMin + "/" + this.params.mappedDataMax);
}if (!this.params.rangeDefined) if (!Float.isNaN (red) && !Float.isNaN (blue)) {
if (red == 0 && blue == 0) {
red = -1;
blue = 1;
}this.params.valueMappedToRed = Math.min (red, blue);
this.params.valueMappedToBlue = Math.max (red, blue);
this.params.isColorReversed = (red > blue);
this.params.rangeDefined = true;
} else {
this.params.valueMappedToRed = 0;
this.params.valueMappedToBlue = 1;
this.params.rangeDefined = true;
}JU.Logger.info ("JVXL read: color red/blue: " + this.params.valueMappedToRed + "/" + this.params.valueMappedToBlue);
}this.jvxlData.valueMappedToRed = this.params.valueMappedToRed;
this.jvxlData.valueMappedToBlue = this.params.valueMappedToBlue;
this.jvxlData.mappedDataMin = this.params.mappedDataMin;
this.jvxlData.mappedDataMax = this.params.mappedDataMax;
this.jvxlData.isColorReversed = this.params.isColorReversed;
this.jvxlData.insideOut = insideOut;
if (this.params.insideOut) this.jvxlData.insideOut = !this.jvxlData.insideOut;
this.params.insideOut = this.jvxlData.insideOut;
}, "~N,~N,~N,~N,~B");
Clazz_overrideMethod (c$, "readSurfaceData", 
function (isMapDataIgnored) {
this.thisInside = !this.params.isContoured;
if (this.readSurfaceDataXML ()) return;
this.tempDataXml = this.xr.getXmlData ("jvxlEdgeData", null, true, false);
this.bsVoxelBitSet = J.jvxl.data.JvxlCoder.jvxlDecodeBitSet (this.xr.getXmlData ("jvxlEdgeData", this.tempDataXml, false, false));
this.readSurfaceDataJXR ();
}, "~B");
Clazz_defineMethod (c$, "readSurfaceDataXML", 
function () {
if (this.vertexDataOnly) {
this.getEncodedVertexData ();
return true;
}if (this.params.thePlane != null) {
this.volumeData.setDataDistanceToPlane (this.params.thePlane);
this.setVolumeDataV (this.volumeData);
this.params.cutoff = 0;
this.jvxlData.setSurfaceInfo (this.params.thePlane, this.params.mapLattice, 0, "");
this.jvxlData.scale3d = this.params.scale3d;
return true;
}return false;
});
Clazz_defineMethod (c$, "readSurfaceDataJXR", 
function () {
this.readSurfaceDataVFR (false);
this.volumeData.setMappingPlane (null);
});
Clazz_defineMethod (c$, "jvxlReadFractionData", 
function (type, nPoints) {
var str;
try {
if (type.equals ("edge")) {
str = J.jvxl.data.JvxlCoder.jvxlDecompressString (J.jvxl.readers.XmlReader.getXmlAttrib (this.tempDataXml, "data"));
} else {
var data = this.xr.getXmlData ("jvxlColorData", null, true, false);
this.jvxlData.isJvxlPrecisionColor = J.jvxl.readers.JvxlXmlReader.getEncoding (data).endsWith ("2");
str = J.jvxl.data.JvxlCoder.jvxlDecompressString (J.jvxl.readers.XmlReader.getXmlAttrib (data, "data"));
}} catch (e) {
if (Clazz_exceptionOf (e, Exception)) {
JU.Logger.error ("Error reading " + type + " data " + e);
throw  new NullPointerException ();
} else {
throw e;
}
}
return str;
}, "~S,~N");
Clazz_overrideMethod (c$, "getVoxelBitSet", 
function (nPoints) {
if (this.bsVoxelBitSet != null) return this.bsVoxelBitSet;
var bs =  new JU.BS ();
var bsVoxelPtr = 0;
if (this.surfaceDataCount <= 0) return bs;
var nThisValue = 0;
while (bsVoxelPtr < nPoints) {
nThisValue = this.parseInt ();
if (nThisValue == -2147483648) {
this.readLine ();
if (this.line == null || (nThisValue = this.parseIntStr (this.line)) == -2147483648) {
if (!this.endOfData) JU.Logger.error ("end of file in JvxlReader?" + " line=" + this.line);
this.endOfData = true;
nThisValue = 10000;
}}this.thisInside = !this.thisInside;
++this.jvxlNSurfaceInts;
if (this.thisInside) bs.setBits (bsVoxelPtr, bsVoxelPtr + nThisValue);
bsVoxelPtr += nThisValue;
}
return bs;
}, "~N");
Clazz_overrideMethod (c$, "getSurfacePointAndFraction", 
function (cutoff, isCutoffAbsolute, valueA, valueB, pointA, edgeVector, x, y, z, vA, vB, fReturn, ptReturn) {
if (this.edgeDataCount <= 0) return this.getSPFv (cutoff, isCutoffAbsolute, valueA, valueB, pointA, edgeVector, x, y, z, vA, vB, fReturn, ptReturn);
ptReturn.scaleAdd2 (fReturn[0] = this.jvxlGetNextFraction (this.edgeFractionBase, this.edgeFractionRange, 0.5), edgeVector, pointA);
if (Float.isNaN (this.valueMin)) this.setValueMinMax ();
return (this.valueCount == 0 || this.includeValueNaN && Float.isNaN (fReturn[0]) ? fReturn[0] : this.getNextValue ());
}, "~N,~B,~N,~N,JU.T3,JU.V3,~N,~N,~N,~N,~N,~A,JU.T3");
Clazz_defineMethod (c$, "getNextValue", 
 function () {
var fraction = NaN;
while (this.colorPtr < this.valueCount && Float.isNaN (fraction)) {
if (this.jvxlData.isJvxlPrecisionColor) {
fraction = J.jvxl.data.JvxlCoder.jvxlFractionFromCharacter2 (this.jvxlColorDataRead.charCodeAt (this.colorPtr), this.jvxlColorDataRead.charCodeAt ((this.colorPtr++) + this.valueCount), this.colorFractionBase, this.colorFractionRange);
} else {
fraction = J.jvxl.data.JvxlCoder.jvxlFractionFromCharacter (this.jvxlColorDataRead.charCodeAt (this.colorPtr++), this.colorFractionBase, this.colorFractionRange, 0.5);
}break;
}
return this.valueMin + fraction * this.valueRange;
});
Clazz_defineMethod (c$, "setValueMinMax", 
 function () {
this.valueCount = this.jvxlColorDataRead.length;
if (this.jvxlData.isJvxlPrecisionColor) this.valueCount /= 2;
this.includeValueNaN = (this.valueCount != this.jvxlEdgeDataRead.length);
this.valueMin = (!this.jvxlData.isJvxlPrecisionColor ? this.params.valueMappedToRed : this.params.mappedDataMin == 3.4028235E38 ? 0.0 : this.params.mappedDataMin);
this.valueRange = (!this.jvxlData.isJvxlPrecisionColor ? this.params.valueMappedToBlue : this.params.mappedDataMin == 3.4028235E38 ? 1.0 : this.params.mappedDataMax) - this.valueMin;
this.haveReadColorData = true;
});
Clazz_defineMethod (c$, "jvxlGetNextFraction", 
 function (base, range, fracOffset) {
if (this.fractionPtr >= this.strFractionTemp.length) {
if (!this.endOfData) JU.Logger.error ("end of file reading compressed fraction data");
this.endOfData = true;
this.strFractionTemp = "" + String.fromCharCode (base);
this.fractionPtr = 0;
}return J.jvxl.data.JvxlCoder.jvxlFractionFromCharacter (this.strFractionTemp.charCodeAt (this.fractionPtr++), base, range, fracOffset);
}, "~N,~N,~N");
Clazz_overrideMethod (c$, "readColorData", 
function () {
if (!this.jvxlDataIsColorMapped) return "";
var vertexCount = this.jvxlData.vertexCount = this.meshData.vc;
var colixes = this.meshData.vcs;
var vertexValues = this.meshData.vvs;
if ("none".equals (this.jvxlColorEncodingRead)) {
this.jvxlData.vertexColors =  Clazz_newIntArray (vertexCount, 0);
var nextc =  Clazz_newIntArray (1, 0);
var n = JU.PT.parseIntNext (this.jvxlColorDataRead, nextc);
n = Math.min (n, vertexCount);
var tokens = JU.PT.getTokens (this.jvxlColorDataRead.substring (nextc[0]));
var haveTranslucent = false;
var trans = this.jvxlData.translucency;
var lastColor = 0;
for (var i = 0; i < n; i++) try {
var c = J.jvxl.readers.JvxlXmlReader.getColor (tokens[i]);
if (c == 0) c = lastColor;
 else lastColor = c;
colixes[i] = JU.C.getColixTranslucent (this.jvxlData.vertexColors[i] = c);
if (JU.C.isColixTranslucent (colixes[i])) haveTranslucent = true;
 else if (trans != 0) colixes[i] = JU.C.getColixTranslucent3 (colixes[i], true, trans);
} catch (e) {
if (Clazz_exceptionOf (e, Exception)) {
JU.Logger.info ("JvxlXmlReader: Cannot interpret color code: " + tokens[i]);
} else {
throw e;
}
}

if (haveTranslucent && trans == 0) {
this.jvxlData.translucency = 0.5;
}return "-";
}if (this.params.colorEncoder == null) this.params.colorEncoder =  new JU.ColorEncoder (null);
this.params.colorEncoder.setColorScheme (null, false);
this.params.colorEncoder.setRange (this.params.valueMappedToRed, this.params.valueMappedToBlue, this.params.isColorReversed);
JU.Logger.info ("JVXL reading color data mapped min/max: " + this.params.mappedDataMin + "/" + this.params.mappedDataMax + " for " + vertexCount + " vertices." + " using encoding keys " + this.colorFractionBase + " " + this.colorFractionRange);
JU.Logger.info ("mapping red-->blue for " + this.params.valueMappedToRed + " to " + this.params.valueMappedToBlue + " colorPrecision:" + this.jvxlData.isJvxlPrecisionColor);
var getValues = (Float.isNaN (this.valueMin));
if (getValues) this.setValueMinMax ();
var contourPlaneMinimumValue = 3.4028235E38;
var contourPlaneMaximumValue = -3.4028235E38;
if (colixes == null || colixes.length < vertexCount) this.meshData.vcs = colixes =  Clazz_newShortArray (vertexCount, 0);
var colixNeg = 0;
var colixPos = 0;
if (this.params.colorBySign) {
colixPos = JU.C.getColix (this.params.isColorReversed ? this.params.colorNeg : this.params.colorPos);
colixNeg = JU.C.getColix (this.params.isColorReversed ? this.params.colorPos : this.params.colorNeg);
}var vertexIncrement = this.meshData.vertexIncrement;
var needContourMinMax = (this.params.mappedDataMin == 3.4028235E38);
for (var i = 0; i < vertexCount; i += vertexIncrement) {
var value;
if (getValues) value = vertexValues[i] = this.getNextValue ();
 else value = vertexValues[i];
if (needContourMinMax) {
if (value < contourPlaneMinimumValue) contourPlaneMinimumValue = value;
if (value > contourPlaneMaximumValue) contourPlaneMaximumValue = value;
}}
if (needContourMinMax) {
this.params.mappedDataMin = contourPlaneMinimumValue;
this.params.mappedDataMax = contourPlaneMaximumValue;
}if (this.jvxlData.colorScheme != null) for (var i = 0; i < vertexCount; i += vertexIncrement) {
var value = vertexValues[i];
if (this.marchingSquares != null && this.params.isContoured) {
this.marchingSquares.setContourData (i, value);
continue;
}var colix = (!this.params.colorBySign ? this.params.colorEncoder.getColorIndex (value) : (this.params.isColorReversed ? value > 0 : value <= 0) ? colixNeg : colixPos);
colixes[i] = JU.C.getColixTranslucent3 (colix, true, this.jvxlData.translucency);
}
return this.jvxlColorDataRead + "\n";
});
c$.getColor = Clazz_defineMethod (c$, "getColor", 
 function (c) {
var n = 0;
try {
switch (c.charAt (0)) {
case '[':
n = JU.CU.getArgbFromString (c);
break;
case '0':
n = JU.PT.parseIntRadix (c.substring (2), 16);
break;
default:
n = JU.PT.parseIntRadix (c, 10);
}
} catch (e) {
if (Clazz_exceptionOf (e, Exception)) {
} else {
throw e;
}
}
return n;
}, "~S");
Clazz_defineMethod (c$, "getEncodedVertexData", 
function () {
var sdata = this.xr.getXmlData ("jvxlSurfaceData", null, true, false);
this.jvxlDecodeVertexData (this.xr.getXmlData ("jvxlVertexData", sdata, true, false), false);
var tData = this.xr.getXmlData ("jvxlTriangleData", sdata, true, false);
var edgeData = this.xr.getXmlData ("jvxlTriangleEdgeData", sdata, true, false);
var polygonColorData = this.xr.getXmlData ("jvxlPolygonColorData", sdata, false, false);
this.jvxlDecodeTriangleData (tData, edgeData, polygonColorData);
var cData = this.xr.getXmlData ("jvxlColorData", sdata, true, false);
this.jvxlColorEncodingRead = J.jvxl.readers.JvxlXmlReader.getEncoding (cData);
this.jvxlData.isJvxlPrecisionColor = this.jvxlColorEncodingRead.endsWith ("2");
cData = this.getData (cData, "jvxlColorData");
this.jvxlColorDataRead = (this.jvxlColorEncodingRead.equals ("none") ? cData : J.jvxl.data.JvxlCoder.jvxlDecompressString (cData));
this.jvxlDataIsColorMapped = ((this.params.colorRgb == -2147483648 || this.params.colorRgb == 2147483647) && this.jvxlColorDataRead.length > 0);
if (this.haveContourData) this.jvxlDecodeContourData (this.jvxlData, this.xr.getXmlData ("jvxlContourData", null, false, false));
});
Clazz_defineMethod (c$, "getData", 
 function (sdata, name) {
var data = J.jvxl.readers.XmlReader.getXmlAttrib (sdata, "data");
if (data.length == 0) data = this.xr.getXmlData (name, sdata, false, false);
return data;
}, "~S,~S");
c$.getEncoding = Clazz_defineMethod (c$, "getEncoding", 
 function (data) {
if (J.jvxl.readers.XmlReader.getXmlAttrib (data, "len").length > 0) return "";
var s = J.jvxl.readers.XmlReader.getXmlAttrib (data, "encoding");
return (s.length == 0 ? "none" : s);
}, "~S");
Clazz_defineMethod (c$, "jvxlDecodeVertexData", 
function (data, asArray) {
var vertexCount = this.parseIntStr (J.jvxl.readers.XmlReader.getXmlAttrib (data, "count"));
if (!asArray) JU.Logger.info ("Reading " + vertexCount + " vertices");
var ptCount = vertexCount * 3;
var vertices = (asArray ?  new Array (vertexCount) : null);
var fraction;
var vData = J.jvxl.readers.XmlReader.getXmlAttrib (data, "data");
var encoding = J.jvxl.readers.JvxlXmlReader.getEncoding (data);
if ("none".equals (encoding)) {
if (vData.length == 0) vData = this.xr.getXmlData ("jvxlVertexData", data, false, false);
var fdata = JU.PT.parseFloatArray (vData);
if (fdata[0] != vertexCount * 3) JU.Logger.info ("JvxlXmlReader: vertexData count=" + (Clazz_floatToInt (fdata[0])) + "; expected " + (vertexCount * 3));
for (var i = 0, pt = 1; i < vertexCount; i++) {
var p = JU.P3.new3 (fdata[pt++], fdata[pt++], fdata[pt++]);
if (asArray) vertices[i] = p;
 else this.addVertexCopy (p, 0, i, false);
}
} else {
var min = this.xr.getXmlPoint (data, "min");
var range = this.xr.getXmlPoint (data, "max");
range.sub (min);
var colorFractionBase = this.jvxlData.colorFractionBase;
var colorFractionRange = this.jvxlData.colorFractionRange;
var s = J.jvxl.data.JvxlCoder.jvxlDecompressString (vData);
if (s.length == 0) s = this.xr.getXmlData ("jvxlVertexData", data, false, false);
for (var i = 0, pt = -1; i < vertexCount; i++) {
var p =  new JU.P3 ();
fraction = J.jvxl.data.JvxlCoder.jvxlFractionFromCharacter2 (s.charCodeAt (++pt), s.charCodeAt (pt + ptCount), colorFractionBase, colorFractionRange);
p.x = min.x + fraction * range.x;
fraction = J.jvxl.data.JvxlCoder.jvxlFractionFromCharacter2 (s.charCodeAt (++pt), s.charCodeAt (pt + ptCount), colorFractionBase, colorFractionRange);
p.y = min.y + fraction * range.y;
fraction = J.jvxl.data.JvxlCoder.jvxlFractionFromCharacter2 (s.charCodeAt (++pt), s.charCodeAt (pt + ptCount), colorFractionBase, colorFractionRange);
p.z = min.z + fraction * range.z;
if (asArray) vertices[i] = p;
 else this.addVertexCopy (p, 0, i, false);
}
}return vertices;
}, "~S,~B");
Clazz_defineMethod (c$, "jvxlDecodeTriangleData", 
function (tdata, edgeData, colorData) {
var nTriangles = this.parseIntStr (J.jvxl.readers.XmlReader.getXmlAttrib (tdata, "count"));
if (nTriangles < 0) return;
var nextc =  Clazz_newIntArray (1, 0);
var nColors = (colorData == null ? -1 : JU.PT.parseIntNext (colorData, nextc));
var color = 0;
JU.Logger.info ("Reading " + nTriangles + " triangles");
var encoding = J.jvxl.readers.JvxlXmlReader.getEncoding (tdata);
tdata = this.getData (tdata, "jvxlTriangleData");
var edata = this.getData (edgeData, "jvxlTriangleEdgeData");
var vertex =  Clazz_newIntArray (3, 0);
var nextp =  Clazz_newIntArray (1, 0);
var nexte = null;
var edgeMask = 7;
var haveEdgeInfo;
var haveEncoding = !"none".equals (encoding);
if (haveEncoding) {
tdata = J.jvxl.data.JvxlCoder.jvxlDecompressString (tdata);
edata = J.jvxl.data.JvxlCoder.jvxlDecompressString (edata).trim ();
haveEdgeInfo = (edata.length == nTriangles);
} else {
var n = JU.PT.parseIntNext (tdata, nextp);
haveEdgeInfo = (edata.length > 0);
if (haveEdgeInfo) {
nexte =  Clazz_newIntArray (1, 0);
JU.PT.parseIntNext (edata, nexte);
} else if (n > 0) {
JU.Logger.info ("JvxlXmlReader: jvxlTriangleEdgeData count=" + n + "; expected " + nTriangles);
}}for (var i = 0, v = 0, p = 0, pt = -1; i < nTriangles; ) {
if (haveEncoding) {
var ch = tdata.charAt (++pt);
var diff;
switch (ch) {
case '!':
diff = 0;
break;
case '+':
case '.':
case ' ':
case '\n':
case '\r':
case '\t':
case ',':
continue;
case '-':
case '0':
case '1':
case '2':
case '3':
case '4':
case '5':
case '6':
case '7':
case '8':
case '9':
nextp[0] = pt;
diff = JU.PT.parseIntNext (tdata, nextp);
pt = nextp[0] - 1;
break;
default:
diff = ch.charCodeAt (0) - 92;
}
v += diff;
} else {
v = JU.PT.parseIntNext (tdata, nextp) - 1;
}vertex[p] = v;
if (++p == 3) {
p = 0;
if (haveEdgeInfo) {
edgeMask = (nexte == null ? edata.charCodeAt (i) - 48 : JU.PT.parseIntNext (edata, nexte));
if (edgeMask < 0 || edgeMask > 7) edgeMask = 7;
}if (nColors > 0) {
var c = JU.PT.parseIntNext (colorData, nextc);
if (c == -2147483648) nColors = 0;
 else color = c;
nColors--;
}this.addTriangleCheck (vertex[0], vertex[1], vertex[2], edgeMask, 0, false, color);
i++;
}}
}, "~S,~S,~S");
Clazz_defineMethod (c$, "jvxlDecodeContourData", 
function (jvxlData, data) {
var vs =  new JU.Lst ();
var values =  new JU.SB ();
var colors =  new JU.SB ();
var pt = -1;
jvxlData.vContours = null;
if (data == null) return;
while ((pt = data.indexOf ("<jvxlContour", pt + 1)) >= 0) {
var v =  new JU.Lst ();
var s = this.xr.getXmlData ("jvxlContour", data.substring (pt), true, false);
var value = this.parseFloatStr (J.jvxl.readers.XmlReader.getXmlAttrib (s, "value"));
values.append (" ").appendF (value);
var color = J.jvxl.readers.JvxlXmlReader.getColor (J.jvxl.readers.XmlReader.getXmlAttrib (s, "color"));
var colix = JU.C.getColix (color);
colors.append (" ").append (JU.Escape.escapeColor (color));
var fData = J.jvxl.data.JvxlCoder.jvxlDecompressString (J.jvxl.readers.XmlReader.getXmlAttrib (s, "data"));
var bs = J.jvxl.data.JvxlCoder.jvxlDecodeBitSet (this.xr.getXmlData ("jvxlContour", s, false, false));
var n = bs.length ();
J.shapesurface.IsosurfaceMesh.setContourVector (v, n, bs, value, colix, color, JU.SB.newS (fData));
vs.addLast (v);
}
var n = vs.size ();
if (n > 0) {
jvxlData.vContours = JU.AU.createArrayOfArrayList (n);
jvxlData.contourColixes = this.params.contourColixes =  Clazz_newShortArray (n, 0);
jvxlData.contourValues = this.params.contoursDiscrete =  Clazz_newFloatArray (n, 0);
for (var i = 0; i < n; i++) {
jvxlData.vContours[i] = vs.get (i);
jvxlData.contourValues[i] = (jvxlData.vContours[i].get (2)).floatValue ();
jvxlData.contourColixes[i] = (jvxlData.vContours[i].get (3))[0];
}
jvxlData.contourColors = JU.C.getHexCodes (jvxlData.contourColixes);
JU.Logger.info ("JVXL read: " + n + " discrete contours");
JU.Logger.info ("JVXL read: contour values: " + values);
JU.Logger.info ("JVXL read: contour colors: " + colors);
}}, "J.jvxl.data.JvxlData,~S");
Clazz_overrideMethod (c$, "postProcessVertices", 
function () {
var bsInvalid = this.params.bsExcluded[1];
if (bsInvalid != null) {
if (this.meshDataServer != null) this.meshDataServer.fillMeshData (this.meshData, 1, null);
this.meshData.invalidateVertices (bsInvalid);
if (this.meshDataServer != null) {
this.meshDataServer.fillMeshData (this.meshData, 4, null);
this.meshData =  new J.jvxl.data.MeshData ();
}this.updateTriangles ();
}});
});
Clazz_declarePackage ("J.rendersurface");
Clazz_load (["J.render.MeshRenderer", "JU.P3", "$.P3i"], "J.rendersurface.IsosurfaceRenderer", ["java.lang.Boolean", "$.Float", "JU.V3", "JU.C", "$.Normix"], function () {
c$ = Clazz_decorateAsClass (function () {
this.iHideBackground = false;
this.isBicolorMap = false;
this.backgroundColix = 0;
this.nError = 0;
this.vertexValues = null;
this.imesh = null;
this.isosurface = null;
this.isNavigationMode = false;
this.iShowNormals = false;
this.showNumbers = false;
this.$showKey = null;
this.hasColorRange = false;
this.meshScale = -1;
this.mySlabValue = 0;
this.globalSlabValue = 0;
this.ptTemp = null;
this.ptTempi = null;
Clazz_instantialize (this, arguments);
}, J.rendersurface, "IsosurfaceRenderer", J.render.MeshRenderer);
Clazz_prepareFields (c$, function () {
this.ptTemp =  new JU.P3 ();
this.ptTempi =  new JU.P3i ();
});
Clazz_overrideMethod (c$, "render", 
function () {
return this.renderIso ();
});
Clazz_defineMethod (c$, "renderIso", 
function () {
this.setGlobals ();
for (var i = this.isosurface.meshCount; --i >= 0; ) {
this.mesh = this.imesh = this.isosurface.meshes[i];
if (this.imesh.connections != null && !this.vwr.ms.at[this.imesh.connections[0]].checkVisible ()) continue;
this.hasColorRange = false;
if (this.renderMeshSlab ()) {
if (!this.isExport) this.renderInfo ();
if (this.isExport && this.isGhostPass) {
this.exportPass = 1;
this.renderMeshSlab ();
this.exportPass = 2;
}}}
return this.needTranslucent;
});
Clazz_defineMethod (c$, "setGlobals", 
 function () {
this.needTranslucent = false;
this.iShowNormals = this.vwr.getTestFlag (4);
this.showNumbers = this.vwr.getTestFlag (3);
this.isosurface = this.shape;
this.exportPass = (this.isExport ? 2 : 0);
this.isNavigationMode = this.vwr.getBoolean (603979887);
this.$showKey = (this.vwr.getBoolean (603979870) ? Boolean.TRUE : null);
this.isosurface.keyXy = null;
this.meshScale = -1;
this.globalSlabValue = this.g3d.getSlab ();
this.mySlabValue = (this.isNavigationMode ? Clazz_floatToInt (this.tm.getNavigationOffset ().z) : 2147483647);
});
Clazz_defineMethod (c$, "renderInfo", 
function () {
if (this.hasColorRange && this.imesh.colorEncoder != null && Boolean.TRUE === this.$showKey) this.showKey ();
});
Clazz_defineMethod (c$, "showKey", 
 function () {
this.$showKey = Boolean.FALSE;
var colors = null;
var colixes = null;
var vContours = null;
var n = 0;
var type = 0;
if (this.imesh.showContourLines) {
vContours = this.imesh.getContours ();
if (vContours == null) {
colixes = this.imesh.jvxlData.contourColixes;
if (colixes == null) return;
n = colixes.length;
} else {
n = vContours.length;
type = 1;
}} else {
colors = this.imesh.colorEncoder.getColorSchemeArray (this.imesh.colorEncoder.currentPalette);
n = (colors == null ? 0 : colors.length);
type = 2;
}if (n < 2) return;
var factor = (this.g3d.isAntialiased () ? 2 : 1);
var height = this.vwr.getScreenHeight () * factor;
var dy = Clazz_doubleToInt (Clazz_doubleToInt (height / 2) / (n - 1));
var y = Clazz_doubleToInt (height / 4) * 3 - dy;
var x = 10 * factor;
var dx = 20 * factor;
this.isosurface.keyXy = [Clazz_doubleToInt (x / factor), 0, Clazz_doubleToInt ((x + dx) / factor), Clazz_doubleToInt ((y + dy) / factor), Clazz_doubleToInt (dy / factor)];
for (var i = 0; i < n; i++, y -= dy) {
switch (type) {
case 0:
if (!this.g3d.setC (colixes[i])) return;
break;
case 1:
if (!this.g3d.setC ((vContours[i].get (3))[0])) return;
break;
case 2:
this.g3d.setColor (colors[i]);
break;
}
this.g3d.fillRect (x, y, 5, -2147483648, dx, dy);
}
this.isosurface.keyXy[1] = Clazz_doubleToInt ((y + dy) / factor);
});
Clazz_defineMethod (c$, "renderMeshSlab", 
 function () {
this.volumeRender = (this.imesh.jvxlData.colorDensity && this.imesh.jvxlData.allowVolumeRender);
var thisSlabValue = this.mySlabValue;
this.frontOnly = this.mesh.frontOnly || this.shapeID == 26;
if (!this.isNavigationMode) {
this.meshSlabValue = this.imesh.jvxlData.slabValue;
if (this.meshSlabValue != -2147483648 && this.imesh.jvxlData.isSlabbable) {
var points = this.imesh.jvxlData.boundingBox;
var z0 = 3.4028235E38;
var z1 = 1.4E-45;
for (var i = points.length; --i >= 0; ) {
this.pt2f.setT (points[i]);
this.tm.transformPt3f (this.pt2f, this.pt2f);
if (this.pt2f.z < z0) z0 = this.pt2f.z;
if (this.pt2f.z > z1) z1 = this.pt2f.z;
}
thisSlabValue = Math.round (z0 + (z1 - z0) * (100 - this.meshSlabValue) / 100);
this.frontOnly = new Boolean (this.frontOnly & (this.meshSlabValue >= 100)).valueOf ();
}}var tcover = this.g3d.getTranslucentCoverOnly ();
this.g3d.setTranslucentCoverOnly (this.frontOnly || !this.vwr.getBoolean (603979967));
this.thePlane = this.imesh.jvxlData.jvxlPlane;
this.vertexValues = this.mesh.vvs;
var isOK;
if (thisSlabValue != 2147483647 && this.imesh.jvxlData.isSlabbable) {
this.g3d.setSlab (thisSlabValue);
isOK = this.renderMesh (this.mesh);
this.g3d.setSlab (this.globalSlabValue);
} else {
isOK = this.renderMesh (this.mesh);
}this.g3d.setTranslucentCoverOnly (tcover);
return isOK;
});
Clazz_overrideMethod (c$, "render2", 
function (isExport) {
if (this.volumeRender) {
this.renderPoints ();
return;
}switch (this.imesh.dataType) {
case 70:
this.renderLonePair (false);
return;
case 71:
this.renderLonePair (true);
return;
}
this.isBicolorMap = this.imesh.jvxlData.isBicolorMap;
this.render2b (isExport);
if (!this.g3d.setC (4)) return;
if (this.imesh.showContourLines) this.renderContourLines ();
}, "~B");
Clazz_defineMethod (c$, "renderLonePair", 
 function (isRadical) {
this.pt2f.setT (this.vertices[1]);
this.tm.transformPt3f (this.pt2f, this.pt2f);
var r = Clazz_floatToInt (this.vwr.tm.scaleToScreen (Clazz_floatToInt (this.pt2f.z), 100));
if (r < 1) r = 1;
if (!isRadical) {
var v1 =  new JU.V3 ();
var v2 =  new JU.V3 ();
this.pt1f.setT (this.vertices[0]);
this.tm.transformPt3f (this.pt1f, this.pt1f);
v1.sub2 (this.pt2f, this.pt1f);
v2.set (v1.x, v1.y, v1.z + 1);
v2.cross (v2, v1);
v2.normalize ();
var f = this.vwr.tm.scaleToScreen (Clazz_floatToInt (this.pt1f.z), 100);
v2.scale (f);
this.pt1f.add2 (this.pt2f, v2);
this.pt2f.sub (v2);
this.screens[0].set (Math.round (this.pt1f.x), Math.round (this.pt1f.y), Math.round (this.pt1f.z));
this.g3d.fillSphereI (r, this.screens[0]);
}this.screens[1].set (Math.round (this.pt2f.x), Math.round (this.pt2f.y), Math.round (this.pt2f.z));
this.g3d.fillSphereI (r, this.screens[1]);
}, "~B");
Clazz_defineMethod (c$, "renderContourLines", 
 function () {
var vContours = this.imesh.getContours ();
if (vContours == null) {
if (this.imesh.jvxlData.contourValues != null) this.hasColorRange = true;
return;
}this.hasColorRange = (this.mesh.meshColix == 0);
for (var i = vContours.length; --i >= 0; ) {
var v = vContours[i];
if (v.size () < 6) continue;
this.colix = (this.mesh.meshColix == 0 ? (v.get (3))[0] : this.mesh.meshColix);
if (!this.g3d.setC (this.colix)) return;
var n = v.size () - 1;
for (var j = 6; j < n; j++) {
var pt1 = v.get (j);
var pt2 = v.get (++j);
this.tm.transformPtScr (pt1, this.pt1i);
this.tm.transformPtScr (pt2, this.pt2i);
if (Float.isNaN (pt1.x) || Float.isNaN (pt2.x)) break;
this.pt1i.z -= 2;
this.pt2i.z -= 2;
this.g3d.drawLineAB (this.pt1i, this.pt2i);
}
}
});
Clazz_overrideMethod (c$, "renderPoints", 
function () {
try {
if (this.volumeRender) this.g3d.volumeRender (true);
var slabPoints = ((this.volumeRender || this.mesh.pc == 0) && this.selectedPolyOnly);
var incr = this.imesh.vertexIncrement;
var diam;
if (this.mesh.diameter <= 0) {
diam = this.vwr.getInt (553648144);
this.frontOnly = false;
} else {
diam = Clazz_doubleToInt (this.vwr.getScreenDim () / (this.volumeRender ? 50 : 100));
}var ptSize = Math.round (Float.isNaN (this.mesh.volumeRenderPointSize) ? 150 : this.mesh.volumeRenderPointSize * 1000);
if (diam < 1) diam = 1;
var cX = (this.showNumbers ? Clazz_doubleToInt (this.vwr.getScreenWidth () / 2) : 0);
var cY = (this.showNumbers ? Clazz_doubleToInt (this.vwr.getScreenHeight () / 2) : 0);
if (this.showNumbers) this.g3d.setFontFid (this.g3d.getFontFidFS ("Monospaced", 24));
for (var i = (!this.imesh.hasGridPoints || this.imesh.firstRealVertex < 0 ? 0 : this.imesh.firstRealVertex); i < this.vertexCount; i += incr) {
if (this.vertexValues != null && Float.isNaN (this.vertexValues[i]) || this.frontOnly && this.transformedVectors[this.normixes[i]].z < 0 || this.imesh.jvxlData.thisSet >= 0 && this.mesh.vertexSets[i] != this.imesh.jvxlData.thisSet || !this.mesh.isColorSolid && this.mesh.vcs != null && !this.setColix (this.mesh.vcs[i]) || this.haveBsDisplay && !this.mesh.bsDisplay.get (i) || slabPoints && !this.bsPolygons.get (i)) continue;
this.hasColorRange = true;
if (this.showNumbers && this.screens[i].z > 10 && Math.abs (this.screens[i].x - cX) < 150 && Math.abs (this.screens[i].y - cY) < 150) {
var s = i + (this.mesh.isColorSolid ? "" : " " + this.mesh.vvs[i]);
this.g3d.setC (4);
this.g3d.drawStringNoSlab (s, null, this.screens[i].x, this.screens[i].y, this.screens[i].z - 30, 0);
}if (this.volumeRender) {
diam = Clazz_floatToInt (this.vwr.tm.scaleToScreen (this.screens[i].z, ptSize));
if (diam < 1) diam = 1;
this.g3d.volumeRender4 (diam, this.screens[i].x, this.screens[i].y, this.screens[i].z);
} else {
this.g3d.fillSphereI (diam, this.screens[i]);
}}
if (incr == 3) {
this.g3d.setC (this.isTranslucent ? JU.C.getColixTranslucent3 (12, true, 0.5) : 12);
for (var i = 1; i < this.vertexCount; i += 3) this.g3d.fillCylinder (3, Clazz_doubleToInt (diam / 4), this.screens[i], this.screens[i + 1]);

this.g3d.setC (this.isTranslucent ? JU.C.getColixTranslucent3 (21, true, 0.5) : 21);
for (var i = 1; i < this.vertexCount; i += 3) this.g3d.fillSphereI (diam, this.screens[i]);

this.g3d.setC (this.isTranslucent ? JU.C.getColixTranslucent3 (7, true, 0.5) : 7);
for (var i = 2; i < this.vertexCount; i += 3) {
this.g3d.fillSphereI (diam, this.screens[i]);
}
}} catch (e) {
}
if (this.volumeRender) this.g3d.volumeRender (false);
});
Clazz_overrideMethod (c$, "renderTriangles", 
function (fill, iShowTriangles, isExport) {
this.g3d.addRenderer (1073742182);
var polygonIndexes = this.mesh.pis;
this.colix = (this.isGhostPass ? this.mesh.slabColix : !fill && this.mesh.meshColix != 0 ? this.mesh.meshColix : this.mesh.colix);
var vertexColixes = (!fill && this.mesh.meshColix != 0 ? null : this.mesh.vcs);
if (this.isTranslucentInherit) this.colix = JU.C.copyColixTranslucency (this.mesh.slabColix, this.mesh.colix);
this.g3d.setC (this.colix);
var generateSet = isExport;
if (generateSet) {
if (this.frontOnly && fill) this.frontOnly = false;
this.bsPolygonsToExport.clearAll ();
}if (this.exportType == 1) {
this.frontOnly = false;
}var colorSolid = (this.isGhostPass && (!this.isBicolorMap) || vertexColixes == null || this.mesh.isColorSolid);
var noColor = (this.isGhostPass && !this.isBicolorMap || vertexColixes == null || !fill && this.mesh.meshColix != 0);
var isPlane = (this.imesh.jvxlData.jvxlPlane != null);
var colix = this.colix;
if (isPlane && !colorSolid && !fill && this.mesh.fillTriangles) {
colorSolid = true;
colix = 4;
}var colorArrayed = (colorSolid && this.mesh.pcs != null);
if (colorArrayed && !fill && this.mesh.fillTriangles) colorArrayed = false;
var contourColixes = this.imesh.jvxlData.contourColixes;
this.hasColorRange = !colorSolid && !this.isBicolorMap;
var diam;
if (this.mesh.diameter <= 0) {
diam = (this.meshScale < 0 ? this.meshScale = this.vwr.getInt (553648151) : this.meshScale);
if (this.g3d.isAntialiased ()) diam *= 2;
} else {
diam = Clazz_doubleToInt (this.vwr.getScreenDim () / 100);
}if (diam < 1) diam = 1;
for (var i = this.mesh.pc; --i >= 0; ) {
var polygon = polygonIndexes[i];
if (polygon == null || this.selectedPolyOnly && !this.bsPolygons.get (i)) continue;
var iA = polygon[0];
var iB = polygon[1];
var iC = polygon[2];
if (this.imesh.jvxlData.thisSet >= 0 && this.mesh.vertexSets != null && this.mesh.vertexSets[iA] != this.imesh.jvxlData.thisSet) continue;
if (this.haveBsDisplay && (!this.mesh.bsDisplay.get (iA) || !this.mesh.bsDisplay.get (iB) || !this.mesh.bsDisplay.get (iC))) continue;
var nA = this.normixes[iA];
var nB = this.normixes[iB];
var nC = this.normixes[iC];
var check = this.checkNormals (nA, nB, nC);
if (fill && check == 0) continue;
var colixA;
var colixB;
var colixC;
if (colorSolid) {
if (colorArrayed && i < this.mesh.pcs.length) {
var c = this.mesh.pcs[i];
if (c == 0) continue;
colix = c;
}colixA = colixB = colixC = colix;
} else {
colixA = vertexColixes[iA];
colixB = vertexColixes[iB];
colixC = vertexColixes[iC];
if (this.isBicolorMap) {
if (colixA != colixB || colixB != colixC) continue;
if (this.isGhostPass) {
colixA = colixB = colixC = JU.C.copyColixTranslucency (this.mesh.slabColix, colixA);
}}}if (fill) {
if (generateSet) {
this.bsPolygonsToExport.set (i);
continue;
}if (iB == iC) {
this.setColix (colixA);
if (iA == iB) this.g3d.fillSphereI (diam, this.screens[iA]);
 else this.g3d.fillCylinder (3, diam, this.screens[iA], this.screens[iB]);
} else if (iShowTriangles) {
this.g3d.fillTriangle (this.screens[iA], colixA, nA, this.screens[iB], colixB, nB, this.screens[iC], colixC, nC, 0.1);
} else {
if (this.isTranslucentInherit && vertexColixes != null) {
colixA = JU.C.copyColixTranslucency (this.mesh.slabColix, vertexColixes[iA]);
colixB = JU.C.copyColixTranslucency (this.mesh.slabColix, vertexColixes[iB]);
colixC = JU.C.copyColixTranslucency (this.mesh.slabColix, vertexColixes[iC]);
}this.g3d.fillTriangle3CN (this.screens[iA], colixA, nA, this.screens[iB], colixB, nB, this.screens[iC], colixC, nC);
}if (this.iShowNormals) this.renderNormals ();
} else {
check &= polygon[3];
if (iShowTriangles) check = 7;
if (check == 0) continue;
this.pt1i.setT (this.screens[iA]);
this.pt2i.setT (this.screens[iB]);
this.pt3i.setT (this.screens[iC]);
this.pt1i.z -= 2;
this.pt2i.z -= 2;
this.pt3i.z -= 2;
if (noColor) {
} else if (colorArrayed) {
this.g3d.setC (this.mesh.fillTriangles ? 4 : contourColixes[polygon[4] % contourColixes.length]);
} else {
this.drawTriangle (this.pt1i, colixA, this.pt2i, colixB, this.pt3i, colixC, check, diam);
continue;
}this.drawTriangle (this.pt1i, colix, this.pt2i, colix, this.pt3i, colix, check, diam);
}}
if (generateSet) this.exportSurface (colorSolid ? colix : 0);
}, "~B,~B,~B");
Clazz_defineMethod (c$, "renderNormals", 
 function () {
if (!this.g3d.setC (8)) return;
this.g3d.setFontFid (this.g3d.getFontFidFS ("Monospaced", 24));
var vertexVectors = JU.Normix.getVertexVectors ();
for (var i = this.vertexCount; --i >= 0; ) {
if (this.vertexValues != null && Float.isNaN (this.vertexValues[i])) continue;
if (i > 100) continue;
this.ptTemp.setT (this.vertices[i]);
var n = this.mesh.normixes[i];
if (n >= 0) {
this.ptTemp.scaleAdd2 (3, vertexVectors[n], this.ptTemp);
this.tm.transformPtScr (this.ptTemp, this.ptTempi);
this.g3d.drawLineAB (this.screens[i], this.ptTempi);
}}
});
});
})(Clazz
,Clazz.doubleToInt
,Clazz.declarePackage
,Clazz.instanceOf
,Clazz.load
,Clazz.instantialize
,Clazz.decorateAsClass
,Clazz.floatToInt
,Clazz.makeConstructor
,Clazz.defineEnumConstant
,Clazz.exceptionOf
,Clazz.newIntArray
,Clazz.defineStatics
,Clazz.newFloatArray
,Clazz.declareType
,Clazz.prepareFields
,Clazz.superConstructor
,Clazz.newByteArray
,Clazz.declareInterface
,Clazz.p0p
,Clazz.pu$h
,Clazz.newShortArray
,Clazz.innerTypeInstance
,Clazz.isClassDefined
,Clazz.prepareCallback
,Clazz.newArray
,Clazz.castNullAs
,Clazz.floatToShort
,Clazz.superCall
,Clazz.decorateAsType
,Clazz.newBooleanArray
,Clazz.newCharArray
,Clazz.implementOf
,Clazz.newDoubleArray
,Clazz.overrideConstructor
,Clazz.clone
,Clazz.doubleToShort
,Clazz.getInheritedLevel
,Clazz.getParamsType
,Clazz.isAF
,Clazz.isAI
,Clazz.isAS
,Clazz.isASS
,Clazz.isAP
,Clazz.isAFloat
,Clazz.isAII
,Clazz.isAFF
,Clazz.isAFFF
,Clazz.tryToSearchAndExecute
,Clazz.getStackTrace
,Clazz.inheritArgs
,Clazz.alert
,Clazz.defineMethod
,Clazz.overrideMethod
,Clazz.declareAnonymous
//,Clazz.checkPrivateMethod
,Clazz.cloneFinals
);
