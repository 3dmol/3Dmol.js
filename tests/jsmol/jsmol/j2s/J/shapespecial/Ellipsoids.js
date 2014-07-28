Clazz.declarePackage ("J.shapespecial");
Clazz.load (["J.shape.Shape", "java.util.Hashtable"], "J.shapespecial.Ellipsoids", ["JU.BS", "$.Lst", "$.PT", "$.SB", "$.V3", "J.c.PAL", "J.shapespecial.Ellipsoid", "JU.BSUtil", "$.C", "$.Escape", "$.Txt"], function () {
c$ = Clazz.decorateAsClass (function () {
this.simpleEllipsoids = null;
this.atomEllipsoids = null;
this.typeSelected = "1";
this.selectedAtoms = null;
this.ellipsoidSet = null;
Clazz.instantialize (this, arguments);
}, J.shapespecial, "Ellipsoids", J.shape.Shape);
Clazz.prepareFields (c$, function () {
this.simpleEllipsoids =  new java.util.Hashtable ();
this.atomEllipsoids =  new java.util.Hashtable ();
});
Clazz.defineMethod (c$, "isActive", 
function () {
return !this.atomEllipsoids.isEmpty () || !this.simpleEllipsoids.isEmpty ();
});
Clazz.overrideMethod (c$, "getIndexFromName", 
function (thisID) {
return (this.checkID (thisID) ? 1 : -1);
}, "~S");
Clazz.overrideMethod (c$, "setSize", 
function (size, bsSelected) {
if (this.ms.at == null || size == 0 && this.ms.atomTensors == null) return;
var isAll = (bsSelected == null);
if (!isAll && this.selectedAtoms != null) bsSelected = this.selectedAtoms;
var tensors = this.vwr.ms.getAllAtomTensors (this.typeSelected);
if (tensors == null) return;
var atoms = this.ms.at;
for (var i = tensors.size (); --i >= 0; ) {
var t = tensors.get (i);
if (isAll || t.isSelected (bsSelected, -1)) {
var e = this.atomEllipsoids.get (t);
var isNew = (size != 0 && e == null);
if (isNew) this.atomEllipsoids.put (t, e = J.shapespecial.Ellipsoid.getEllipsoidForAtomTensor (t, atoms[t.atomIndex1]));
if (e != null && (isNew || size != 2147483647)) {
e.setScale (size, true);
}}}
}, "~N,JU.BS");
Clazz.overrideMethod (c$, "getPropertyData", 
function (property, data) {
if (property === "checkID") {
return (this.checkID (data[0]));
}return false;
}, "~S,~A");
Clazz.defineMethod (c$, "checkID", 
 function (thisID) {
this.ellipsoidSet =  new JU.Lst ();
if (thisID == null) return false;
thisID = thisID.toLowerCase ();
if (JU.Txt.isWild (thisID)) {
for (var e, $e = this.simpleEllipsoids.entrySet ().iterator (); $e.hasNext () && ((e = $e.next ()) || true);) {
var key = e.getKey ().toLowerCase ();
if (JU.Txt.isMatch (key, thisID, true, true)) this.ellipsoidSet.addLast (e.getValue ());
}
}var e = this.simpleEllipsoids.get (thisID);
if (e != null) this.ellipsoidSet.addLast (e);
return (this.ellipsoidSet.size () > 0);
}, "~S");
Clazz.defineMethod (c$, "initEllipsoids", 
 function (value) {
var haveID = (value != null);
this.checkID (value);
if (haveID) this.typeSelected = null;
this.selectedAtoms = null;
return haveID;
}, "~O");
Clazz.overrideMethod (c$, "initShape", 
function () {
this.setProperty ("thisID", null, null);
});
Clazz.overrideMethod (c$, "setProperty", 
function (propertyName, value, bs) {
if (propertyName === "thisID") {
if (this.initEllipsoids (value) && this.ellipsoidSet.size () == 0) {
var id = value;
var e = J.shapespecial.Ellipsoid.getEmptyEllipsoid (id, this.vwr.am.cmi);
this.ellipsoidSet.addLast (e);
this.simpleEllipsoids.put (id, e);
}return;
}if ("atoms" === propertyName) {
this.selectedAtoms = value;
return;
}if (propertyName === "deleteModelAtoms") {
var modelIndex = ((value)[2])[0];
var e = this.simpleEllipsoids.values ().iterator ();
while (e.hasNext ()) if (e.next ().tensor.modelIndex == modelIndex) e.remove ();

e = this.atomEllipsoids.values ().iterator ();
while (e.hasNext ()) if (e.next ().modelIndex == modelIndex) e.remove ();

this.ellipsoidSet.clear ();
return;
}var mode = "ax ce co de eq mo on op sc tr".indexOf ((propertyName + "  ").substring (0, 2));
if (this.ellipsoidSet.size () > 0) {
if ("translucentLevel" === propertyName) {
this.setPropS (propertyName, value, bs);
return;
}if (mode >= 0) for (var i = this.ellipsoidSet.size (); --i >= 0; ) this.setProp (this.ellipsoidSet.get (i), Clazz.doubleToInt (mode / 3), value);

return;
}if ("color" === propertyName) {
var colix = JU.C.getColixO (value);
var pid = J.c.PAL.pidOf (value);
if (this.selectedAtoms != null) bs = this.selectedAtoms;
for (var e, $e = this.atomEllipsoids.values ().iterator (); $e.hasNext () && ((e = $e.next ()) || true);) if (e.tensor.type.equals (this.typeSelected) && e.tensor.isSelected (bs, -1)) {
e.colix = this.getColixI (colix, pid, e.tensor.atomIndex1);
e.pid = pid;
}
return;
}if ("on" === propertyName) {
var isOn = (value).booleanValue ();
if (this.selectedAtoms != null) bs = this.selectedAtoms;
if (isOn) this.setSize (2147483647, bs);
for (var e, $e = this.atomEllipsoids.values ().iterator (); $e.hasNext () && ((e = $e.next ()) || true);) {
var t = e.tensor;
if ((t.type.equals (this.typeSelected) || this.typeSelected.equals (t.altType)) && t.isSelected (bs, -1)) {
e.isOn = isOn;
}}
return;
}if ("options" === propertyName) {
var options = (value).toLowerCase ().trim ();
if (options.length == 0) options = null;
if (this.selectedAtoms != null) bs = this.selectedAtoms;
if (options != null) this.setSize (2147483647, bs);
for (var e, $e = this.atomEllipsoids.values ().iterator (); $e.hasNext () && ((e = $e.next ()) || true);) if (e.tensor.type.equals (this.typeSelected) && e.tensor.isSelected (bs, -1)) e.options = options;

return;
}if ("params" === propertyName) {
var data = value;
data[2] = null;
this.typeSelected = "0";
this.setSize (50, bs);
}if ("points" === propertyName) {
return;
}if ("scale" === propertyName) {
this.setSize (Clazz.floatToInt ((value).floatValue () * 100), bs);
return;
}if ("select" === propertyName) {
this.typeSelected = (value).toLowerCase ();
return;
}if ("translucency" === propertyName) {
var isTranslucent = (value.equals ("translucent"));
for (var e, $e = this.atomEllipsoids.values ().iterator (); $e.hasNext () && ((e = $e.next ()) || true);) if (e.tensor.type.equals (this.typeSelected) && e.tensor.isSelected (bs, -1)) e.colix = JU.C.getColixTranslucent3 (e.colix, isTranslucent, this.translucentLevel);

return;
}this.setPropS (propertyName, value, bs);
}, "~S,~O,JU.BS");
Clazz.defineMethod (c$, "setProp", 
 function (e, mode, value) {
switch (mode) {
case 0:
e.setAxes (value);
return;
case 1:
e.setCenter (value);
return;
case 2:
e.colix = JU.C.getColixO (value);
return;
case 3:
this.simpleEllipsoids.remove (e.id);
return;
case 4:
e.setEquation (value);
return;
case 5:
e.tensor.modelIndex = (value).intValue ();
return;
case 6:
e.isOn = (value).booleanValue ();
return;
case 7:
e.options = (value).toLowerCase ();
return;
case 8:
e.setScale ((value).floatValue (), false);
return;
case 9:
e.colix = JU.C.getColixTranslucent3 (e.colix, value.equals ("translucent"), this.translucentLevel);
return;
}
return;
}, "J.shapespecial.Ellipsoid,~N,~O");
Clazz.overrideMethod (c$, "getShapeState", 
function () {
if (!this.isActive ()) return "";
var sb =  new JU.SB ();
sb.append ("\n");
if (!this.simpleEllipsoids.isEmpty ()) this.getStateID (sb);
if (!this.atomEllipsoids.isEmpty ()) this.getStateAtoms (sb);
return sb.toString ();
});
Clazz.defineMethod (c$, "getStateID", 
 function (sb) {
var v1 =  new JU.V3 ();
for (var ellipsoid, $ellipsoid = this.simpleEllipsoids.values ().iterator (); $ellipsoid.hasNext () && ((ellipsoid = $ellipsoid.next ()) || true);) {
var t = ellipsoid.tensor;
if (!ellipsoid.isValid || t == null) continue;
sb.append ("  Ellipsoid ID ").append (ellipsoid.id).append (" modelIndex ").appendI (t.modelIndex).append (" center ").append (JU.Escape.eP (ellipsoid.center)).append (" axes");
for (var i = 0; i < 3; i++) {
v1.setT (t.eigenVectors[i]);
v1.scale (ellipsoid.lengths[i]);
sb.append (" ").append (JU.Escape.eP (v1));
}
sb.append (" " + J.shape.Shape.getColorCommandUnk ("", ellipsoid.colix, this.translucentAllowed));
if (ellipsoid.options != null) sb.append (" options ").append (JU.PT.esc (ellipsoid.options));
if (!ellipsoid.isOn) sb.append (" off");
sb.append (";\n");
}
}, "JU.SB");
Clazz.defineMethod (c$, "getStateAtoms", 
 function (sb) {
var bsDone =  new JU.BS ();
var temp =  new java.util.Hashtable ();
var temp2 =  new java.util.Hashtable ();
for (var e, $e = this.atomEllipsoids.values ().iterator (); $e.hasNext () && ((e = $e.next ()) || true);) {
var iType = e.tensor.iType;
if (bsDone.get (iType + 1)) continue;
bsDone.set (iType + 1);
var isADP = (e.tensor.iType == 1);
var cmd = (isADP ? null : "Ellipsoids set " + JU.PT.esc (e.tensor.type));
for (var e2, $e2 = this.atomEllipsoids.values ().iterator (); $e2.hasNext () && ((e2 = $e2.next ()) || true);) {
if (e2.tensor.iType != iType || isADP && !e2.isOn) continue;
var i = e2.tensor.atomIndex1;
JU.BSUtil.setMapBitSet (temp, i, i, (isADP ? "Ellipsoids " + e2.percent : cmd + " scale " + e2.scale + (e2.options == null ? "" : " options " + JU.PT.esc (e2.options)) + (e2.isOn ? " ON" : " OFF")));
if (e2.colix != 0) JU.BSUtil.setMapBitSet (temp2, i, i, J.shape.Shape.getColorCommand (cmd, e2.pid, e2.colix, this.translucentAllowed));
}
}
sb.append (this.vwr.getCommands (temp, temp2, "select"));
}, "JU.SB");
Clazz.overrideMethod (c$, "setVisibilityFlags", 
function (bsModels) {
if (!this.isActive ()) return;
var atoms = this.vwr.ms.at;
this.setVis (this.simpleEllipsoids, bsModels, atoms);
if (this.atomEllipsoids != null) for (var i = atoms.length; --i >= 0; ) atoms[i].setShapeVisibility (this.vf, false);

this.setVis (this.atomEllipsoids, bsModels, atoms);
}, "JU.BS");
Clazz.defineMethod (c$, "setVis", 
 function (ellipsoids, bs, atoms) {
for (var e, $e = ellipsoids.values ().iterator (); $e.hasNext () && ((e = $e.next ()) || true);) {
var t = e.tensor;
var isOK = (t != null && e.isValid && e.isOn);
if (isOK && t.atomIndex1 >= 0) {
if (t.iType == 1) {
var isModTensor = t.isModulated;
var isUnmodTensor = t.isUnmodulated;
var isModAtom = this.ms.isModulated (t.atomIndex1);
isOK = (!isModTensor && !isUnmodTensor || isModTensor == isModAtom);
}atoms[t.atomIndex1].setShapeVisibility (this.vf, true);
}e.visible = isOK && (e.modelIndex < 0 || bs.get (e.modelIndex));
}
}, "java.util.Map,JU.BS,~A");
Clazz.overrideMethod (c$, "setModelClickability", 
function () {
if (this.atomEllipsoids.isEmpty ()) return;
for (var e, $e = this.atomEllipsoids.values ().iterator (); $e.hasNext () && ((e = $e.next ()) || true);) {
var i = e.tensor.atomIndex1;
var atom = this.ms.at[i];
if ((atom.shapeVisibilityFlags & this.vf) == 0 || this.ms.isAtomHidden (i)) continue;
atom.setClickable (this.vf);
}
});
Clazz.defineStatics (c$,
"PROPERTY_MODES", "ax ce co de eq mo on op sc tr");
});
