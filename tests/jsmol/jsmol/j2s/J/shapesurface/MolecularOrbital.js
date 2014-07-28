Clazz.declarePackage ("J.shapesurface");
Clazz.load (["J.shapesurface.Isosurface"], "J.shapesurface.MolecularOrbital", ["java.lang.Boolean", "$.Float", "java.util.Hashtable", "JU.AU", "$.Lst", "$.PT", "$.SB", "J.c.QS", "J.jvxl.data.JvxlCoder", "JU.Escape"], function () {
c$ = Clazz.decorateAsClass (function () {
this.moTranslucency = null;
this.moTranslucentLevel = null;
this.moPlane = null;
this.moCutoff = null;
this.moResolution = null;
this.moScale = null;
this.moColorPos = null;
this.moColorNeg = null;
this.moMonteCarloCount = null;
this.moIsPositiveOnly = false;
this.moSquareData = null;
this.moSquareLinear = null;
this.moRandomSeed = null;
this.moFill = 1073742046;
this.moMesh = 1073742018;
this.moDots = 1073742042;
this.moFrontOnly = 1073741960;
this.moTitleFormat = null;
this.moDebug = false;
this.myColorPt = 0;
this.strID = null;
this.$moNumber = 0;
this.$moLinearCombination = null;
this.htModels = null;
this.thisModel = null;
this.moSlab = null;
this.moSlabValue = null;
Clazz.instantialize (this, arguments);
}, J.shapesurface, "MolecularOrbital", J.shapesurface.Isosurface);
Clazz.defineMethod (c$, "initShape", 
function () {
Clazz.superCall (this, J.shapesurface.MolecularOrbital, "initShape", []);
this.myType = "mo";
this.setPropI ("thisID", "mo", null);
});
Clazz.overrideMethod (c$, "setProperty", 
function (propertyName, value, bs) {
if ("init" === propertyName) {
this.myColorPt = 0;
this.moDebug = false;
var modelIndex = (value).intValue ();
this.strID = this.getId (modelIndex);
this.setPropI ("init", null, null);
this.setPropI ("modelIndex", Integer.$valueOf (modelIndex), null);
if (this.htModels == null) this.htModels =  new java.util.Hashtable ();
if (!this.htModels.containsKey (this.strID)) this.htModels.put (this.strID,  new java.util.Hashtable ());
this.thisModel = this.htModels.get (this.strID);
this.$moNumber = (!this.thisModel.containsKey ("moNumber") ? 0 : (this.thisModel.get ("moNumber")).intValue ());
this.$moLinearCombination = this.thisModel.get ("moLinearCombination");
this.moSquareData = this.moSquareLinear = null;
return;
}if ("slab" === propertyName) {
if (Clazz.instanceOf (value, Integer)) {
this.thisModel.put ("slabValue", value);
} else {
var slabInfo = value;
var tok = (slabInfo[0]).intValue ();
this.moSlab = this.thisModel.get ("slab");
if (this.moSlab == null) this.thisModel.put ("slab", this.moSlab =  new JU.Lst ());
if (tok == 1048587) {
this.moSlab = null;
this.thisModel.remove ("slab");
return;
}this.moSlab.addLast (value);
}return;
}if ("cutoff" === propertyName) {
this.thisModel.put ("moCutoff", value);
this.thisModel.put ("moIsPositiveOnly", Boolean.FALSE);
return;
}if ("scale" === propertyName) {
this.thisModel.put ("moScale", value);
return;
}if ("squareData" === propertyName) {
this.thisModel.put ("moSquareData", Boolean.TRUE);
this.moSquareData = Boolean.TRUE;
return;
}if ("squareLinear" === propertyName) {
this.thisModel.put ("moSquareLinear", Boolean.TRUE);
this.moSquareLinear = Boolean.TRUE;
return;
}if ("cutoffPositive" === propertyName) {
this.thisModel.put ("moCutoff", value);
this.thisModel.put ("moIsPositiveOnly", Boolean.TRUE);
return;
}if ("resolution" === propertyName) {
this.thisModel.put ("moResolution", value);
return;
}if ("titleFormat" === propertyName) {
this.moTitleFormat = value;
return;
}if ("color" === propertyName) {
if (!(Clazz.instanceOf (value, Integer))) return;
this.thisModel.remove ("moTranslucency");
this.setPropI ("color", value, bs);
propertyName = "colorRGB";
this.myColorPt = 0;
}if ("colorRGB" === propertyName) {
this.moColorPos = value;
if (this.myColorPt++ == 0) this.moColorNeg = this.moColorPos;
this.thisModel.put ("moColorNeg", this.moColorNeg);
this.thisModel.put ("moColorPos", this.moColorPos);
return;
}if ("plane" === propertyName) {
if (value == null) this.thisModel.remove ("moPlane");
 else this.thisModel.put ("moPlane", value);
return;
}if ("monteCarloCount" === propertyName) {
this.thisModel.put ("monteCarloCount", value);
return;
}if ("randomSeed" === propertyName) {
if (value == null) this.thisModel.remove ("randomSeed");
 else this.thisModel.put ("randomSeed", value);
return;
}if ("molecularOrbital" === propertyName) {
if (Clazz.instanceOf (value, Integer)) {
this.$moNumber = (value).intValue ();
this.thisModel.put ("moNumber", value);
this.thisModel.remove ("moLinearCombination");
this.$moLinearCombination = null;
} else {
this.$moNumber = 0;
this.$moLinearCombination = value;
this.thisModel.put ("moNumber", Integer.$valueOf (0));
this.thisModel.put ("moLinearCombination", this.$moLinearCombination);
}if (this.moSquareData === Boolean.TRUE) this.thisModel.put ("moSquareData", Boolean.TRUE);
 else this.thisModel.remove ("moSquareData");
if (this.moSquareLinear === Boolean.TRUE) this.thisModel.put ("moSquareLinear", Boolean.TRUE);
 else this.thisModel.remove ("moSquareLinear");
this.setOrbital (this.$moNumber, this.$moLinearCombination);
return;
}if ("translucentLevel" === propertyName) {
if (this.thisModel == null) {
if (this.currentMesh == null) return;
this.thisModel = this.htModels.get (this.currentMesh.thisID);
}this.thisModel.put ("moTranslucentLevel", value);
}if ("delete" === propertyName) {
this.htModels.remove (this.strID);
this.$moNumber = 0;
this.$moLinearCombination = null;
}if ("token" === propertyName) {
var tok = (value).intValue ();
switch (tok) {
case 1113198595:
case 1073742042:
this.moDots = tok;
break;
case 1073741938:
case 1073742046:
this.moFill = tok;
break;
case 1073742018:
case 1073742052:
this.moMesh = tok;
break;
case 1073741960:
case 1073742058:
this.moFrontOnly = tok;
break;
}
}if ("translucency" === propertyName) {
if (this.thisModel == null) {
if (this.currentMesh == null) return;
this.thisModel = this.htModels.get (this.currentMesh.thisID);
}this.thisModel.put ("moTranslucency", value);
}if (propertyName === "deleteModelAtoms") {
var modelIndex = ((value)[2])[0];
var htModelsNew =  new java.util.Hashtable ();
for (var i = this.meshCount; --i >= 0; ) {
if (this.meshes[i] == null) continue;
if (this.meshes[i].modelIndex == modelIndex) {
this.meshCount--;
if (this.meshes[i] === this.currentMesh) {
this.currentMesh = null;
this.thisModel = null;
}this.meshes = JU.AU.deleteElements (this.meshes, i, 1);
continue;
}var htModel = this.htModels.get (this.meshes[i].thisID);
if (this.meshes[i].modelIndex > modelIndex) {
this.meshes[i].modelIndex--;
this.meshes[i].thisID = this.getId (this.meshes[i].modelIndex);
}htModelsNew.put (this.meshes[i].thisID, htModel);
}
this.htModels = htModelsNew;
return;
}this.setPropI (propertyName, value, bs);
}, "~S,~O,JU.BS");
Clazz.defineMethod (c$, "getId", 
 function (modelIndex) {
return "mo_model" + this.vwr.getModelNumberDotted (modelIndex);
}, "~N");
Clazz.overrideMethod (c$, "getProperty", 
function (propertyName, param) {
if (propertyName.equals ("list")) {
var s = this.getPropI ("list");
if (s.length > 1) s += "cutoff = " + this.getPropI ("cutoff") + "\n";
return this.vwr.getMoInfo (-1) + "\n" + s;
}if (propertyName === "moNumber") return Integer.$valueOf (this.$moNumber);
if (propertyName === "moLinearCombination") return this.$moLinearCombination;
if (propertyName === "showMO") {
var str =  new JU.SB ();
var mos = (this.sg.getMoData ().get ("mos"));
var nOrb = (mos == null ? 0 : mos.size ());
var thisMO = param;
var currentMO = this.$moNumber;
var isShowCurrent = (thisMO == -2147483648);
if (thisMO == 2147483647) {
thisMO = currentMO;
}if (nOrb == 0 || isShowCurrent && currentMO == 0) return "";
var doOneMo = (thisMO != 0);
if (currentMO == 0) thisMO = 0;
var haveHeader = false;
var nTotal = (thisMO > 0 ? 1 : nOrb);
var i0 = (nTotal == 1 && currentMO > 0 ? currentMO : 1);
for (var i = i0; i <= nOrb; i++) if (thisMO == 0 || thisMO == i || !doOneMo && i == currentMO) {
if (!doOneMo) {
var params = this.sg.getParams ();
this.setPropI ("init", params, null);
this.setOrbital (i, null);
}this.jvxlData.moleculeXml = this.vwr.getModelCml (this.vwr.getModelUndeletedAtomsBitSet (this.thisMesh.modelIndex), 100, true);
if (!haveHeader) {
str.append (J.jvxl.data.JvxlCoder.jvxlGetFile (this.jvxlData, null, null, "HEADERONLY", true, nTotal, null, null));
haveHeader = true;
}str.append (J.jvxl.data.JvxlCoder.jvxlGetFile (this.jvxlData, null, this.jvxlData.title, null, false, 1, this.thisMesh.getState ("mo"), (this.thisMesh.scriptCommand == null ? "" : this.thisMesh.scriptCommand)));
if (!doOneMo) this.setPropI ("delete", "mo_show", null);
if (nTotal == 1) break;
}
str.append (J.jvxl.data.JvxlCoder.jvxlGetFile (this.jvxlData, null, null, "TRAILERONLY", true, 0, null, null));
return str.toString ();
}return this.getPropI (propertyName);
}, "~S,~N");
Clazz.overrideMethod (c$, "clearSg", 
function () {
});
Clazz.defineMethod (c$, "getSettings", 
 function (strID) {
this.thisModel = this.htModels.get (strID);
if (this.thisModel == null || this.thisModel.get ("moNumber") == null) return false;
this.moTranslucency = this.thisModel.get ("moTranslucency");
this.moTranslucentLevel = this.thisModel.get ("moTranslucentLevel");
this.moPlane = this.thisModel.get ("moPlane");
this.moCutoff = this.thisModel.get ("moCutoff");
if (this.moCutoff == null) this.moCutoff = this.sg.getMoData ().get ("defaultCutoff");
if (this.moCutoff == null) {
this.moCutoff = Float.$valueOf (0.05);
}this.thisModel.put ("moCutoff", Float.$valueOf (this.moCutoff.floatValue ()));
this.moResolution = this.thisModel.get ("moResolution");
this.moScale = this.thisModel.get ("moScale");
this.moColorPos = this.thisModel.get ("moColorPos");
this.moColorNeg = this.thisModel.get ("moColorNeg");
this.moSquareData = this.thisModel.get ("moSquareData");
this.moSquareLinear = this.thisModel.get ("moSquareLinear");
this.moMonteCarloCount = this.thisModel.get ("monteCarloCount");
this.moRandomSeed = this.thisModel.get ("randomSeed");
this.moSlabValue = this.thisModel.get ("slabValue");
this.moSlab = this.thisModel.get ("slab");
if (this.moRandomSeed == null) this.thisModel.put ("randomSeed", this.moRandomSeed = Integer.$valueOf ((-System.currentTimeMillis ()) % 10000));
this.$moNumber = (this.thisModel.get ("moNumber")).intValue ();
this.$moLinearCombination = this.thisModel.get ("moLinearCombination");
var b = this.thisModel.get ("moIsPositiveOnly");
this.moIsPositiveOnly = (b != null && ((b)).booleanValue ());
return true;
}, "~S");
Clazz.defineMethod (c$, "setOrbital", 
 function (moNumber, linearCombination) {
this.setPropI ("reset", this.strID, null);
if (this.moDebug) this.setPropI ("debug", Boolean.TRUE, null);
this.getSettings (this.strID);
if (this.moScale != null) this.setPropI ("scale", this.moScale, null);
if (this.moResolution != null) this.setPropI ("resolution", this.moResolution, null);
if (this.moPlane != null) {
this.setPropI ("plane", this.moPlane, null);
if (this.moCutoff != null) {
this.setPropI ("red", Float.$valueOf (-this.moCutoff.floatValue ()), null);
this.setPropI ("blue", this.moCutoff, null);
}} else {
if (this.moCutoff != null) this.setPropI ((this.moIsPositiveOnly ? "cutoffPositive" : "cutoff"), this.moCutoff, null);
if (this.moColorNeg != null) this.setPropI ("colorRGB", this.moColorNeg, null);
if (this.moColorPos != null) this.setPropI ("colorRGB", this.moColorPos, null);
if (this.moMonteCarloCount != null) {
this.setPropI ("randomSeed", this.moRandomSeed, null);
this.setPropI ("monteCarloCount", this.moMonteCarloCount, null);
}}this.setPropI ("squareData", this.moSquareData, null);
this.setPropI ("squareLinear", this.moSquareLinear, null);
this.setPropI ("title", this.moTitleFormat, null);
this.setPropI ("fileName", this.vwr.getFileName (), null);
this.setPropI ("molecularOrbital", linearCombination == null ? Integer.$valueOf (moNumber) : linearCombination, null);
if (this.moPlane != null && this.moColorNeg != null) this.setPropI ("colorRGB", this.moColorNeg, null);
if (this.moPlane != null && this.moColorPos != null) this.setPropI ("colorRGB", this.moColorPos, null);
this.currentMesh.isColorSolid = false;
if (this.moSlabValue != null) this.setPropI ("slab", this.moSlabValue, null);
if (this.moSlab != null) for (var i = 0; i < this.moSlab.size (); i++) this.setPropI ("slab", this.moSlab.get (i), null);

if (this.moTranslucentLevel != null) this.setPropI ("translucentLevel", this.moTranslucentLevel, null);
if (this.moTranslucency != null) this.setPropI ("translucency", this.moTranslucency, null);
this.setPropI ("token", Integer.$valueOf (this.moFill), null);
this.setPropI ("token", Integer.$valueOf (this.moMesh), null);
this.setPropI ("token", Integer.$valueOf (this.moDots), null);
this.setPropI ("token", Integer.$valueOf (this.moFrontOnly), null);
this.thisModel.put ("mesh", this.currentMesh);
return;
}, "~N,~A");
Clazz.overrideMethod (c$, "getShapeState", 
function () {
if (this.htModels == null) return "";
var s =  new JU.SB ();
var modelCount = this.vwr.getModelCount ();
for (var i = 0; i < modelCount; i++) s.append (this.getMoState (i));

return s.toString ();
});
Clazz.defineMethod (c$, "getMoState", 
 function (modelIndex) {
this.strID = this.getId (modelIndex);
if (!this.getSettings (this.strID)) return "";
var s =  new JU.SB ();
var modelCount = this.vwr.getModelCount ();
if (modelCount > 1) J.shape.Shape.appendCmd (s, "frame " + this.vwr.getModelNumberDotted (modelIndex));
if (this.moCutoff != null) J.shape.Shape.appendCmd (s, "mo cutoff " + (this.sg.getIsPositiveOnly () ? "+" : "") + this.moCutoff);
if (this.moScale != null) J.shape.Shape.appendCmd (s, "mo scale " + this.moScale);
if (this.moMonteCarloCount != null) J.shape.Shape.appendCmd (s, "mo points " + this.moMonteCarloCount + " " + this.moRandomSeed);
if (this.moResolution != null) J.shape.Shape.appendCmd (s, "mo resolution " + this.moResolution);
if (this.moPlane != null) J.shape.Shape.appendCmd (s, "mo plane {" + this.moPlane.x + " " + this.moPlane.y + " " + this.moPlane.z + " " + this.moPlane.w + "}");
if (this.moTitleFormat != null) J.shape.Shape.appendCmd (s, "mo titleFormat " + JU.PT.esc (this.moTitleFormat));
if (this.moColorNeg != null) J.shape.Shape.appendCmd (s, "mo color " + JU.Escape.escapeColor (this.moColorNeg.intValue ()) + (this.moColorNeg.equals (this.moColorPos) ? "" : " " + JU.Escape.escapeColor (this.moColorPos.intValue ())));
if (this.moSlab != null) {
if (this.thisMesh.slabOptions != null) J.shape.Shape.appendCmd (s, this.thisMesh.slabOptions.toString ());
if (this.thisMesh.jvxlData.slabValue != -2147483648) J.shape.Shape.appendCmd (s, "mo slab " + this.thisMesh.jvxlData.slabValue);
}if (this.$moLinearCombination == null) {
J.shape.Shape.appendCmd (s, "mo " + (this.moSquareData === Boolean.TRUE ? "squared " : "") + this.$moNumber);
} else {
J.shape.Shape.appendCmd (s, "mo " + J.c.QS.getMOString (this.$moLinearCombination) + (this.moSquareLinear === Boolean.TRUE ? " squared" : ""));
}if (this.moTranslucency != null) J.shape.Shape.appendCmd (s, "mo translucent " + this.moTranslucentLevel);
J.shape.Shape.appendCmd (s, (this.thisModel.get ("mesh")).getState ("mo"));
return s.toString ();
}, "~N");
Clazz.defineMethod (c$, "merge", 
function (shape) {
var mo = shape;
this.moColorNeg = mo.moColorNeg;
this.moColorPos = mo.moColorPos;
this.moCutoff = mo.moCutoff;
this.moPlane = mo.moPlane;
this.moResolution = mo.moResolution;
this.moScale = mo.moScale;
this.moSlab = mo.moSlab;
this.moSlabValue = mo.moSlabValue;
this.moTitleFormat = mo.moTitleFormat;
this.moTranslucency = mo.moTranslucency;
if (this.htModels == null) this.htModels =  new java.util.Hashtable ();
var ht = mo.htModels;
if (ht != null) {
for (var entry, $entry = ht.entrySet ().iterator (); $entry.hasNext () && ((entry = $entry.next ()) || true);) {
var key = entry.getKey ();
this.htModels.put (key, entry.getValue ());
}
}Clazz.superCall (this, J.shapesurface.MolecularOrbital, "merge", [shape]);
}, "J.shape.Shape");
});
