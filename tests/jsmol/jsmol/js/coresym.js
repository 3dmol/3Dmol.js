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
Clazz_declarePackage ("J.adapter.smarter");
Clazz_load (["JU.P3"], "J.adapter.smarter.XtalSymmetry", ["java.lang.Boolean", "$.Float", "java.util.Hashtable", "JU.BS", "$.Lst", "$.M3", "$.M4", "$.P3i", "$.PT", "$.SB", "$.V3", "J.adapter.smarter.Atom", "J.api.Interface", "JU.BSUtil", "$.Logger", "$.Vibration"], function () {
c$ = Clazz_decorateAsClass (function () {
this.asc = null;
this.symmetry = null;
this.notionalUnitCell = null;
this.symmetryRange = 0;
this.doCentroidUnitCell = false;
this.centroidPacked = false;
this.packingError = 0;
this.filterSymop = null;
this.applySymmetryToBonds = false;
this.latticeCells = null;
this.ptSupercell = null;
this.fmatSupercell = null;
this.matSupercell = null;
this.trajectoryUnitCells = null;
this.doNormalize = true;
this.doPackUnitCell = false;
this.vibsFractional = false;
this.rminx = 0;
this.rminy = 0;
this.rminz = 0;
this.rmaxx = 0;
this.rmaxy = 0;
this.rmaxz = 0;
this.ptOffset = null;
this.unitCellOffset = null;
this.minXYZ = null;
this.maxXYZ = null;
this.minXYZ0 = null;
this.maxXYZ0 = null;
this.checkAll = false;
this.bondCount0 = 0;
this.dtype = 3;
this.unitCellTranslations = null;
this.latticeOp = 0;
this.latticeOnly = false;
this.noSymmetryCount = 0;
this.firstSymmetryAtom = 0;
this.ptTemp = null;
this.mTemp = null;
this.nVib = 0;
Clazz_instantialize (this, arguments);
}, J.adapter.smarter, "XtalSymmetry");
Clazz_prepareFields (c$, function () {
this.notionalUnitCell =  Clazz_newFloatArray (6, 0);
this.ptOffset =  new JU.P3 ();
});
Clazz_makeConstructor (c$, 
function () {
});
Clazz_defineMethod (c$, "set", 
function (asc) {
this.asc = asc;
this.getSymmetry ();
return this;
}, "J.adapter.smarter.AtomSetCollection");
Clazz_defineMethod (c$, "getSymmetry", 
function () {
return (this.symmetry == null ? (this.symmetry = J.api.Interface.getSymmetry ()) : this.symmetry);
});
Clazz_defineMethod (c$, "setSymmetry", 
function (symmetry) {
return (this.symmetry = symmetry);
}, "J.api.SymmetryInterface");
Clazz_defineMethod (c$, "setSymmetryRange", 
 function (factor) {
this.symmetryRange = factor;
this.asc.setInfo ("symmetryRange", Float.$valueOf (factor));
}, "~N");
Clazz_defineMethod (c$, "setLatticeCells", 
 function (acr) {
this.latticeCells = acr.latticeCells;
var isLatticeRange = (this.latticeCells[0] <= 555 && this.latticeCells[1] >= 555 && (this.latticeCells[2] == 0 || this.latticeCells[2] == 1 || this.latticeCells[2] == -1));
this.doNormalize = this.latticeCells[0] != 0 && (!isLatticeRange || this.latticeCells[2] == 1);
this.applySymmetryToBonds = acr.applySymmetryToBonds;
this.doPackUnitCell = acr.doPackUnitCell;
this.doCentroidUnitCell = acr.doCentroidUnitCell;
this.centroidPacked = acr.centroidPacked;
this.filterSymop = acr.filterSymop;
if (acr.strSupercell == null) this.setSupercellFromPoint (acr.ptSupercell);
 else this.setSuperCellFromString (acr.strSupercell);
}, "J.adapter.smarter.AtomSetCollectionReader");
Clazz_defineMethod (c$, "setSupercellFromPoint", 
function (pt) {
this.ptSupercell = pt;
if (pt == null) {
this.matSupercell = null;
return;
}this.matSupercell =  new JU.M4 ();
this.matSupercell.m00 = pt.x;
this.matSupercell.m11 = pt.y;
this.matSupercell.m22 = pt.z;
this.matSupercell.m33 = 1;
JU.Logger.info ("Using supercell \n" + this.matSupercell);
}, "JU.P3");
Clazz_defineMethod (c$, "setSuperCellFromString", 
 function (supercell) {
if (this.fmatSupercell != null) return;
this.fmatSupercell =  Clazz_newFloatArray (16, 0);
if (this.symmetry.getMatrixFromString (supercell, this.fmatSupercell, true, 0) == null) {
this.fmatSupercell = null;
this.matSupercell = null;
return;
}this.matSupercell = JU.M4.newA16 (this.fmatSupercell);
JU.Logger.info ("Using supercell \n" + this.matSupercell);
}, "~S");
Clazz_defineMethod (c$, "setNotionalUnitCell", 
 function (info, matUnitCellOrientation, unitCellOffset) {
this.notionalUnitCell =  Clazz_newFloatArray (info.length, 0);
this.unitCellOffset = unitCellOffset;
for (var i = 0; i < info.length; i++) this.notionalUnitCell[i] = info[i];

this.asc.haveUnitCell = true;
this.asc.setAtomSetAuxiliaryInfo ("notionalUnitcell", this.notionalUnitCell);
if (this.asc.isTrajectory) {
if (this.trajectoryUnitCells == null) {
this.trajectoryUnitCells =  new JU.Lst ();
this.asc.setInfo ("unitCells", this.trajectoryUnitCells);
}this.trajectoryUnitCells.addLast (this.notionalUnitCell);
}this.asc.setGlobalBoolean (2);
this.getSymmetry ().setUnitCell (this.notionalUnitCell, false);
if (unitCellOffset != null) {
this.symmetry.setOffsetPt (unitCellOffset);
this.asc.setAtomSetAuxiliaryInfo ("unitCellOffset", unitCellOffset);
}if (matUnitCellOrientation != null) {
this.symmetry.initializeOrientation (matUnitCellOrientation);
this.asc.setAtomSetAuxiliaryInfo ("matUnitCellOrientation", matUnitCellOrientation);
}}, "~A,JU.M3,JU.P3");
Clazz_defineMethod (c$, "addSpaceGroupOperation", 
function (acr, xyz) {
if (acr != null) this.setLatticeCells (acr);
this.symmetry.setSpaceGroup (this.doNormalize);
return this.symmetry.addSpaceGroupOperation (xyz, 0);
}, "J.adapter.smarter.AtomSetCollectionReader,~S");
Clazz_defineMethod (c$, "setLatticeParameter", 
function (latt) {
this.symmetry.setSpaceGroup (this.doNormalize);
this.symmetry.setLattice (latt);
}, "~N");
Clazz_defineMethod (c$, "applySymmetryFromReader", 
function (acr, readerSymmetry) {
this.asc.setCoordinatesAreFractional (acr.iHaveFractionalCoordinates);
this.setNotionalUnitCell (acr.notionalUnitCell, acr.matUnitCellOrientation, acr.unitCellOffset);
this.asc.setAtomSetSpaceGroupName (acr.sgName);
this.setSymmetryRange (acr.symmetryRange);
if (acr.doConvertToFractional || acr.fileCoordinatesAreFractional) {
this.setLatticeCells (acr);
var doApplySymmetry = true;
if (acr.ignoreFileSpaceGroupName || !acr.iHaveSymmetryOperators) {
if (!acr.merging || readerSymmetry == null) readerSymmetry = acr.getNewSymmetry ();
doApplySymmetry = readerSymmetry.createSpaceGroup (acr.desiredSpaceGroupIndex, (acr.sgName.indexOf ("!") >= 0 ? "P1" : acr.sgName), acr.notionalUnitCell);
} else {
acr.doPreSymmetry ();
this.vibsFractional = acr.vibsFractional;
readerSymmetry = null;
}if (doApplySymmetry) {
if (readerSymmetry != null) this.symmetry.setSpaceGroupFrom (readerSymmetry);
this.packingError = acr.packingError;
this.applySymmetryLattice (acr.ms, acr.altCell);
if (readerSymmetry != null && this.filterSymop == null) this.asc.setAtomSetSpaceGroupName (readerSymmetry.getSpaceGroupName ());
}}if (acr.iHaveFractionalCoordinates && acr.merging && readerSymmetry != null) {
this.asc.toCartesian (readerSymmetry);
this.asc.setCoordinatesAreFractional (false);
acr.addVibrations = false;
}return this.symmetry;
}, "J.adapter.smarter.AtomSetCollectionReader,J.api.SymmetryInterface");
Clazz_defineMethod (c$, "applySymmetryLattice", 
 function (ms, altCell) {
if (!this.asc.coordinatesAreFractional || this.symmetry.getSpaceGroup () == null) return;
var maxX = this.latticeCells[0];
var maxY = this.latticeCells[1];
var maxZ = Math.abs (this.latticeCells[2]);
this.firstSymmetryAtom = this.asc.getLastAtomSetAtomIndex ();
var bsAtoms = null;
this.rminx = this.rminy = this.rminz = 3.4028235E38;
this.rmaxx = this.rmaxy = this.rmaxz = -3.4028235E38;
var oabc = null;
var offset = null;
var pt0 = null;
this.nVib = 0;
var va = null;
var vb = null;
var vc = null;
if (altCell != null && altCell.indexOf (",") >= 0) {
oabc = this.symmetry.getV0abc (altCell);
if (oabc != null) {
this.minXYZ =  new JU.P3i ();
this.maxXYZ = JU.P3i.new3 (maxX, maxY, maxZ);
this.symmetry.setMinMaxLatticeParameters (this.minXYZ, this.maxXYZ);
pt0 = JU.P3.newP (oabc[0]);
this.symmetry.toFractional (pt0, true);
va = JU.P3.newP (oabc[1]);
vb = JU.P3.newP (oabc[2]);
vc = JU.P3.newP (oabc[3]);
oabc[0].scaleAdd2 (this.minXYZ.x, va, oabc[0]);
oabc[0].scaleAdd2 (this.minXYZ.y, vb, oabc[0]);
oabc[0].scaleAdd2 (this.minXYZ.z, vc, oabc[0]);
var pt = JU.P3.newP (oabc[0]);
this.symmetry.toFractional (pt, true);
this.setSymmetryMinMax (pt);
oabc[1].scale (this.maxXYZ.x - this.minXYZ.x);
oabc[2].scale (this.maxXYZ.y - this.minXYZ.y);
oabc[3].scale (this.maxXYZ.z - this.minXYZ.z);
for (var i = 0; i < 3; i++) {
for (var j = i + 1; j < 4; j++) {
pt.add2 (oabc[i], oabc[j]);
if (i != 0) pt.add (oabc[0]);
this.symmetry.toFractional (pt, false);
this.setSymmetryMinMax (pt);
}
}
this.symmetry.toCartesian (pt, false);
pt.add (oabc[1]);
this.symmetry.toFractional (pt, false);
this.setSymmetryMinMax (pt);
this.minXYZ = JU.P3i.new3 (Clazz_doubleToInt (Math.floor (this.rminx + 0.001)), Clazz_doubleToInt (Math.floor (this.rminy + 0.001)), Clazz_doubleToInt (Math.floor (this.rminz + 0.001)));
this.maxXYZ = JU.P3i.new3 (Clazz_doubleToInt (Math.ceil (this.rmaxx - 0.001)), Clazz_doubleToInt (Math.ceil (this.rmaxy - 0.001)), Clazz_doubleToInt (Math.ceil (this.rmaxz - 0.001)));
}} else if (this.fmatSupercell != null) {
var pt =  new JU.P3 ();
for (var i = 0; i <= 1; i++) for (var j = 0; j <= 1; j++) for (var k = 0; k <= 1; k++) {
pt.set (i, j, k);
this.setSym (pt);
}


offset = this.asc.getAtomSetAuxiliaryInfoValue (-1, "unitCellOffset");
this.minXYZ = JU.P3i.new3 (Clazz_floatToInt (this.rminx), Clazz_floatToInt (this.rminy), Clazz_floatToInt (this.rminz));
this.maxXYZ = JU.P3i.new3 (Clazz_floatToInt (this.rmaxx), Clazz_floatToInt (this.rmaxy), Clazz_floatToInt (this.rmaxz));
va = this.setSym (JU.P3.new3 (1, 0, 0));
vb = this.setSym (JU.P3.new3 (0, 1, 0));
vc = this.setSym (JU.P3.new3 (0, 0, 1));
}if (this.rminx == 3.4028235E38) {
this.fmatSupercell = null;
this.matSupercell = null;
altCell = null;
oabc = null;
} else {
JU.Logger.info ("setting min/max for original lattice to " + this.minXYZ + " and " + this.maxXYZ);
var doPack0 = this.doPackUnitCell;
this.doPackUnitCell = (oabc != null);
if (this.asc.bsAtoms == null) this.asc.bsAtoms = JU.BSUtil.setAll (this.asc.ac);
bsAtoms = this.asc.bsAtoms;
this.applyAllSymmetry (ms, null, false);
this.doPackUnitCell = doPack0;
this.setVibVectors ();
var atoms = this.asc.atoms;
var atomCount = this.asc.ac;
var iAtomFirst = this.asc.getLastAtomSetAtomIndex ();
for (var i = iAtomFirst; i < atomCount; i++) {
this.symmetry.toCartesian (atoms[i], true);
}
this.symmetry = null;
this.symmetry = this.getSymmetry ();
this.setNotionalUnitCell ([0, 0, 0, 0, 0, 0, va.x, va.y, va.z, vb.x, vb.y, vb.z, vc.x, vc.y, vc.z], null, offset);
this.asc.setAtomSetSpaceGroupName (oabc == null ? "P1" : "cell=" + altCell);
this.symmetry.setSpaceGroup (this.doNormalize);
this.symmetry.addSpaceGroupOperation ("x,y,z", 0);
for (var i = iAtomFirst; i < atomCount; i++) {
this.symmetry.toFractional (atoms[i], true);
if (pt0 != null) atoms[i].sub (pt0);
}
this.asc.haveAnisou = false;
this.asc.setAtomSetAuxiliaryInfo ("matUnitCellOrientation", null);
}this.minXYZ =  new JU.P3i ();
this.maxXYZ = JU.P3i.new3 (maxX, maxY, maxZ);
if (oabc == null) {
this.applyAllSymmetry (ms, bsAtoms, this.fmatSupercell != null);
this.fmatSupercell = null;
} else {
this.trimAtomSet ();
}}, "J.adapter.smarter.MSInterface,~S");
Clazz_defineMethod (c$, "setSym", 
 function (pt) {
this.matSupercell.rotate (pt);
this.setSymmetryMinMax (pt);
this.symmetry.toCartesian (pt, false);
return pt;
}, "JU.P3");
Clazz_defineMethod (c$, "setSymmetryMinMax", 
 function (c) {
if (this.rminx > c.x) this.rminx = c.x;
if (this.rminy > c.y) this.rminy = c.y;
if (this.rminz > c.z) this.rminz = c.z;
if (this.rmaxx < c.x) this.rmaxx = c.x;
if (this.rmaxy < c.y) this.rmaxy = c.y;
if (this.rmaxz < c.z) this.rmaxz = c.z;
}, "JU.P3");
Clazz_defineMethod (c$, "isInSymmetryRange", 
 function (c) {
return (c.x >= this.rminx && c.y >= this.rminy && c.z >= this.rminz && c.x <= this.rmaxx && c.y <= this.rmaxy && c.z <= this.rmaxz);
}, "JU.P3");
Clazz_defineMethod (c$, "isWithinCell", 
function (dtype, pt, minX, maxX, minY, maxY, minZ, maxZ, slop) {
return (pt.x > minX - slop && pt.x < maxX + slop && (dtype < 2 || pt.y > minY - slop && pt.y < maxY + slop) && (dtype < 3 || pt.z > minZ - slop && pt.z < maxZ + slop));
}, "~N,JU.P3,~N,~N,~N,~N,~N,~N,~N");
Clazz_defineMethod (c$, "applyAllSymmetry", 
 function (ms, bsAtoms, disableSymmetry) {
if (this.asc.ac == 0) return;
this.noSymmetryCount = (this.asc.baseSymmetryAtomCount == 0 ? this.asc.getLastAtomSetAtomCount () : this.asc.baseSymmetryAtomCount);
this.asc.setTensors ();
this.bondCount0 = this.asc.bondCount;
this.finalizeSymmetry (this.symmetry);
var operationCount = this.symmetry.getSpaceGroupOperationCount ();
this.dtype = Clazz_floatToInt (this.symmetry.getUnitCellInfoType (6));
this.symmetry.setMinMaxLatticeParameters (this.minXYZ, this.maxXYZ);
if (this.doCentroidUnitCell) this.asc.setInfo ("centroidMinMax", [this.minXYZ.x, this.minXYZ.y, this.minXYZ.z, this.maxXYZ.x, this.maxXYZ.y, this.maxXYZ.z, (this.centroidPacked ? 1 : 0)]);
if (this.ptSupercell != null) {
this.asc.setAtomSetAuxiliaryInfo ("supercell", this.ptSupercell);
switch (this.dtype) {
case 3:
this.minXYZ.z *= Clazz_floatToInt (Math.abs (this.ptSupercell.z));
this.maxXYZ.z *= Clazz_floatToInt (Math.abs (this.ptSupercell.z));
case 2:
this.minXYZ.y *= Clazz_floatToInt (Math.abs (this.ptSupercell.y));
this.maxXYZ.y *= Clazz_floatToInt (Math.abs (this.ptSupercell.y));
case 1:
this.minXYZ.x *= Clazz_floatToInt (Math.abs (this.ptSupercell.x));
this.maxXYZ.x *= Clazz_floatToInt (Math.abs (this.ptSupercell.x));
}
}if (this.doCentroidUnitCell || this.doPackUnitCell || this.symmetryRange != 0 && this.maxXYZ.x - this.minXYZ.x == 1 && this.maxXYZ.y - this.minXYZ.y == 1 && this.maxXYZ.z - this.minXYZ.z == 1) {
this.minXYZ0 = JU.P3.new3 (this.minXYZ.x, this.minXYZ.y, this.minXYZ.z);
this.maxXYZ0 = JU.P3.new3 (this.maxXYZ.x, this.maxXYZ.y, this.maxXYZ.z);
if (ms != null) {
ms.setMinMax0 (this.minXYZ0, this.maxXYZ0);
this.minXYZ.set (Clazz_floatToInt (this.minXYZ0.x), Clazz_floatToInt (this.minXYZ0.y), Clazz_floatToInt (this.minXYZ0.z));
this.maxXYZ.set (Clazz_floatToInt (this.maxXYZ0.x), Clazz_floatToInt (this.maxXYZ0.y), Clazz_floatToInt (this.maxXYZ0.z));
}switch (this.dtype) {
case 3:
this.minXYZ.z--;
this.maxXYZ.z++;
case 2:
this.minXYZ.y--;
this.maxXYZ.y++;
case 1:
this.minXYZ.x--;
this.maxXYZ.x++;
}
}var nCells = (this.maxXYZ.x - this.minXYZ.x) * (this.maxXYZ.y - this.minXYZ.y) * (this.maxXYZ.z - this.minXYZ.z);
var cartesianCount = (this.asc.checkSpecial ? this.noSymmetryCount * operationCount * nCells : this.symmetryRange > 0 ? this.noSymmetryCount * operationCount : this.symmetryRange < 0 ? 1 : 1);
var cartesians =  new Array (cartesianCount);
for (var i = 0; i < this.noSymmetryCount; i++) this.asc.atoms[i + this.firstSymmetryAtom].bsSymmetry = JU.BS.newN (operationCount * (nCells + 1));

var pt = 0;
var unitCells =  Clazz_newIntArray (nCells, 0);
this.unitCellTranslations =  new Array (nCells);
var iCell = 0;
var cell555Count = 0;
var absRange = Math.abs (this.symmetryRange);
var checkCartesianRange = (this.symmetryRange != 0);
var checkRangeNoSymmetry = (this.symmetryRange < 0);
var checkRange111 = (this.symmetryRange > 0);
if (checkCartesianRange) {
this.rminx = this.rminy = this.rminz = 3.4028235E38;
this.rmaxx = this.rmaxy = this.rmaxz = -3.4028235E38;
}var symmetry = this.symmetry;
var lastSymmetry = symmetry;
this.latticeOp = symmetry.getLatticeOp ();
this.checkAll = (disableSymmetry || this.asc.atomSetCount == 1 && this.asc.checkSpecial && this.latticeOp >= 0);
this.latticeOnly = (this.asc.checkLatticeOnly && this.latticeOp >= 0);
var pttemp = null;
var op = symmetry.getSpaceGroupOperation (0);
if (this.doPackUnitCell) {
pttemp =  new JU.P3 ();
this.ptOffset.set (0, 0, 0);
}for (var tx = this.minXYZ.x; tx < this.maxXYZ.x; tx++) for (var ty = this.minXYZ.y; ty < this.maxXYZ.y; ty++) for (var tz = this.minXYZ.z; tz < this.maxXYZ.z; tz++) {
this.unitCellTranslations[iCell] = JU.V3.new3 (tx, ty, tz);
unitCells[iCell++] = 555 + tx * 100 + ty * 10 + tz;
if (tx != 0 || ty != 0 || tz != 0 || cartesians.length == 0) continue;
for (pt = 0; pt < this.noSymmetryCount; pt++) {
var atom = this.asc.atoms[this.firstSymmetryAtom + pt];
if (ms != null) {
symmetry = ms.getAtomSymmetry (atom, this.symmetry);
if (symmetry !== lastSymmetry) {
if (symmetry.getSpaceGroupOperationCount () == 0) this.finalizeSymmetry (lastSymmetry = symmetry);
op = symmetry.getSpaceGroupOperation (0);
}}var c = JU.P3.newP (atom);
op.rotTrans (c);
symmetry.toCartesian (c, false);
if (this.doPackUnitCell) {
symmetry.toUnitCell (c, this.ptOffset);
pttemp.setT (c);
symmetry.toFractional (pttemp, false);
if (bsAtoms == null) atom.setT (pttemp);
 else if (atom.distance (pttemp) < 0.0001) bsAtoms.set (atom.index);
 else {
bsAtoms.clear (atom.index);
continue;
}}if (bsAtoms != null) atom.bsSymmetry.clearAll ();
atom.bsSymmetry.set (iCell * operationCount);
atom.bsSymmetry.set (0);
if (checkCartesianRange) this.setSymmetryMinMax (c);
if (pt < cartesianCount) cartesians[pt] = c;
}
if (checkRangeNoSymmetry) {
this.rminx -= absRange;
this.rminy -= absRange;
this.rminz -= absRange;
this.rmaxx += absRange;
this.rmaxy += absRange;
this.rmaxz += absRange;
}cell555Count = pt = this.symmetryAddAtoms (0, 0, 0, 0, pt, iCell * operationCount, cartesians, ms, false);
}


if (checkRange111) {
this.rminx -= absRange;
this.rminy -= absRange;
this.rminz -= absRange;
this.rmaxx += absRange;
this.rmaxy += absRange;
this.rmaxz += absRange;
}iCell = 0;
for (var tx = this.minXYZ.x; tx < this.maxXYZ.x; tx++) for (var ty = this.minXYZ.y; ty < this.maxXYZ.y; ty++) for (var tz = this.minXYZ.z; tz < this.maxXYZ.z; tz++) {
iCell++;
if (tx != 0 || ty != 0 || tz != 0) pt = this.symmetryAddAtoms (tx, ty, tz, cell555Count, pt, iCell * operationCount, cartesians, ms, disableSymmetry);
}


if (iCell * this.noSymmetryCount == this.asc.ac - this.firstSymmetryAtom) this.appendAtomProperties (iCell);
this.setSymmetryOps ();
this.asc.setAtomSetAuxiliaryInfo ("presymmetryAtomIndex", Integer.$valueOf (this.firstSymmetryAtom));
this.asc.setAtomSetAuxiliaryInfo ("presymmetryAtomCount", Integer.$valueOf (this.noSymmetryCount));
this.asc.setAtomSetAuxiliaryInfo ("latticeDesignation", symmetry.getLatticeDesignation ());
this.asc.setAtomSetAuxiliaryInfo ("unitCellRange", unitCells);
this.asc.setAtomSetAuxiliaryInfo ("unitCellTranslations", this.unitCellTranslations);
this.notionalUnitCell =  Clazz_newFloatArray (6, 0);
this.reset ();
}, "J.adapter.smarter.MSInterface,JU.BS,~B");
Clazz_defineMethod (c$, "symmetryAddAtoms", 
 function (transX, transY, transZ, baseCount, pt, iCellOpPt, cartesians, ms, disableSymmetry) {
var isBaseCell = (baseCount == 0);
var addBonds = (this.bondCount0 > this.asc.bondIndex0 && this.applySymmetryToBonds);
var atomMap = (addBonds ?  Clazz_newIntArray (this.noSymmetryCount, 0) : null);
if (this.doPackUnitCell) this.ptOffset.set (transX, transY, transZ);
var range2 = this.symmetryRange * this.symmetryRange;
var checkRangeNoSymmetry = (this.symmetryRange < 0);
var checkRange111 = (this.symmetryRange > 0);
var checkSymmetryMinMax = (isBaseCell && checkRange111);
checkRange111 = new Boolean (checkRange111 & !isBaseCell).valueOf ();
var nOperations = this.symmetry.getSpaceGroupOperationCount ();
var checkSpecial = (nOperations == 1 && !this.doPackUnitCell ? false : this.asc.checkSpecial);
var checkSymmetryRange = (checkRangeNoSymmetry || checkRange111);
var checkDistances = (checkSpecial || checkSymmetryRange);
var addCartesian = (checkSpecial || checkSymmetryMinMax);
var symmetry = this.symmetry;
if (checkRangeNoSymmetry) baseCount = this.noSymmetryCount;
var atomMax = this.firstSymmetryAtom + this.noSymmetryCount;
var ptAtom =  new JU.P3 ();
var code = null;
var subSystemId = '\u0000';
for (var iSym = 0; iSym < nOperations; iSym++) {
if (isBaseCell && iSym == 0 || this.latticeOnly && iSym > 0 && iSym != this.latticeOp) continue;
var pt0 = (disableSymmetry ? 0 : checkSpecial ? pt : checkRange111 ? baseCount : 0);
var spinOp = (this.asc.vibScale == 0 ? symmetry.getSpinOp (iSym) : this.asc.vibScale);
for (var i = this.firstSymmetryAtom; i < atomMax; i++) {
var a = this.asc.atoms[i];
if (a.ignoreSymmetry) continue;
if (this.asc.bsAtoms != null && !this.asc.bsAtoms.get (i)) continue;
if (ms == null) {
symmetry.newSpaceGroupPoint (iSym, a, ptAtom, transX, transY, transZ);
} else {
symmetry = ms.getAtomSymmetry (a, this.symmetry);
symmetry.newSpaceGroupPoint (iSym, a, ptAtom, transX, transY, transZ);
code = symmetry.getSpaceGroupOperationCode (iSym);
if (code != null) {
subSystemId = code.charAt (0);
symmetry = ms.getSymmetryFromCode (code);
if (symmetry.getSpaceGroupOperationCount () == 0) this.finalizeSymmetry (symmetry);
}}var cartesian = JU.P3.newP (ptAtom);
symmetry.toCartesian (cartesian, false);
if (this.doPackUnitCell) {
symmetry.toUnitCell (cartesian, this.ptOffset);
ptAtom.setT (cartesian);
symmetry.toFractional (ptAtom, false);
if (!this.isWithinCell (this.dtype, ptAtom, this.minXYZ0.x, this.maxXYZ0.x, this.minXYZ0.y, this.maxXYZ0.y, this.minXYZ0.z, this.maxXYZ0.z, this.packingError)) continue;
}if (checkSymmetryMinMax) this.setSymmetryMinMax (cartesian);
var special = null;
if (checkDistances) {
var minDist2 = 3.4028235E38;
if (checkSymmetryRange && !this.isInSymmetryRange (cartesian)) continue;
var j0 = (this.checkAll ? this.asc.ac : pt0);
var name = a.atomName;
var id = (code == null ? a.altLoc : subSystemId);
for (var j = j0; --j >= 0; ) {
var pc = cartesians[j];
if (pc == null) continue;
var d2 = cartesian.distanceSquared (pc);
if (checkSpecial && d2 < 0.0001) {
special = this.asc.atoms[this.firstSymmetryAtom + j];
if ((special.atomName == null || special.atomName.equals (name)) && special.altLoc == id) break;
special = null;
}if (checkRange111 && j < baseCount && d2 < minDist2) minDist2 = d2;
}
if (checkRange111 && minDist2 > range2) continue;
}var atomSite = a.atomSite;
if (special != null) {
if (addBonds) atomMap[atomSite] = special.index;
special.bsSymmetry.set (iCellOpPt + iSym);
special.bsSymmetry.set (iSym);
} else {
if (addBonds) atomMap[atomSite] = this.asc.ac;
var atom1 = this.asc.newCloneAtom (a);
if (this.asc.bsAtoms != null) this.asc.bsAtoms.set (atom1.index);
atom1.setT (ptAtom);
if (spinOp != 0 && atom1.vib != null) {
symmetry.getSpaceGroupOperation (iSym).rotate (atom1.vib);
atom1.vib.scale (spinOp);
System.out.println ("vib for iSym " + iSym + " " + atom1 + " " + atom1.vib);
}atom1.atomSite = atomSite;
if (code != null) atom1.altLoc = subSystemId;
atom1.bsSymmetry = JU.BSUtil.newAndSetBit (iCellOpPt + iSym);
atom1.bsSymmetry.set (iSym);
if (addCartesian) cartesians[pt++] = cartesian;
var tensors = a.tensors;
if (tensors != null) {
atom1.tensors = null;
for (var j = tensors.size (); --j >= 0; ) {
var t = tensors.get (j);
if (t == null) continue;
if (nOperations == 1) atom1.addTensor (t.copyTensor (), null, false);
 else this.addRotatedTensor (atom1, t, iSym, false, symmetry);
}
}}}
if (addBonds) {
var bonds = this.asc.bonds;
var atoms = this.asc.atoms;
for (var bondNum = this.asc.bondIndex0; bondNum < this.bondCount0; bondNum++) {
var bond = bonds[bondNum];
var atom1 = atoms[bond.atomIndex1];
var atom2 = atoms[bond.atomIndex2];
if (atom1 == null || atom2 == null) continue;
var iAtom1 = atomMap[atom1.atomSite];
var iAtom2 = atomMap[atom2.atomSite];
if (iAtom1 >= atomMax || iAtom2 >= atomMax) this.asc.addNewBondWithOrder (iAtom1, iAtom2, bond.order);
}
}}
return pt;
}, "~N,~N,~N,~N,~N,~N,~A,J.adapter.smarter.MSInterface,~B");
Clazz_defineMethod (c$, "appendAtomProperties", 
 function (nTimes) {
var p = this.asc.getAtomSetAuxiliaryInfoValue (-1, "atomProperties");
if (p == null) {
return;
}for (var entry, $entry = p.entrySet ().iterator (); $entry.hasNext () && ((entry = $entry.next ()) || true);) {
var key = entry.getKey ();
var data = entry.getValue ();
var s =  new JU.SB ();
for (var i = nTimes; --i >= 0; ) s.append (data);

p.put (key, s.toString ());
}
}, "~N");
Clazz_defineMethod (c$, "finalizeSymmetry", 
 function (symmetry) {
var name = this.asc.getAtomSetAuxiliaryInfoValue (-1, "spaceGroup");
symmetry.setFinalOperations (name, this.asc.atoms, this.firstSymmetryAtom, this.noSymmetryCount, this.doNormalize, this.filterSymop);
if (this.filterSymop != null || name == null || name.equals ("unspecified!")) this.asc.setAtomSetSpaceGroupName (symmetry.getSpaceGroupName ());
}, "J.api.SymmetryInterface");
Clazz_defineMethod (c$, "setSymmetryOps", 
 function () {
var operationCount = this.symmetry.getSpaceGroupOperationCount ();
if (operationCount > 0) {
var symmetryList =  new Array (operationCount);
for (var i = 0; i < operationCount; i++) symmetryList[i] = "" + this.symmetry.getSpaceGroupXyz (i, this.doNormalize);

this.asc.setAtomSetAuxiliaryInfo ("symmetryOperations", symmetryList);
this.asc.setAtomSetAuxiliaryInfo ("symmetryOps", this.symmetry.getSymmetryOperations ());
}this.asc.setAtomSetAuxiliaryInfo ("symmetryCount", Integer.$valueOf (operationCount));
});
Clazz_defineMethod (c$, "applySymmetryBio", 
function (thisBiomolecule, notionalUnitCell, applySymmetryToBonds, filter) {
if (this.latticeCells != null && this.latticeCells[0] != 0) {
JU.Logger.error ("Cannot apply biomolecule when lattice cells are indicated");
return;
}var particleMode = (filter.indexOf ("BYCHAIN") >= 0 ? 1 : filter.indexOf ("BYSYMOP") >= 0 ? 2 : 0);
this.doNormalize = false;
var biomts = thisBiomolecule.get ("biomts");
if (biomts.size () < 2) return;
this.symmetry = null;
if (!Float.isNaN (notionalUnitCell[0])) this.setNotionalUnitCell (notionalUnitCell, null, this.unitCellOffset);
this.getSymmetry ().setSpaceGroup (this.doNormalize);
this.addSpaceGroupOperation (null, "x,y,z");
var name = thisBiomolecule.get ("name");
this.asc.setAtomSetSpaceGroupName (name);
var len = biomts.size ();
this.applySymmetryToBonds = applySymmetryToBonds;
this.bondCount0 = this.asc.bondCount;
var addBonds = (this.bondCount0 > this.asc.bondIndex0 && applySymmetryToBonds);
var atomMap = (addBonds ?  Clazz_newIntArray (this.asc.ac, 0) : null);
this.firstSymmetryAtom = this.asc.getLastAtomSetAtomIndex ();
var atomMax = this.asc.ac;
var ht =  new java.util.Hashtable ();
var nChain = 0;
var atoms = this.asc.atoms;
switch (particleMode) {
case 1:
for (var i = atomMax; --i >= this.firstSymmetryAtom; ) {
var id = Integer.$valueOf (atoms[i].chainID);
var bs = ht.get (id);
if (bs == null) {
nChain++;
ht.put (id, bs =  new JU.BS ());
}bs.set (i);
}
this.asc.bsAtoms =  new JU.BS ();
for (var i = 0; i < nChain; i++) {
this.asc.bsAtoms.set (atomMax + i);
var a =  new J.adapter.smarter.Atom ();
a.set (0, 0, 0);
a.radius = 16;
this.asc.addAtom (a);
}
var ichain = 0;
for (var e, $e = ht.entrySet ().iterator (); $e.hasNext () && ((e = $e.next ()) || true);) {
var a = atoms[atomMax + ichain++];
var bs = e.getValue ();
for (var i = bs.nextSetBit (0); i >= 0; i = bs.nextSetBit (i + 1)) a.add (atoms[i]);

a.scale (1 / bs.cardinality ());
a.atomName = "Pt" + ichain;
a.chainID = e.getKey ().intValue ();
}
this.firstSymmetryAtom = atomMax;
atomMax += nChain;
break;
case 2:
this.asc.bsAtoms =  new JU.BS ();
this.asc.bsAtoms.set (atomMax);
var a = atoms[atomMax] =  new J.adapter.smarter.Atom ();
a.set (0, 0, 0);
for (var i = atomMax; --i >= this.firstSymmetryAtom; ) a.add (atoms[i]);

a.scale (1 / (atomMax - this.firstSymmetryAtom));
a.atomName = "Pt";
a.radius = 16;
this.firstSymmetryAtom = atomMax++;
break;
}
if (filter.indexOf ("#<") >= 0) {
len = Math.min (len, JU.PT.parseInt (filter.substring (filter.indexOf ("#<") + 2)) - 1);
filter = JU.PT.rep (filter, "#<", "_<");
}for (var iAtom = this.firstSymmetryAtom; iAtom < atomMax; iAtom++) atoms[iAtom].bsSymmetry = JU.BSUtil.newAndSetBit (0);

for (var i = 1; i < len; i++) {
if (filter.indexOf ("!#") >= 0) {
if (filter.indexOf ("!#" + (i + 1) + ";") >= 0) continue;
} else if (filter.indexOf ("#") >= 0 && filter.indexOf ("#" + (i + 1) + ";") < 0) {
continue;
}var mat = biomts.get (i);
for (var iAtom = this.firstSymmetryAtom; iAtom < atomMax; iAtom++) {
if (this.asc.bsAtoms != null && !this.asc.bsAtoms.get (iAtom)) continue;
try {
var atomSite = atoms[iAtom].atomSite;
var atom1;
if (addBonds) atomMap[atomSite] = this.asc.ac;
atom1 = this.asc.newCloneAtom (atoms[iAtom]);
if (this.asc.bsAtoms != null) this.asc.bsAtoms.set (atom1.index);
atom1.atomSite = atomSite;
mat.rotTrans (atom1);
atom1.bsSymmetry = JU.BSUtil.newAndSetBit (i);
if (addBonds) {
for (var bondNum = this.asc.bondIndex0; bondNum < this.bondCount0; bondNum++) {
var bond = this.asc.bonds[bondNum];
var iAtom1 = atomMap[atoms[bond.atomIndex1].atomSite];
var iAtom2 = atomMap[atoms[bond.atomIndex2].atomSite];
if (iAtom1 >= atomMax || iAtom2 >= atomMax) this.asc.addNewBondWithOrder (iAtom1, iAtom2, bond.order);
}
}} catch (e) {
if (Clazz_exceptionOf (e, Exception)) {
this.asc.errorMessage = "appendAtomCollection error: " + e;
} else {
throw e;
}
}
}
if (i > 0) this.symmetry.addBioMoleculeOperation (mat, false);
}
this.noSymmetryCount = atomMax - this.firstSymmetryAtom;
this.asc.setAtomSetAuxiliaryInfo ("presymmetryAtomIndex", Integer.$valueOf (this.firstSymmetryAtom));
this.asc.setAtomSetAuxiliaryInfo ("presymmetryAtomCount", Integer.$valueOf (this.noSymmetryCount));
this.asc.setAtomSetAuxiliaryInfo ("biosymmetryCount", Integer.$valueOf (len));
this.asc.setAtomSetAuxiliaryInfo ("biosymmetry", this.symmetry);
this.finalizeSymmetry (this.symmetry);
this.setSymmetryOps ();
this.reset ();
}, "java.util.Map,~A,~B,~S");
Clazz_defineMethod (c$, "reset", 
 function () {
this.asc.coordinatesAreFractional = false;
this.asc.setAtomSetAuxiliaryInfo ("hasSymmetry", Boolean.TRUE);
this.asc.setGlobalBoolean (1);
});
Clazz_defineMethod (c$, "addRotatedTensor", 
function (a, t, iSym, reset, symmetry) {
if (this.ptTemp == null) {
this.ptTemp =  new JU.P3 ();
this.mTemp =  new JU.M3 ();
}return a.addTensor ((J.api.Interface.getUtil ("Tensor")).setFromEigenVectors (symmetry.rotateAxes (iSym, t.eigenVectors, this.ptTemp, this.mTemp), t.eigenValues, t.isIsotropic ? "iso" : t.type, t.id, t), null, reset);
}, "J.adapter.smarter.Atom,JU.Tensor,~N,~B,J.api.SymmetryInterface");
Clazz_defineMethod (c$, "setTensors", 
function () {
var n = this.asc.ac;
for (var i = this.asc.getLastAtomSetAtomIndex (); i < n; i++) {
var a = this.asc.atoms[i];
if (a.anisoBorU == null) continue;
a.addTensor (this.symmetry.getTensor (a.anisoBorU), null, false);
if (Float.isNaN (a.bfactor)) a.bfactor = a.anisoBorU[7] * 100;
a.anisoBorU = null;
}
});
Clazz_defineMethod (c$, "setTimeReversal", 
function (op, timeRev) {
this.symmetry.setTimeReversal (op, timeRev);
}, "~N,~N");
Clazz_defineMethod (c$, "rotateToSuperCell", 
function (t) {
if (this.matSupercell != null) this.matSupercell.rotate (t);
}, "JU.V3");
Clazz_defineMethod (c$, "setVibVectors", 
function () {
if (this.nVib > 0 || this.asc.iSet < 0 || !this.vibsFractional) return this.nVib;
var i0 = this.asc.getAtomSetAtomIndex (this.asc.iSet);
for (var i = this.asc.ac; --i >= i0; ) {
if (this.asc.atoms[i].vib != null) {
var v =  new JU.Vibration ();
v.setT (this.asc.atoms[i].vib);
v.modDim = -2;
this.symmetry.toCartesian (v, true);
this.asc.atoms[i].vib = v;
this.nVib++;
}}
return this.nVib;
});
Clazz_defineMethod (c$, "trimAtomSet", 
 function () {
this.symmetry.setMinMaxLatticeParameters (this.minXYZ, this.maxXYZ);
var bs = this.asc.bsAtoms;
var atoms = this.asc.atoms;
if (bs == null) bs = this.asc.bsAtoms = JU.BSUtil.newBitSet2 (0, this.asc.ac);
for (var i = bs.nextSetBit (0); i >= 0; i = bs.nextSetBit (i + 1)) {
if (!this.isWithinCell (this.dtype, atoms[i], this.minXYZ.x, this.maxXYZ.x, this.minXYZ.y, this.maxXYZ.y, this.minXYZ.z, this.maxXYZ.z, this.packingError)) {
bs.clear (i);
}}
});
Clazz_defineStatics (c$,
"PARTICLE_NONE", 0,
"PARTICLE_CHAIN", 1,
"PARTICLE_SYMOP", 2);
});
Clazz_declarePackage ("J.api");
Clazz_declareInterface (J.api, "SymmetryInterface");
Clazz_declarePackage ("JS");
Clazz_load (["J.api.SymmetryInterface"], "JS.Symmetry", ["java.lang.Float", "java.util.Hashtable", "JU.BS", "$.Lst", "$.M3", "$.M4", "$.P3", "$.SB", "$.V3", "JM.Atom", "JS.PointGroup", "$.SpaceGroup", "$.SymmetryInfo", "$.SymmetryOperation", "$.UnitCell", "JU.Escape", "$.Logger", "$.SimpleUnitCell"], function () {
c$ = Clazz_decorateAsClass (function () {
this.pointGroup = null;
this.spaceGroup = null;
this.symmetryInfo = null;
this.unitCell = null;
this.$isBio = false;
Clazz_instantialize (this, arguments);
}, JS, "Symmetry", null, J.api.SymmetryInterface);
Clazz_defineMethod (c$, "isBio", 
function () {
return this.$isBio;
});
Clazz_makeConstructor (c$, 
function () {
});
Clazz_overrideMethod (c$, "setPointGroup", 
function (siLast, atomset, bsAtoms, haveVibration, distanceTolerance, linearTolerance) {
this.pointGroup = JS.PointGroup.getPointGroup (siLast == null ? null : (siLast).pointGroup, atomset, bsAtoms, haveVibration, distanceTolerance, linearTolerance);
return this;
}, "J.api.SymmetryInterface,~A,JU.BS,~B,~N,~N");
Clazz_overrideMethod (c$, "getPointGroupName", 
function () {
return this.pointGroup.getName ();
});
Clazz_overrideMethod (c$, "getPointGroupInfo", 
function (modelIndex, asDraw, asInfo, type, index, scale) {
if (!asDraw && !asInfo && this.pointGroup.textInfo != null) return this.pointGroup.textInfo;
 else if (asDraw && this.pointGroup.isDrawType (type, index, scale)) return this.pointGroup.drawInfo;
 else if (asInfo && this.pointGroup.info != null) return this.pointGroup.info;
return this.pointGroup.getInfo (modelIndex, asDraw, asInfo, type, index, scale);
}, "~N,~B,~B,~S,~N,~N");
Clazz_defineMethod (c$, "setSpaceGroup", 
function (doNormalize) {
if (this.spaceGroup == null) this.spaceGroup = (JS.SpaceGroup.getNull (true)).set (doNormalize);
}, "~B");
Clazz_defineMethod (c$, "addSpaceGroupOperation", 
function (xyz, opId) {
return this.spaceGroup.addSymmetry (xyz, opId, false);
}, "~S,~N");
Clazz_defineMethod (c$, "addBioMoleculeOperation", 
function (mat, isReverse) {
this.$isBio = this.spaceGroup.isBio = true;
return this.spaceGroup.addSymmetry ((isReverse ? "!" : "") + "[[bio" + mat, 0, false);
}, "JU.M4,~B");
Clazz_overrideMethod (c$, "setLattice", 
function (latt) {
this.spaceGroup.setLatticeParam (latt);
}, "~N");
Clazz_defineMethod (c$, "getSpaceGroup", 
function () {
return this.spaceGroup;
});
Clazz_overrideMethod (c$, "setSpaceGroupFrom", 
function (symmetry) {
this.spaceGroup = symmetry.getSpaceGroup ();
}, "J.api.SymmetryInterface");
Clazz_overrideMethod (c$, "createSpaceGroup", 
function (desiredSpaceGroupIndex, name, object) {
this.spaceGroup = JS.SpaceGroup.createSpaceGroup (desiredSpaceGroupIndex, name, object);
if (this.spaceGroup != null && JU.Logger.debugging) JU.Logger.debug ("using generated space group " + this.spaceGroup.dumpInfo (null));
return this.spaceGroup != null;
}, "~N,~S,~O");
Clazz_overrideMethod (c$, "getSpaceGroupInfoStr", 
function (name, cellInfo) {
return JS.SpaceGroup.getInfo (name, cellInfo);
}, "~S,J.api.SymmetryInterface");
Clazz_overrideMethod (c$, "getLatticeDesignation", 
function () {
return this.spaceGroup.getLatticeDesignation ();
});
Clazz_overrideMethod (c$, "setFinalOperations", 
function (name, atoms, iAtomFirst, noSymmetryCount, doNormalize, filterSymop) {
if (name != null && (name.startsWith ("bio") || name.indexOf (" *(") >= 0)) this.spaceGroup.name = name;
if (filterSymop != null) {
var lst =  new JU.Lst ();
lst.addLast (this.spaceGroup.operations[0]);
for (var i = 1; i < this.spaceGroup.operationCount; i++) if (filterSymop.contains (" " + (i + 1) + " ")) lst.addLast (this.spaceGroup.operations[i]);

this.spaceGroup = JS.SpaceGroup.createSpaceGroup (-1, name + " *(" + filterSymop.trim () + ")", lst);
}this.spaceGroup.setFinalOperations (atoms, iAtomFirst, noSymmetryCount, doNormalize);
}, "~S,~A,~N,~N,~B,~S");
Clazz_defineMethod (c$, "getSpaceGroupOperation", 
function (i) {
return (i >= this.spaceGroup.operations.length ? null : this.spaceGroup.finalOperations == null ? this.spaceGroup.operations[i] : this.spaceGroup.finalOperations[i]);
}, "~N");
Clazz_overrideMethod (c$, "getSpaceGroupXyz", 
function (i, doNormalize) {
return this.spaceGroup.getXyz (i, doNormalize);
}, "~N,~B");
Clazz_defineMethod (c$, "newSpaceGroupPoint", 
function (i, atom1, atom2, transX, transY, transZ) {
if (this.spaceGroup.finalOperations == null) {
if (!this.spaceGroup.operations[i].isFinalized) this.spaceGroup.operations[i].doFinalize ();
this.spaceGroup.operations[i].newPoint (atom1, atom2, transX, transY, transZ);
return;
}this.spaceGroup.finalOperations[i].newPoint (atom1, atom2, transX, transY, transZ);
}, "~N,JU.P3,JU.P3,~N,~N,~N");
Clazz_overrideMethod (c$, "rotateAxes", 
function (iop, axes, ptTemp, mTemp) {
return (iop == 0 ? axes : this.spaceGroup.finalOperations[iop].rotateAxes (axes, this.unitCell, ptTemp, mTemp));
}, "~N,~A,JU.P3,JU.M3");
Clazz_overrideMethod (c$, "getSpaceGroupOperationCode", 
function (iOp) {
return this.spaceGroup.operations[iOp].subsystemCode;
}, "~N");
Clazz_overrideMethod (c$, "setTimeReversal", 
function (op, val) {
this.spaceGroup.operations[op].setTimeReversal (val);
}, "~N,~N");
Clazz_overrideMethod (c$, "getSpinOp", 
function (op) {
return this.spaceGroup.operations[op].getSpinOp ();
}, "~N");
Clazz_overrideMethod (c$, "addLatticeVectors", 
function (lattvecs) {
return this.spaceGroup.addLatticeVectors (lattvecs);
}, "JU.Lst");
Clazz_overrideMethod (c$, "getLatticeOp", 
function () {
return this.spaceGroup.latticeOp;
});
Clazz_overrideMethod (c$, "getOperationRsVs", 
function (iop) {
return (this.spaceGroup.finalOperations == null ? this.spaceGroup.operations : this.spaceGroup.finalOperations)[iop].rsvs;
}, "~N");
Clazz_overrideMethod (c$, "getSiteMultiplicity", 
function (pt) {
return this.spaceGroup.getSiteMultiplicity (pt, this.unitCell);
}, "JU.P3");
Clazz_overrideMethod (c$, "addOp", 
function (code, rs, vs, sigma) {
this.spaceGroup.isSSG = true;
var s = JS.SymmetryOperation.getXYZFromRsVs (rs, vs, false);
var i = this.spaceGroup.addSymmetry (s, -1, true);
this.spaceGroup.operations[i].setSigma (code, sigma);
return s;
}, "~S,JU.Matrix,JU.Matrix,JU.Matrix");
Clazz_overrideMethod (c$, "getMatrixFromString", 
function (xyz, rotTransMatrix, allowScaling, modDim) {
return JS.SymmetryOperation.getMatrixFromString (null, xyz, rotTransMatrix, allowScaling);
}, "~S,~A,~B,~N");
Clazz_defineMethod (c$, "getSpaceGroupName", 
function () {
return (this.symmetryInfo != null ? this.symmetryInfo.sgName : this.spaceGroup != null ? this.spaceGroup.getName () : this.unitCell != null && this.unitCell.name.length > 0 ? "cell=" + this.unitCell.name : "");
});
Clazz_overrideMethod (c$, "getSpaceGroupOperationCount", 
function () {
return (this.symmetryInfo != null ? this.symmetryInfo.symmetryOperations.length : this.spaceGroup != null && this.spaceGroup.finalOperations != null ? this.spaceGroup.finalOperations.length : 0);
});
Clazz_overrideMethod (c$, "getCoordinatesAreFractional", 
function () {
return this.symmetryInfo == null || this.symmetryInfo.coordinatesAreFractional;
});
Clazz_overrideMethod (c$, "getCellRange", 
function () {
return this.symmetryInfo.cellRange;
});
Clazz_overrideMethod (c$, "getSymmetryInfoStr", 
function () {
return this.symmetryInfo.infoStr;
});
Clazz_defineMethod (c$, "getSymmetryOperations", 
function () {
return this.symmetryInfo == null ? this.spaceGroup.finalOperations : this.symmetryInfo.symmetryOperations;
});
Clazz_overrideMethod (c$, "isPeriodic", 
function () {
return (this.symmetryInfo == null ? false : this.symmetryInfo.isPeriodic ());
});
Clazz_overrideMethod (c$, "setSymmetryInfo", 
function (modelIndex, modelAuxiliaryInfo, notionalCell) {
this.symmetryInfo =  new JS.SymmetryInfo ();
var notionalUnitcell = this.symmetryInfo.setSymmetryInfo (modelAuxiliaryInfo, notionalCell);
if (notionalUnitcell == null) return;
this.setUnitCell (notionalUnitcell, modelAuxiliaryInfo.containsKey ("jmolData"));
this.unitCell.moreInfo = modelAuxiliaryInfo.get ("moreUnitCellInfo");
modelAuxiliaryInfo.put ("infoUnitCell", this.getUnitCellAsArray (false));
this.setOffsetPt (modelAuxiliaryInfo.get ("unitCellOffset"));
var matUnitCellOrientation = modelAuxiliaryInfo.get ("matUnitCellOrientation");
if (matUnitCellOrientation != null) this.initializeOrientation (matUnitCellOrientation);
if (JU.Logger.debugging) JU.Logger.debug ("symmetryInfos[" + modelIndex + "]:\n" + this.unitCell.dumpInfo (true));
}, "~N,java.util.Map,~A");
Clazz_overrideMethod (c$, "haveUnitCell", 
function () {
return (this.unitCell != null);
});
Clazz_overrideMethod (c$, "checkUnitCell", 
function (uc, cell, ptTemp, isAbsolute) {
uc.toFractional (ptTemp, isAbsolute);
var slop = 0.02;
return (ptTemp.x >= cell.x - 1 - slop && ptTemp.x <= cell.x + slop && ptTemp.y >= cell.y - 1 - slop && ptTemp.y <= cell.y + slop && ptTemp.z >= cell.z - 1 - slop && ptTemp.z <= cell.z + slop);
}, "J.api.SymmetryInterface,JU.P3,JU.P3,~B");
Clazz_defineMethod (c$, "setUnitCell", 
function (notionalUnitCell, setRelative) {
this.unitCell = JS.UnitCell.newA (notionalUnitCell, setRelative);
}, "~A,~B");
Clazz_overrideMethod (c$, "unitCellEquals", 
function (uc2) {
return ((uc2)).unitCell.isSameAs (this.unitCell);
}, "J.api.SymmetryInterface");
Clazz_overrideMethod (c$, "getUnitCellState", 
function () {
return (this.unitCell == null ? "" : this.unitCell.getState ());
});
Clazz_overrideMethod (c$, "getMoreInfo", 
function () {
return this.unitCell.moreInfo;
});
Clazz_defineMethod (c$, "getUnitsymmetryInfo", 
function () {
return this.unitCell.dumpInfo (false);
});
Clazz_overrideMethod (c$, "initializeOrientation", 
function (mat) {
this.unitCell.initOrientation (mat);
}, "JU.M3");
Clazz_overrideMethod (c$, "unitize", 
function (ptFrac) {
this.unitCell.unitize (ptFrac);
}, "JU.P3");
Clazz_overrideMethod (c$, "toUnitCell", 
function (pt, offset) {
this.unitCell.toUnitCell (pt, offset);
}, "JU.P3,JU.P3");
Clazz_defineMethod (c$, "toCartesian", 
function (fpt, ignoreOffset) {
if (!this.$isBio) this.unitCell.toCartesian (fpt, ignoreOffset);
}, "JU.T3,~B");
Clazz_overrideMethod (c$, "toSupercell", 
function (fpt) {
return this.unitCell.toSupercell (fpt);
}, "JU.P3");
Clazz_defineMethod (c$, "toFractional", 
function (pt, isAbsolute) {
if (!this.$isBio) this.unitCell.toFractional (pt, isAbsolute);
}, "JU.T3,~B");
Clazz_defineMethod (c$, "getNotionalUnitCell", 
function () {
return this.unitCell.getNotionalUnitCell ();
});
Clazz_overrideMethod (c$, "getUnitCellAsArray", 
function (vectorsOnly) {
return this.unitCell.getUnitCellAsArray (vectorsOnly);
}, "~B");
Clazz_overrideMethod (c$, "getTensor", 
function (parBorU) {
if (parBorU == null) return null;
if (this.unitCell == null) this.unitCell = JS.UnitCell.newA ([1, 1, 1, 90, 90, 90], true);
return this.unitCell.getTensor (parBorU);
}, "~A");
Clazz_overrideMethod (c$, "getUnitCellVertices", 
function () {
return this.unitCell.getVertices ();
});
Clazz_overrideMethod (c$, "getCartesianOffset", 
function () {
return this.unitCell.getCartesianOffset ();
});
Clazz_overrideMethod (c$, "getFractionalOffset", 
function () {
return this.unitCell.getFractionalOffset ();
});
Clazz_overrideMethod (c$, "setOffsetPt", 
function (pt) {
this.unitCell.setOffset (pt);
}, "JU.T3");
Clazz_overrideMethod (c$, "setOffset", 
function (nnn) {
var pt =  new JU.P3 ();
JU.SimpleUnitCell.ijkToPoint3f (nnn, pt, 0);
this.unitCell.setOffset (pt);
}, "~N");
Clazz_overrideMethod (c$, "getUnitCellMultiplier", 
function () {
return this.unitCell.getUnitCellMultiplier ();
});
Clazz_overrideMethod (c$, "getCanonicalCopy", 
function (scale, withOffset) {
return this.unitCell.getCanonicalCopy (scale, withOffset);
}, "~N,~B");
Clazz_overrideMethod (c$, "getUnitCellInfoType", 
function (infoType) {
return this.unitCell.getInfo (infoType);
}, "~N");
Clazz_overrideMethod (c$, "getUnitCellInfo", 
function () {
return this.unitCell.dumpInfo (false);
});
Clazz_overrideMethod (c$, "isSlab", 
function () {
return this.unitCell.isSlab ();
});
Clazz_overrideMethod (c$, "isPolymer", 
function () {
return this.unitCell.isPolymer ();
});
Clazz_overrideMethod (c$, "setMinMaxLatticeParameters", 
function (minXYZ, maxXYZ) {
this.unitCell.setMinMaxLatticeParameters (minXYZ, maxXYZ);
}, "JU.P3i,JU.P3i");
Clazz_overrideMethod (c$, "checkDistance", 
function (f1, f2, distance, dx, iRange, jRange, kRange, ptOffset) {
return this.unitCell.checkDistance (f1, f2, distance, dx, iRange, jRange, kRange, ptOffset);
}, "JU.P3,JU.P3,~N,~N,~N,~N,~N,JU.P3");
Clazz_overrideMethod (c$, "getUnitCellVectors", 
function () {
return this.unitCell.getUnitCellVectors ();
});
Clazz_overrideMethod (c$, "getUnitCell", 
function (points, setRelative, name) {
this.unitCell = JS.UnitCell.newP (points, setRelative);
if (name != null) this.unitCell.name = name;
return this;
}, "~A,~B,~S");
Clazz_overrideMethod (c$, "isSupercell", 
function () {
return this.unitCell.isSupercell ();
});
Clazz_overrideMethod (c$, "notInCentroid", 
function (modelSet, bsAtoms, minmax) {
try {
var bsDelete =  new JU.BS ();
var iAtom0 = bsAtoms.nextSetBit (0);
var molecules = modelSet.getMolecules ();
var moleculeCount = molecules.length;
var atoms = modelSet.at;
var isOneMolecule = (molecules[moleculeCount - 1].firstAtomIndex == modelSet.am[atoms[iAtom0].mi].firstAtomIndex);
var center =  new JU.P3 ();
var centroidPacked = (minmax[6] == 1);
nextMol : for (var i = moleculeCount; --i >= 0 && bsAtoms.get (molecules[i].firstAtomIndex); ) {
var bs = molecules[i].atomList;
center.set (0, 0, 0);
var n = 0;
for (var j = bs.nextSetBit (0); j >= 0; j = bs.nextSetBit (j + 1)) {
if (isOneMolecule || centroidPacked) {
center.setT (atoms[j]);
if (this.isNotCentroid (center, 1, minmax, centroidPacked)) {
if (isOneMolecule) bsDelete.set (j);
} else if (!isOneMolecule) {
continue nextMol;
}} else {
center.add (atoms[j]);
n++;
}}
if (centroidPacked || n > 0 && this.isNotCentroid (center, n, minmax, false)) bsDelete.or (bs);
}
return bsDelete;
} catch (e) {
if (Clazz_exceptionOf (e, Exception)) {
return null;
} else {
throw e;
}
}
}, "JM.ModelSet,JU.BS,~A");
Clazz_defineMethod (c$, "isNotCentroid", 
 function (center, n, minmax, centroidPacked) {
center.scale (1 / n);
this.toFractional (center, false);
if (centroidPacked) return (center.x + 0.000005 <= minmax[0] || center.x - 0.000005 > minmax[3] || center.y + 0.000005 <= minmax[1] || center.y - 0.000005 > minmax[4] || center.z + 0.000005 <= minmax[2] || center.z - 0.000005 > minmax[5]);
return (center.x + 0.000005 <= minmax[0] || center.x + 0.00005 > minmax[3] || center.y + 0.000005 <= minmax[1] || center.y + 0.00005 > minmax[4] || center.z + 0.000005 <= minmax[2] || center.z + 0.00005 > minmax[5]);
}, "JU.P3,~N,~A,~B");
Clazz_overrideMethod (c$, "getSpaceGroupInfo", 
function (modelSet, modelIndex, sgName, symOp, pt1, pt2, drawID, type) {
var strOperations = null;
var info = null;
var cellInfo = null;
var infolist = null;
var isStandard = true;
if (sgName == null) {
if (modelIndex <= 0) modelIndex = (Clazz_instanceOf (pt1, JM.Atom) ? (pt1).mi : modelSet.vwr.am.cmi);
var isBio = false;
if (modelIndex < 0) strOperations = "no single current model";
 else if (!(isBio = (cellInfo = modelSet.am[modelIndex].biosymmetry) != null) && (cellInfo = modelSet.getUnitCell (modelIndex)) == null) strOperations = "not applicable";
if (strOperations != null) {
info =  new java.util.Hashtable ();
info.put ("spaceGroupInfo", strOperations);
info.put ("symmetryInfo", "");
} else if (pt1 == null && drawID == null && symOp != 0) {
info = modelSet.getInfo (modelIndex, "spaceGroupInfo");
}if (info != null) return info;
info =  new java.util.Hashtable ();
if (pt1 == null && drawID == null && symOp == 0) modelSet.setInfo (modelIndex, "spaceGroupInfo", info);
sgName = cellInfo.getSpaceGroupName ();
var ops = cellInfo.getSymmetryOperations ();
var sg = (isBio ? (cellInfo).spaceGroup : null);
var jf = "";
if (ops == null) {
strOperations = "\n no symmetry operations";
} else {
isStandard = !isBio;
if (isBio) this.spaceGroup = (JS.SpaceGroup.getNull (false)).set (false);
 else this.setSpaceGroup (false);
strOperations = "\n" + ops.length + " symmetry operations:";
infolist =  new Array (ops.length);
var centering = null;
for (var i = 0; i < ops.length; i++) {
var op = (ops[i]);
var xyz = op.xyz;
var iop = (isBio ? this.addBioMoleculeOperation (sg.finalOperations[i], false) : this.addSpaceGroupOperation ("=" + xyz, i + 1));
if (iop < 0) continue;
op = this.getSpaceGroupOperation (i);
if (op.timeReversal != 0 || op.modDim > 0) isStandard = false;
centering = op.setCentering (centering, false);
jf += ";" + xyz;
infolist[i] = (symOp > 0 && symOp - 1 != iop ? null : op.getDescription (modelSet, cellInfo, pt1, pt2, drawID));
if (infolist[i] != null) strOperations += "\n" + (i + 1) + "\t" + infolist[i][0] + "\t" + infolist[i][2];
}
}jf = jf.substring (jf.indexOf (";") + 1);
if (sgName.indexOf ("[--]") >= 0) sgName = jf;
} else {
info =  new java.util.Hashtable ();
}info.put ("spaceGroupName", sgName);
if (infolist != null) {
info.put ("operations", infolist);
info.put ("symmetryInfo", strOperations);
}var data;
if (isStandard) {
data = this.getSpaceGroupInfoStr (sgName, cellInfo);
if (data == null || data.equals ("?")) data = "could not identify space group from name: " + sgName + "\nformat: show spacegroup \"2\" or \"P 2c\" " + "or \"C m m m\" or \"x, y, z;-x ,-y, -z\"";
} else {
data = sgName;
}info.put ("spaceGroupInfo", data);
return info;
}, "JM.ModelSet,~N,~S,~N,JU.P3,JU.P3,~S,~S");
Clazz_overrideMethod (c$, "getSymmetryInfoString", 
function (modelSet, modelIndex, symOp, pt1, pt2, drawID, type) {
var sginfo = this.getSpaceGroupInfo (modelSet, modelIndex, null, symOp, pt1, pt2, drawID, type);
if (sginfo == null) return "";
var labelOnly = "label".equals (type);
var prettyMat = "fmatrix".equals (type);
var infolist = sginfo.get ("operations");
if (infolist == null) return "";
var sb =  new JU.SB ();
symOp--;
for (var i = 0; i < infolist.length; i++) {
if (infolist[i] == null || symOp >= 0 && symOp != i) continue;
if (drawID != null) return infolist[i][3];
if (sb.length () > 0) sb.appendC ('\n');
if (prettyMat) {
sb.append (JS.SymmetryOperation.cleanMatrix (infolist[i][10])).append ("\t");
} else if (!labelOnly) {
if (symOp < 0) sb.appendI (i + 1).append ("\t");
sb.append (infolist[i][0]).append ("\t");
}sb.append (infolist[i][2]);
}
if (sb.length () == 0 && drawID != null) sb.append ("draw " + drawID + "* delete");
return sb.toString ();
}, "JM.ModelSet,~N,~N,JU.P3,JU.P3,~S,~S");
Clazz_overrideMethod (c$, "getSymmetryInfo", 
function (modelSet, iModel, iAtom, uc, xyz, op, pt, pt2, id, type) {
if (pt2 != null) return this.getSymmetryInfoString (modelSet, iModel, op, pt, pt2, (id == null ? "sym" : id), type == 1826248716 ? "label" : null);
var isBio = uc.isBio ();
var sym = uc;
var iop = op;
var centering = null;
if (xyz == null) {
var ops = sym.getSymmetryOperations ();
if (ops == null || op == 0 || Math.abs (op) > ops.length) {
return (type == 135176 ? "draw ID sym_* delete" : "");
}if (op > 0) {
xyz = ops[iop = op - 1].xyz;
} else {
xyz = ops[iop = -1 - op].xyz;
}centering = ops[iop].centering;
} else {
iop = op = 0;
}var symTemp = modelSet.getSymTemp (true);
symTemp.setSpaceGroup (false);
var i = (isBio ? symTemp.addBioMoleculeOperation (sym.spaceGroup.finalOperations[iop], op < 0) : symTemp.addSpaceGroupOperation ((op < 0 ? "!" : "=") + xyz, Math.abs (op)));
if (i < 0) return "";
var opTemp = symTemp.getSpaceGroupOperation (i);
if (!isBio) opTemp.centering = centering;
var info;
if (pt != null || iAtom >= 0) pt = JU.P3.newP (pt == null ? modelSet.at[iAtom] : pt);
if (type == 135266320) {
if (isBio) return "";
symTemp.setUnitCell (uc.getNotionalUnitCell (), false);
uc.toFractional (pt, false);
if (Float.isNaN (pt.x)) return "";
var sympt =  new JU.P3 ();
symTemp.newSpaceGroupPoint (i, pt, sympt, 0, 0, 0);
symTemp.toCartesian (sympt, false);
return sympt;
}info = opTemp.getDescription (modelSet, uc, pt, pt2, (id == null ? "sym" : id));
var ang = (info[9]).intValue ();
switch (type) {
case 135266306:
return info;
case 1073742001:
var sinfo = [info[0], info[1], info[2], JU.Escape.eP (info[4]), JU.Escape.eP (info[5]), JU.Escape.eP (info[6]), JU.Escape.eP (info[7]), JU.Escape.eP (info[8]), "" + info[9], "" + JU.Escape.e (info[10])];
return sinfo;
case 1073741982:
return info[0];
default:
case 1826248716:
return info[2];
case 135176:
return info[3];
case 1073742178:
return info[5];
case 12289:
return info[6];
case 135266320:
return info[7];
case 1073741854:
case 135266319:
return ((ang == 0) == (type == 135266319) ? info[8] : null);
case 135266305:
return info[9];
case 12:
return info[10];
}
}, "JM.ModelSet,~N,~N,J.api.SymmetryInterface,~S,~N,JU.P3,JU.P3,~S,~N");
Clazz_overrideMethod (c$, "fcoord", 
function (p) {
return JS.SymmetryOperation.fcoord (p);
}, "JU.T3");
Clazz_overrideMethod (c$, "getV0abc", 
function (def) {
if (this.unitCell == null) return null;
var m;
var isRev = false;
if (Clazz_instanceOf (def, String)) {
var sdef = def;
if (sdef.indexOf (";") < 0) sdef += ";0,0,0";
isRev = sdef.startsWith ("!");
if (isRev) sdef = sdef.substring (1);
var symTemp =  new JS.Symmetry ();
symTemp.setSpaceGroup (false);
var i = symTemp.addSpaceGroupOperation ("=" + sdef, 0);
if (i < 0) return null;
m = symTemp.getSpaceGroupOperation (i);
(m).doFinalize ();
} else {
m = (Clazz_instanceOf (def, JU.M3) ? JU.M4.newMV (def,  new JU.P3 ()) : def);
}var pts =  new Array (4);
var pt =  new JU.P3 ();
var m3 =  new JU.M3 ();
m.getRotationScale (m3);
m.getTranslation (pt);
if (isRev) {
m3.invert ();
m3.transpose ();
m3.rotate (pt);
pt.scale (-1);
} else {
m3.transpose ();
}this.unitCell.toCartesian (pt, false);
pts[0] = JU.V3.newV (pt);
pts[1] = JU.V3.new3 (1, 0, 0);
pts[2] = JU.V3.new3 (0, 1, 0);
pts[3] = JU.V3.new3 (0, 0, 1);
for (var i = 1; i < 4; i++) {
m3.rotate (pts[i]);
this.unitCell.toCartesian (pts[i], true);
}
return pts;
}, "~O");
Clazz_overrideMethod (c$, "getQuaternionRotation", 
function (abc) {
return (this.unitCell == null ? null : this.unitCell.getQuaternionRotation (abc));
}, "~S");
});
Clazz_declarePackage ("JS");
Clazz_load (["JU.P3", "$.V3"], "JS.PointGroup", ["java.lang.Float", "java.util.Hashtable", "JU.Lst", "$.Quat", "$.SB", "JU.BSUtil", "$.Escape", "$.Logger", "$.Txt"], function () {
c$ = Clazz_decorateAsClass (function () {
this.drawInfo = null;
this.info = null;
this.textInfo = null;
this.drawType = "";
this.drawIndex = 0;
this.scale = NaN;
this.nAxes = null;
this.axes = null;
this.nAtoms = 0;
this.radius = 0;
this.distanceTolerance = 0.2;
this.linearTolerance = 8;
this.cosTolerance = 0.99;
this.name = "C_1?";
this.principalAxis = null;
this.principalPlane = null;
this.vTemp = null;
this.centerAtomIndex = -1;
this.haveInversionCenter = false;
this.center = null;
this.points = null;
this.atoms = null;
this.elements = null;
this.bsAtoms = null;
this.maxElement = 0;
this.eCounts = null;
this.nOps = 0;
if (!Clazz_isClassDefined ("JS.PointGroup.Operation")) {
JS.PointGroup.$PointGroup$Operation$ ();
}
Clazz_instantialize (this, arguments);
}, JS, "PointGroup");
Clazz_prepareFields (c$, function () {
this.nAxes =  Clazz_newIntArray (JS.PointGroup.maxAxis, 0);
this.axes =  new Array (JS.PointGroup.maxAxis);
this.vTemp =  new JU.V3 ();
this.center =  new JU.P3 ();
});
Clazz_defineMethod (c$, "getName", 
function () {
return this.name;
});
c$.getPointGroup = Clazz_defineMethod (c$, "getPointGroup", 
function (pgLast, atomset, bsAtoms, haveVibration, distanceTolerance, linearTolerance) {
var pg =  new JS.PointGroup ();
return (pg.set (pgLast, atomset, bsAtoms, haveVibration, distanceTolerance, linearTolerance) ? pg : pgLast);
}, "JS.PointGroup,~A,JU.BS,~B,~N,~N");
Clazz_makeConstructor (c$, 
 function () {
});
Clazz_defineMethod (c$, "isEqual", 
 function (pg) {
if (pg == null) return false;
if (this.linearTolerance != pg.linearTolerance || this.distanceTolerance != pg.distanceTolerance || this.nAtoms != pg.nAtoms || !this.bsAtoms.equals (pg.bsAtoms)) return false;
for (var i = 0; i < this.nAtoms; i++) {
if (this.elements[i] != pg.elements[i] || this.points[i].distance (pg.points[i]) != 0) return false;
}
return true;
}, "JS.PointGroup");
Clazz_defineMethod (c$, "set", 
 function (pgLast, atomset, bsAtoms, haveVibration, distanceTolerance, linearTolerance) {
this.distanceTolerance = distanceTolerance;
this.linearTolerance = linearTolerance;
this.bsAtoms = bsAtoms;
this.cosTolerance = (Math.cos (linearTolerance / 180 * 3.141592653589793));
if (!this.getAtomsAndElements (atomset, bsAtoms)) {
JU.Logger.error ("Too many atoms for point group calculation");
this.name = "point group not determined -- ac > 100 -- select fewer atoms and try again.";
return true;
}this.getElementCounts ();
if (haveVibration) {
var atomVibs =  new Array (this.points.length);
for (var i = this.points.length; --i >= 0; ) {
atomVibs[i] = JU.P3.newP (this.points[i]);
var v = this.atoms[i].getVibrationVector ();
if (v != null) atomVibs[i].add (v);
}
this.points = atomVibs;
}if (this.isEqual (pgLast)) return false;
this.findInversionCenter ();
if (this.isLinear (this.points)) {
if (this.haveInversionCenter) {
this.name = "D(infinity)h";
} else {
this.name = "C(infinity)v";
}this.vTemp.sub2 (this.points[1], this.points[0]);
this.addAxis (16, this.vTemp);
this.principalAxis = this.axes[16][0];
if (this.haveInversionCenter) {
this.axes[0] =  new Array (1);
this.principalPlane = this.axes[0][this.nAxes[0]++] = Clazz_innerTypeInstance (JS.PointGroup.Operation, this, null, this.vTemp);
}return true;
}this.axes[0] =  new Array (15);
var nPlanes = 0;
this.findCAxes ();
nPlanes = this.findPlanes ();
this.findAdditionalAxes (nPlanes);
var n = this.getHighestOrder ();
if (this.nAxes[17] > 1) {
if (this.nAxes[19] > 1) {
if (this.haveInversionCenter) {
this.name = "Ih";
} else {
this.name = "I";
}} else if (this.nAxes[18] > 1) {
if (this.haveInversionCenter) {
this.name = "Oh";
} else {
this.name = "O";
}} else {
if (nPlanes > 0) {
if (this.haveInversionCenter) {
this.name = "Th";
} else {
this.name = "Td";
}} else {
this.name = "T";
}}} else {
if (n < 2) {
if (nPlanes == 1) {
this.name = "Cs";
return true;
}if (this.haveInversionCenter) {
this.name = "Ci";
return true;
}this.name = "C1";
} else if ((n % 2) == 1 && this.nAxes[16] > 0 || (n % 2) == 0 && this.nAxes[16] > 1) {
this.principalAxis = this.setPrincipalAxis (n, nPlanes);
if (nPlanes == 0) {
if (n < 14) {
this.name = "S" + n;
} else {
this.name = "D" + (n - 14);
}} else {
if (n < 14) n = Clazz_doubleToInt (n / 2);
 else n -= 14;
if (nPlanes == n) {
this.name = "D" + n + "d";
} else {
this.name = "D" + n + "h";
}}} else if (nPlanes == 0) {
this.principalAxis = this.axes[n][0];
if (n < 14) {
this.name = "S" + n;
} else {
this.name = "C" + (n - 14);
}} else if (nPlanes == n - 14) {
this.principalAxis = this.axes[n][0];
this.name = "C" + nPlanes + "v";
} else {
this.principalAxis = this.axes[n < 14 ? n + 14 : n][0];
this.principalPlane = this.axes[0][0];
if (n < 14) n /= 2;
 else n -= 14;
this.name = "C" + n + "h";
}}return true;
}, "JS.PointGroup,~A,JU.BS,~B,~N,~N");
Clazz_defineMethod (c$, "setPrincipalAxis", 
 function (n, nPlanes) {
var principalPlane = this.setPrincipalPlane (n, nPlanes);
if (nPlanes == 0 && n < 14 || this.nAxes[n] == 1) {
if (nPlanes > 0 && n < 14) n = 14 + Clazz_doubleToInt (n / 2);
return this.axes[n][0];
}if (principalPlane == null) return null;
for (var i = 0; i < this.nAxes[16]; i++) if (this.isParallel (principalPlane.normalOrAxis, this.axes[16][i].normalOrAxis)) {
if (i != 0) {
var o = this.axes[16][0];
this.axes[16][0] = this.axes[16][i];
this.axes[16][i] = o;
}return this.axes[16][0];
}
return null;
}, "~N,~N");
Clazz_defineMethod (c$, "setPrincipalPlane", 
 function (n, nPlanes) {
if (nPlanes == 1) return this.principalPlane = this.axes[0][0];
if (nPlanes == 0 || nPlanes == n - 14) return null;
for (var i = 0; i < nPlanes; i++) for (var j = 0, nPerp = 0; j < nPlanes; j++) if (this.isPerpendicular (this.axes[0][i].normalOrAxis, this.axes[0][j].normalOrAxis) && ++nPerp > 2) {
if (i != 0) {
var o = this.axes[0][0];
this.axes[0][0] = this.axes[0][i];
this.axes[0][i] = o;
}return this.principalPlane = this.axes[0][0];
}

return null;
}, "~N,~N");
Clazz_defineMethod (c$, "getAtomsAndElements", 
 function (atomset, bsAtoms) {
var ac = JU.BSUtil.cardinalityOf (bsAtoms);
if (ac > 100) return false;
this.points =  new Array (ac);
this.atoms =  new Array (ac);
this.elements =  Clazz_newIntArray (ac, 0);
if (ac == 0) return true;
this.nAtoms = 0;
for (var i = bsAtoms.nextSetBit (0); i >= 0; i = bsAtoms.nextSetBit (i + 1)) {
this.points[this.nAtoms] = JU.P3.newP (atomset[i]);
this.atoms[this.nAtoms] = atomset[i];
var bondIndex = 1 + Math.max (3, atomset[i].getCovalentBondCount ());
this.elements[this.nAtoms] = atomset[i].getElementNumber () * bondIndex;
this.center.add (this.points[this.nAtoms++]);
}
this.center.scale (1 / this.nAtoms);
for (var i = this.nAtoms; --i >= 0; ) {
var r = this.center.distance (this.points[i]);
if (r < this.distanceTolerance) this.centerAtomIndex = i;
this.radius = Math.max (this.radius, r);
}
return true;
}, "~A,JU.BS");
Clazz_defineMethod (c$, "findInversionCenter", 
 function () {
this.haveInversionCenter = this.checkOperation (null, this.center, -1);
if (this.haveInversionCenter) {
this.axes[1] =  new Array (1);
this.axes[1][0] = Clazz_innerTypeInstance (JS.PointGroup.Operation, this, null);
}});
Clazz_defineMethod (c$, "checkOperation", 
 function (q, center, iOrder) {
var pt =  new JU.P3 ();
var nFound = 0;
var isInversion = (iOrder < 14);
out : for (var i = this.points.length; --i >= 0 && nFound < this.points.length; ) if (i == this.centerAtomIndex) {
nFound++;
} else {
var a1 = this.points[i];
var e1 = this.elements[i];
if (q != null) {
pt.sub2 (a1, center);
q.transformP2 (pt, pt).add (center);
} else {
pt.setT (a1);
}if (isInversion) {
this.vTemp.sub2 (center, pt);
pt.scaleAdd2 (2, this.vTemp, pt);
}if ((q != null || isInversion) && pt.distance (a1) < this.distanceTolerance) {
nFound++;
continue;
}for (var j = this.points.length; --j >= 0; ) {
if (j == i || this.elements[j] != e1) continue;
var a2 = this.points[j];
if (pt.distance (a2) < this.distanceTolerance) {
nFound++;
continue out;
}}
}
return nFound == this.points.length;
}, "JU.Quat,JU.P3,~N");
Clazz_defineMethod (c$, "isLinear", 
 function (atoms) {
var v1 = null;
if (atoms.length < 2) return false;
for (var i = atoms.length; --i >= 0; ) {
if (i == this.centerAtomIndex) continue;
if (v1 == null) {
v1 =  new JU.V3 ();
v1.sub2 (atoms[i], this.center);
v1.normalize ();
this.vTemp.setT (v1);
continue;
}this.vTemp.sub2 (atoms[i], this.center);
this.vTemp.normalize ();
if (!this.isParallel (v1, this.vTemp)) return false;
}
return true;
}, "~A");
Clazz_defineMethod (c$, "isParallel", 
 function (v1, v2) {
return (Math.abs (v1.dot (v2)) >= this.cosTolerance);
}, "JU.V3,JU.V3");
Clazz_defineMethod (c$, "isPerpendicular", 
 function (v1, v2) {
return (Math.abs (v1.dot (v2)) <= 1 - this.cosTolerance);
}, "JU.V3,JU.V3");
Clazz_defineMethod (c$, "getElementCounts", 
 function () {
for (var i = this.points.length; --i >= 0; ) {
var e1 = this.elements[i];
if (e1 > this.maxElement) this.maxElement = e1;
}
this.eCounts =  Clazz_newIntArray (++this.maxElement, 0);
for (var i = this.points.length; --i >= 0; ) this.eCounts[this.elements[i]]++;

});
Clazz_defineMethod (c$, "findCAxes", 
 function () {
var v1 =  new JU.V3 ();
var v2 =  new JU.V3 ();
var v3 =  new JU.V3 ();
for (var i = this.points.length; --i >= 0; ) {
if (i == this.centerAtomIndex) continue;
var a1 = this.points[i];
var e1 = this.elements[i];
for (var j = this.points.length; --j > i; ) {
var a2 = this.points[j];
if (this.elements[j] != e1) continue;
v1.sub2 (a1, this.center);
v2.sub2 (a2, this.center);
v1.normalize ();
v2.normalize ();
if (this.isParallel (v1, v2)) {
this.getAllAxes (v1);
continue;
}if (this.nAxes[16] < JS.PointGroup.axesMaxN[16]) {
v3.ave (a1, a2);
v3.sub (this.center);
this.getAllAxes (v3);
}var order = (6.283185307179586 / v1.angle (v2));
var iOrder = Clazz_doubleToInt (Math.floor (order + 0.01));
var isIntegerOrder = (order - iOrder <= 0.02);
if (!isIntegerOrder || (iOrder = iOrder + 14) >= JS.PointGroup.maxAxis) continue;
if (this.nAxes[iOrder] < JS.PointGroup.axesMaxN[iOrder]) {
v3.cross (v1, v2);
this.checkAxisOrder (iOrder, v3, this.center);
}}
}
var vs =  new Array (this.nAxes[16] * 2);
for (var i = 0; i < vs.length; i++) vs[i] =  new JU.V3 ();

var n = 0;
for (var i = 0; i < this.nAxes[16]; i++) {
vs[n++].setT (this.axes[16][i].normalOrAxis);
vs[n].setT (this.axes[16][i].normalOrAxis);
vs[n++].scale (-1);
}
for (var i = vs.length; --i >= 2; ) for (var j = i; --j >= 1; ) for (var k = j; --k >= 0; ) {
v3.add2 (vs[i], vs[j]);
v3.add (vs[k]);
if (v3.length () < 1.0) continue;
this.checkAxisOrder (17, v3, this.center);
}


var nMin = 2147483647;
var iMin = -1;
for (var i = 0; i < this.maxElement; i++) {
if (this.eCounts[i] < nMin && this.eCounts[i] > 2) {
nMin = this.eCounts[i];
iMin = i;
}}
out : for (var i = 0; i < this.points.length - 2; i++) if (this.elements[i] == iMin) for (var j = i + 1; j < this.points.length - 1; j++) if (this.elements[j] == iMin) for (var k = j + 1; k < this.points.length; k++) if (this.elements[k] == iMin) {
v1.sub2 (this.points[i], this.points[j]);
v2.sub2 (this.points[i], this.points[k]);
v1.normalize ();
v2.normalize ();
v3.cross (v1, v2);
this.getAllAxes (v3);
v1.add2 (this.points[i], this.points[j]);
v1.add (this.points[k]);
v1.normalize ();
if (!this.isParallel (v1, v3)) this.getAllAxes (v1);
if (this.nAxes[19] == JS.PointGroup.axesMaxN[19]) break out;
}


vs =  new Array (this.maxElement);
for (var i = this.points.length; --i >= 0; ) {
var e1 = this.elements[i];
if (vs[e1] == null) vs[e1] =  new JU.V3 ();
 else if (this.haveInversionCenter) continue;
vs[e1].add (this.points[i]);
}
if (!this.haveInversionCenter) for (var i = 0; i < this.maxElement; i++) if (vs[i] != null) vs[i].scale (1 / this.eCounts[i]);

for (var i = 0; i < this.maxElement; i++) if (vs[i] != null) for (var j = 0; j < this.maxElement; j++) {
if (i == j || vs[j] == null) continue;
if (this.haveInversionCenter) v1.cross (vs[i], vs[j]);
 else v1.sub2 (vs[i], vs[j]);
this.checkAxisOrder (16, v1, this.center);
}

return this.getHighestOrder ();
});
Clazz_defineMethod (c$, "getAllAxes", 
 function (v3) {
for (var o = 16; o < JS.PointGroup.maxAxis; o++) if (this.nAxes[o] < JS.PointGroup.axesMaxN[o]) this.checkAxisOrder (o, v3, this.center);

}, "JU.V3");
Clazz_defineMethod (c$, "getHighestOrder", 
 function () {
var n = 0;
for (n = 14; --n > 1 && this.nAxes[n] == 0; ) {
}
if (n > 1) return (n + 14 < JS.PointGroup.maxAxis && this.nAxes[n + 14] > 0 ? n + 14 : n);
for (n = JS.PointGroup.maxAxis; --n > 1 && this.nAxes[n] == 0; ) {
}
return n;
});
Clazz_defineMethod (c$, "checkAxisOrder", 
 function (iOrder, v, center) {
switch (iOrder) {
case 22:
if (this.nAxes[17] > 0) return false;
case 20:
case 18:
if (this.nAxes[19] > 0) return false;
break;
case 17:
if (this.nAxes[22] > 0) return false;
break;
case 19:
if (this.nAxes[18] > 0 || this.nAxes[20] > 0 || this.nAxes[22] > 0) return false;
break;
}
v.normalize ();
if (this.haveAxis (iOrder, v)) return false;
var q = JU.Quat.newVA (v, (iOrder < 14 ? 180 : 0) + Clazz_doubleToInt (360 / (iOrder % 14)));
if (!this.checkOperation (q, center, iOrder)) return false;
this.addAxis (iOrder, v);
switch (iOrder) {
case 16:
this.checkAxisOrder (4, v, center);
break;
case 17:
this.checkAxisOrder (3, v, center);
if (this.haveInversionCenter) this.addAxis (6, v);
break;
case 18:
this.addAxis (16, v);
this.checkAxisOrder (4, v, center);
this.checkAxisOrder (8, v, center);
break;
case 19:
this.checkAxisOrder (5, v, center);
if (this.haveInversionCenter) this.addAxis (10, v);
break;
case 20:
this.addAxis (16, v);
this.addAxis (17, v);
this.checkAxisOrder (3, v, center);
this.checkAxisOrder (6, v, center);
this.checkAxisOrder (12, v, center);
break;
case 22:
this.addAxis (16, v);
this.addAxis (18, v);
break;
}
return true;
}, "~N,JU.V3,JU.P3");
Clazz_defineMethod (c$, "addAxis", 
 function (iOrder, v) {
if (this.haveAxis (iOrder, v)) return;
if (this.axes[iOrder] == null) this.axes[iOrder] =  new Array (JS.PointGroup.axesMaxN[iOrder]);
this.axes[iOrder][this.nAxes[iOrder]++] = Clazz_innerTypeInstance (JS.PointGroup.Operation, this, null, v, iOrder);
}, "~N,JU.V3");
Clazz_defineMethod (c$, "haveAxis", 
 function (iOrder, v) {
if (this.nAxes[iOrder] == JS.PointGroup.axesMaxN[iOrder]) {
return true;
}if (this.nAxes[iOrder] > 0) for (var i = this.nAxes[iOrder]; --i >= 0; ) {
if (this.isParallel (v, this.axes[iOrder][i].normalOrAxis)) return true;
}
return false;
}, "~N,JU.V3");
Clazz_defineMethod (c$, "findPlanes", 
 function () {
var pt =  new JU.P3 ();
var v1 =  new JU.V3 ();
var v2 =  new JU.V3 ();
var v3 =  new JU.V3 ();
var nPlanes = 0;
var haveAxes = (this.getHighestOrder () > 1);
for (var i = this.points.length; --i >= 0; ) {
if (i == this.centerAtomIndex) continue;
var a1 = this.points[i];
var e1 = this.elements[i];
for (var j = this.points.length; --j > i; ) {
if (haveAxes && this.elements[j] != e1) continue;
var a2 = this.points[j];
pt.add2 (a1, a2);
pt.scale (0.5);
v1.sub2 (a1, this.center);
v2.sub2 (a2, this.center);
if (!this.isParallel (v1, v2)) {
v3.cross (v1, v2);
v3.normalize ();
nPlanes = this.getPlane (v3);
}v3.sub2 (a2, a1);
v3.normalize ();
nPlanes = this.getPlane (v3);
if (nPlanes == JS.PointGroup.axesMaxN[0]) return nPlanes;
}
}
if (haveAxes) for (var i = 16; i < JS.PointGroup.maxAxis; i++) for (var j = 0; j < this.nAxes[i]; j++) nPlanes = this.getPlane (this.axes[i][j].normalOrAxis);


return nPlanes;
});
Clazz_defineMethod (c$, "getPlane", 
 function (v3) {
if (!this.haveAxis (0, v3) && this.checkOperation (JU.Quat.newVA (v3, 180), this.center, -1)) this.axes[0][this.nAxes[0]++] = Clazz_innerTypeInstance (JS.PointGroup.Operation, this, null, v3);
return this.nAxes[0];
}, "JU.V3");
Clazz_defineMethod (c$, "findAdditionalAxes", 
 function (nPlanes) {
var planes = this.axes[0];
var Cn = 0;
if (nPlanes > 1 && ((Cn = nPlanes + 14) < JS.PointGroup.maxAxis) && this.nAxes[Cn] == 0) {
this.vTemp.cross (planes[0].normalOrAxis, planes[1].normalOrAxis);
if (!this.checkAxisOrder (Cn, this.vTemp, this.center) && nPlanes > 2) {
this.vTemp.cross (planes[1].normalOrAxis, planes[2].normalOrAxis);
this.checkAxisOrder (Cn - 1, this.vTemp, this.center);
}}if (this.nAxes[16] == 0 && nPlanes > 2) {
for (var i = 0; i < nPlanes - 1; i++) {
for (var j = i + 1; j < nPlanes; j++) {
this.vTemp.add2 (planes[1].normalOrAxis, planes[2].normalOrAxis);
this.checkAxisOrder (16, this.vTemp, this.center);
}
}
}}, "~N");
Clazz_defineMethod (c$, "getInfo", 
function (modelIndex, asDraw, asInfo, type, index, scaleFactor) {
this.info = (asInfo ?  new java.util.Hashtable () : null);
var v =  new JU.V3 ();
var op;
if (scaleFactor == 0) scaleFactor = 1;
this.scale = scaleFactor;
var nType =  Clazz_newIntArray (4, 2, 0);
for (var i = 1; i < JS.PointGroup.maxAxis; i++) for (var j = this.nAxes[i]; --j >= 0; ) nType[this.axes[i][j].type][0]++;


var sb =  new JU.SB ().append ("# ").appendI (this.nAtoms).append (" atoms\n");
if (asDraw) {
var haveType = (type != null && type.length > 0);
this.drawType = type = (haveType ? type : "");
this.drawIndex = index;
var anyProperAxis = (type.equalsIgnoreCase ("Cn"));
var anyImproperAxis = (type.equalsIgnoreCase ("Sn"));
sb.append ("set perspectivedepth off;\n");
var m = "_" + modelIndex + "_";
if (!haveType) sb.append ("draw pg0").append (m).append ("* delete;draw pgva").append (m).append ("* delete;draw pgvp").append (m).append ("* delete;");
if (!haveType || type.equalsIgnoreCase ("Ci")) sb.append ("draw pg0").append (m).append (this.haveInversionCenter ? "inv " : " ").append (JU.Escape.eP (this.center)).append (this.haveInversionCenter ? "\"i\";\n" : ";\n");
var offset = 0.1;
for (var i = 2; i < JS.PointGroup.maxAxis; i++) {
if (i == 14) offset = 0.1;
if (this.nAxes[i] == 0) continue;
var label = this.axes[i][0].getLabel ();
offset += 0.25;
var scale = scaleFactor * this.radius + offset;
if (!haveType || type.equalsIgnoreCase (label) || anyProperAxis && i >= 14 || anyImproperAxis && i < 14) for (var j = 0; j < this.nAxes[i]; j++) {
if (index > 0 && j + 1 != index) continue;
op = this.axes[i][j];
v.add2 (op.normalOrAxis, this.center);
if (op.type == 2) scale = -scale;
sb.append ("draw pgva").append (m).append (label).append ("_").appendI (j + 1).append (" width 0.05 scale ").appendF (scale).append (" ").append (JU.Escape.eP (v));
v.scaleAdd2 (-2, op.normalOrAxis, v);
var isPA = (this.principalAxis != null && op.index == this.principalAxis.index);
sb.append (JU.Escape.eP (v)).append ("\"").append (label).append (isPA ? "*" : "").append ("\" color ").append (isPA ? "red" : op.type == 2 ? "blue" : "yellow").append (";\n");
}
}
if (!haveType || type.equalsIgnoreCase ("Cs")) for (var j = 0; j < this.nAxes[0]; j++) {
if (index > 0 && j + 1 != index) continue;
op = this.axes[0][j];
sb.append ("draw pgvp").append (m).appendI (j + 1).append ("disk scale ").appendF (scaleFactor * this.radius * 2).append (" CIRCLE PLANE ").append (JU.Escape.eP (this.center));
v.add2 (op.normalOrAxis, this.center);
sb.append (JU.Escape.eP (v)).append (" color translucent yellow;\n");
v.add2 (op.normalOrAxis, this.center);
sb.append ("draw pgvp").append (m).appendI (j + 1).append ("ring width 0.05 scale ").appendF (scaleFactor * this.radius * 2).append (" arc ").append (JU.Escape.eP (v));
v.scaleAdd2 (-2, op.normalOrAxis, v);
sb.append (JU.Escape.eP (v));
v.add3 (0.011, 0.012, 0.013);
sb.append (JU.Escape.eP (v)).append ("{0 360 0.5} color ").append (this.principalPlane != null && op.index == this.principalPlane.index ? "red" : "blue").append (";\n");
}
sb.append ("# name=").append (this.name);
sb.append (", nCi=").appendI (this.haveInversionCenter ? 1 : 0);
sb.append (", nCs=").appendI (this.nAxes[0]);
sb.append (", nCn=").appendI (nType[1][0]);
sb.append (", nSn=").appendI (nType[2][0]);
sb.append (": ");
for (var i = JS.PointGroup.maxAxis; --i >= 2; ) if (this.nAxes[i] > 0) {
sb.append (" n").append (i < 14 ? "S" : "C").appendI (i % 14);
sb.append ("=").appendI (this.nAxes[i]);
}
sb.append (";\n");
this.drawInfo = sb.toString ();
return this.drawInfo;
}var n = 0;
var nTotal = 1;
var ctype = (this.haveInversionCenter ? "Ci" : "center");
if (this.haveInversionCenter) nTotal++;
if (this.info == null) sb.append ("\n\n").append (this.name).append ("\t").append (ctype).append ("\t").append (JU.Escape.eP (this.center));
 else this.info.put (ctype, this.center);
for (var i = JS.PointGroup.maxAxis; --i >= 0; ) {
if (this.nAxes[i] > 0) {
n = JS.PointGroup.nUnique[i];
var label = this.axes[i][0].getLabel ();
if (this.info == null) sb.append ("\n\n").append (this.name).append ("\tn").append (label).append ("\t").appendI (this.nAxes[i]).append ("\t").appendI (n);
 else this.info.put ("n" + label, Integer.$valueOf (this.nAxes[i]));
n *= this.nAxes[i];
nTotal += n;
nType[this.axes[i][0].type][1] += n;
var vinfo = (this.info == null ? null :  new JU.Lst ());
for (var j = 0; j < this.nAxes[i]; j++) {
if (vinfo == null) sb.append ("\n").append (this.name).append ("\t").append (label).append ("_").appendI (j + 1).append ("\t").appendO (this.axes[i][j].normalOrAxis);
 else vinfo.addLast (this.axes[i][j].normalOrAxis);
}
if (this.info != null) this.info.put (label, vinfo);
}}
if (this.info == null) {
sb.append ("\n");
sb.append ("\n").append (this.name).append ("\ttype\tnType\tnUnique");
sb.append ("\n").append (this.name).append ("\tE\t  1\t  1");
n = (this.haveInversionCenter ? 1 : 0);
sb.append ("\n").append (this.name).append ("\tCi\t  ").appendI (n).append ("\t  ").appendI (n);
sb.append ("\n").append (this.name).append ("\tCs\t");
JU.Txt.rightJustify (sb, "    ", this.nAxes[0] + "\t");
JU.Txt.rightJustify (sb, "    ", this.nAxes[0] + "\n");
sb.append (this.name).append ("\tCn\t");
JU.Txt.rightJustify (sb, "    ", nType[1][0] + "\t");
JU.Txt.rightJustify (sb, "    ", nType[1][1] + "\n");
sb.append (this.name).append ("\tSn\t");
JU.Txt.rightJustify (sb, "    ", nType[2][0] + "\t");
JU.Txt.rightJustify (sb, "    ", nType[2][1] + "\n");
sb.append (this.name).append ("\t\tTOTAL\t");
JU.Txt.rightJustify (sb, "    ", nTotal + "\n");
this.textInfo = sb.toString ();
return this.textInfo;
}this.info.put ("name", this.name);
this.info.put ("nAtoms", Integer.$valueOf (this.nAtoms));
this.info.put ("nTotal", Integer.$valueOf (nTotal));
this.info.put ("nCi", Integer.$valueOf (this.haveInversionCenter ? 1 : 0));
this.info.put ("nCs", Integer.$valueOf (this.nAxes[0]));
this.info.put ("nCn", Integer.$valueOf (nType[1][0]));
this.info.put ("nSn", Integer.$valueOf (nType[2][0]));
this.info.put ("distanceTolerance", Float.$valueOf (this.distanceTolerance));
this.info.put ("linearTolerance", Float.$valueOf (this.linearTolerance));
this.info.put ("detail", sb.toString ().$replace ('\n', ';'));
if (this.principalAxis != null && this.principalAxis.index > 0) this.info.put ("principalAxis", this.principalAxis.normalOrAxis);
if (this.principalPlane != null && this.principalPlane.index > 0) this.info.put ("principalPlane", this.principalPlane.normalOrAxis);
return this.info;
}, "~N,~B,~B,~S,~N,~N");
Clazz_defineMethod (c$, "isDrawType", 
function (type, index, scale) {
return (this.drawInfo != null && this.drawType.equals (type == null ? "" : type) && this.drawIndex == index && this.scale == scale);
}, "~S,~N,~N");
c$.$PointGroup$Operation$ = function () {
Clazz_pu$h(self.c$);
c$ = Clazz_decorateAsClass (function () {
Clazz_prepareCallback (this, arguments);
this.type = 0;
this.order = 0;
this.index = 0;
this.normalOrAxis = null;
Clazz_instantialize (this, arguments);
}, JS.PointGroup, "Operation");
Clazz_makeConstructor (c$, 
function () {
this.index = ++this.b$["JS.PointGroup"].nOps;
this.type = 3;
this.order = 1;
if (JU.Logger.debugging) JU.Logger.debug ("new operation -- " + JS.PointGroup.typeNames[this.type]);
});
Clazz_makeConstructor (c$, 
function (a, b) {
this.index = ++this.b$["JS.PointGroup"].nOps;
this.type = (b < 14 ? 2 : 1);
this.order = b % 14;
this.normalOrAxis = JU.Quat.newVA (a, 180).getNormal ();
if (JU.Logger.debugging) JU.Logger.debug ("new operation -- " + (this.order == b ? "S" : "C") + this.order + " " + this.normalOrAxis);
}, "JU.V3,~N");
Clazz_makeConstructor (c$, 
function (a) {
if (a == null) return;
this.index = ++this.b$["JS.PointGroup"].nOps;
this.type = 0;
this.normalOrAxis = JU.Quat.newVA (a, 180).getNormal ();
if (JU.Logger.debugging) JU.Logger.debug ("new operation -- plane " + this.normalOrAxis);
}, "JU.V3");
Clazz_defineMethod (c$, "getLabel", 
function () {
switch (this.type) {
case 0:
return "Cs";
case 2:
return "S" + this.order;
default:
return "C" + this.order;
}
});
c$ = Clazz_p0p ();
};
Clazz_defineStatics (c$,
"axesMaxN", [15, 0, 0, 1, 3, 1, 10, 0, 1, 0, 6, 0, 1, 0, 0, 0, 15, 10, 6, 6, 10, 0, 1],
"nUnique", [1, 0, 0, 2, 2, 4, 2, 0, 4, 0, 4, 0, 4, 0, 0, 0, 1, 2, 2, 4, 2, 0, 4],
"s3", 3,
"s4", 4,
"s5", 5,
"s6", 6,
"s8", 8,
"s10", 10,
"s12", 12,
"firstProper", 14,
"c2", 16,
"c3", 17,
"c4", 18,
"c5", 19,
"c6", 20,
"c8", 22);
c$.maxAxis = c$.prototype.maxAxis = JS.PointGroup.axesMaxN.length;
Clazz_defineStatics (c$,
"ATOM_COUNT_MAX", 100,
"OPERATION_PLANE", 0,
"OPERATION_PROPER_AXIS", 1,
"OPERATION_IMPROPER_AXIS", 2,
"OPERATION_INVERSION_CENTER", 3,
"typeNames", ["plane", "proper axis", "improper axis", "center of inversion"]);
});
Clazz_declarePackage ("JS");
Clazz_load (["java.util.Hashtable"], "JS.SpaceGroup", ["java.lang.Character", "$.Float", "java.util.Arrays", "JU.AU", "$.Lst", "$.M4", "$.P3", "$.PT", "$.SB", "JS.HallInfo", "$.HallTranslation", "$.SymmetryOperation", "JU.Logger"], function () {
c$ = Clazz_decorateAsClass (function () {
this.index = 0;
this.isSSG = false;
this.name = "unknown!";
this.hallSymbol = null;
this.hmSymbol = null;
this.hmSymbolFull = null;
this.hmSymbolExt = null;
this.hmSymbolAbbr = null;
this.hmSymbolAlternative = null;
this.hmSymbolAbbrShort = null;
this.ambiguityType = '\0';
this.uniqueAxis = '\0';
this.axisChoice = '\0';
this.intlTableNumber = null;
this.intlTableNumberFull = null;
this.intlTableNumberExt = null;
this.hallInfo = null;
this.latticeParameter = 0;
this.latticeCode = '\0';
this.operations = null;
this.finalOperations = null;
this.operationCount = 0;
this.latticeOp = -1;
this.xyzList = null;
this.modDim = 0;
this.doNormalize = true;
this.isBio = false;
this.isBilbao = false;
this.latticeOps = null;
Clazz_instantialize (this, arguments);
}, JS, "SpaceGroup");
Clazz_prepareFields (c$, function () {
this.xyzList =  new java.util.Hashtable ();
});
c$.getNull = Clazz_defineMethod (c$, "getNull", 
function (doInit) {
JS.SpaceGroup.getSpaceGroups ();
return  new JS.SpaceGroup (null, doInit);
}, "~B");
Clazz_makeConstructor (c$, 
 function (cifLine, doInit) {
this.index = ++JS.SpaceGroup.sgIndex;
if (!doInit) return;
if (cifLine == null) {
this.addSymmetry ("x,y,z", 0, false);
} else {
this.buildSpaceGroup (cifLine);
}}, "~S,~B");
Clazz_defineMethod (c$, "set", 
function (doNormalize) {
this.doNormalize = doNormalize;
return this;
}, "~B");
Clazz_defineMethod (c$, "init", 
 function () {
this.xyzList =  new java.util.Hashtable ();
this.operationCount = 0;
this.addSymmetry ("x,y,z", 0, false);
});
c$.createSpaceGroup = Clazz_defineMethod (c$, "createSpaceGroup", 
function (desiredSpaceGroupIndex, name, data) {
var sg = null;
if (desiredSpaceGroupIndex >= 0) {
sg = JS.SpaceGroup.getSpaceGroups ()[desiredSpaceGroupIndex];
} else {
if (Clazz_instanceOf (data, JU.Lst)) sg = JS.SpaceGroup.createSGFromList (name, data);
 else sg = JS.SpaceGroup.determineSpaceGroupNA (name, data);
if (sg == null) sg = JS.SpaceGroup.createSpaceGroupN (name);
}if (sg != null) sg.generateAllOperators (null);
return sg;
}, "~N,~S,~O");
c$.createSGFromList = Clazz_defineMethod (c$, "createSGFromList", 
 function (name, data) {
var sg =  new JS.SpaceGroup ("0;--;--;--", true);
sg.doNormalize = false;
sg.name = name;
var n = data.size ();
for (var i = 0; i < n; i++) {
var operation = data.get (i);
if (Clazz_instanceOf (operation, JS.SymmetryOperation)) {
var op = operation;
var iop = sg.addOp (op, op.xyz, false);
sg.operations[iop].timeReversal = op.timeReversal;
} else {
sg.addSymmetrySM ("xyz matrix:" + operation, operation);
}}
var sgn = sg.getDerivedSpaceGroup ();
if (sgn != null) sg = sgn;
return sg;
}, "~S,JU.Lst");
Clazz_defineMethod (c$, "addSymmetry", 
function (xyz, opId, allowScaling) {
xyz = xyz.toLowerCase ();
return (xyz.indexOf ("[[") < 0 && xyz.indexOf ("x4") < 0 && xyz.indexOf (";") < 0 && (xyz.indexOf ("x") < 0 || xyz.indexOf ("y") < 0 || xyz.indexOf ("z") < 0) ? -1 : this.addOperation (xyz, opId, allowScaling));
}, "~S,~N,~B");
Clazz_defineMethod (c$, "setFinalOperations", 
function (atoms, atomIndex, count, doNormalize) {
if (this.hallInfo == null && this.latticeParameter != 0) {
var h =  new JS.HallInfo (JS.HallTranslation.getHallLatticeEquivalent (this.latticeParameter));
this.generateAllOperators (h);
}this.finalOperations = null;
this.isBio = (this.name.indexOf ("bio") >= 0);
if (this.index >= JS.SpaceGroup.getSpaceGroups ().length && !this.isBio && this.name.indexOf ("SSG:") < 0 && this.name.indexOf ("[subsystem") < 0) {
var sg = this.getDerivedSpaceGroup ();
if (sg != null) this.name = sg.getName ();
}this.finalOperations =  new Array (this.operationCount);
if (doNormalize && count > 0 && atoms != null) {
this.finalOperations[0] =  new JS.SymmetryOperation (this.operations[0], atoms, atomIndex, count, true);
var atom = atoms[atomIndex];
var c = JU.P3.newP (atom);
this.finalOperations[0].rotTrans (c);
if (c.distance (atom) > 0.0001) for (var i = 0; i < count; i++) {
atom = atoms[atomIndex + i];
c.setT (atom);
this.finalOperations[0].rotTrans (c);
atom.setT (c);
}
}var centering = null;
for (var i = 0; i < this.operationCount; i++) {
this.finalOperations[i] =  new JS.SymmetryOperation (this.operations[i], atoms, atomIndex, count, doNormalize);
centering = this.finalOperations[i].setCentering (centering, true);
}
}, "~A,~N,~N,~B");
Clazz_defineMethod (c$, "getOperationCount", 
function () {
return this.finalOperations.length;
});
Clazz_defineMethod (c$, "getOperation", 
function (i) {
return this.finalOperations[i];
}, "~N");
Clazz_defineMethod (c$, "getXyz", 
function (i, doNormalize) {
return (this.finalOperations == null ? this.operations[i].getXyz (doNormalize) : this.finalOperations[i].getXyz (doNormalize));
}, "~N,~B");
Clazz_defineMethod (c$, "newPoint", 
function (i, atom1, atom2, transX, transY, transZ) {
this.finalOperations[i].newPoint (atom1, atom2, transX, transY, transZ);
}, "~N,JU.P3,JU.P3,~N,~N,~N");
c$.getInfo = Clazz_defineMethod (c$, "getInfo", 
function (spaceGroup, cellInfo) {
var sg;
if (cellInfo != null) {
if (spaceGroup.indexOf ("[") >= 0) spaceGroup = spaceGroup.substring (0, spaceGroup.indexOf ("[")).trim ();
if (spaceGroup.equals ("unspecified!")) return "no space group identified in file";
sg = JS.SpaceGroup.determineSpaceGroupNA (spaceGroup, cellInfo.getNotionalUnitCell ());
} else if (spaceGroup.equalsIgnoreCase ("ALL")) {
return JS.SpaceGroup.dumpAll ();
} else if (spaceGroup.equalsIgnoreCase ("ALLSEITZ")) {
return JS.SpaceGroup.dumpAllSeitz ();
} else {
sg = JS.SpaceGroup.determineSpaceGroupN (spaceGroup);
if (sg == null) {
sg = JS.SpaceGroup.createSpaceGroupN (spaceGroup);
} else {
var sb =  new JU.SB ();
while (sg != null) {
sb.append (sg.dumpInfo (null));
sg = JS.SpaceGroup.determineSpaceGroupNS (spaceGroup, sg);
}
return sb.toString ();
}}return sg == null ? "?" : sg.dumpInfo (cellInfo);
}, "~S,J.api.SymmetryInterface");
Clazz_defineMethod (c$, "dumpInfo", 
function (cellInfo) {
var info = this.dumpCanonicalSeitzList ();
if (Clazz_instanceOf (info, JS.SpaceGroup)) return (info).dumpInfo (null);
var sb =  new JU.SB ().append ("\nHermann-Mauguin symbol: ");
sb.append (this.hmSymbol).append (this.hmSymbolExt.length > 0 ? ":" + this.hmSymbolExt : "").append ("\ninternational table number: ").append (this.intlTableNumber).append (this.intlTableNumberExt.length > 0 ? ":" + this.intlTableNumberExt : "").append ("\n\n").appendI (this.operationCount).append (" operators").append (!this.hallInfo.hallSymbol.equals ("--") ? " from Hall symbol " + this.hallInfo.hallSymbol : "").append (": ");
for (var i = 0; i < this.operationCount; i++) {
sb.append ("\n").append (this.operations[i].xyz);
}
sb.append ("\n\n").append (this.hallInfo == null ? "invalid Hall symbol" : this.hallInfo.dumpInfo ());
sb.append ("\n\ncanonical Seitz: ").append (info).append ("\n----------------------------------------------------\n");
return sb.toString ();
}, "J.api.SymmetryInterface");
Clazz_defineMethod (c$, "getName", 
function () {
return this.name;
});
Clazz_defineMethod (c$, "getLatticeDesignation", 
function () {
return this.latticeCode + ": " + JS.HallTranslation.getLatticeDesignation (this.latticeParameter);
});
Clazz_defineMethod (c$, "setLatticeParam", 
function (latticeParameter) {
this.latticeParameter = latticeParameter;
this.latticeCode = JS.HallTranslation.getLatticeCode (latticeParameter);
if (latticeParameter > 10) {
this.latticeParameter = -JS.HallTranslation.getLatticeIndex (this.latticeCode);
}}, "~N");
Clazz_defineMethod (c$, "dumpCanonicalSeitzList", 
 function () {
if (this.hallInfo == null) this.hallInfo =  new JS.HallInfo (this.hallSymbol);
this.generateAllOperators (null);
var s = this.getCanonicalSeitzList ();
if (this.index >= JS.SpaceGroup.SG.length) {
var sgDerived = JS.SpaceGroup.findSpaceGroup (s);
if (sgDerived != null) return sgDerived;
}return (this.index >= 0 && this.index < JS.SpaceGroup.SG.length ? this.hallSymbol + " = " : "") + s;
});
Clazz_defineMethod (c$, "getDerivedSpaceGroup", 
function () {
if (this.index >= 0 && this.index < JS.SpaceGroup.SG.length || this.modDim > 0 || this.operations[0].timeReversal != 0) return this;
if (this.finalOperations != null) this.setFinalOperations (null, 0, 0, false);
var s = this.getCanonicalSeitzList ();
return (s == null ? null : JS.SpaceGroup.findSpaceGroup (s));
});
Clazz_defineMethod (c$, "getCanonicalSeitzList", 
 function () {
var list =  new Array (this.operationCount);
for (var i = 0; i < this.operationCount; i++) list[i] = JS.SymmetryOperation.dumpSeitz (this.operations[i], true);

java.util.Arrays.sort (list, 0, this.operationCount);
var sb =  new JU.SB ().append ("\n[");
for (var i = 0; i < this.operationCount; i++) sb.append (list[i].$replace ('\t', ' ').$replace ('\n', ' ')).append ("; ");

sb.append ("]");
return sb.toString ();
});
c$.findSpaceGroup = Clazz_defineMethod (c$, "findSpaceGroup", 
 function (s) {
JS.SpaceGroup.getSpaceGroups ();
if (JS.SpaceGroup.canonicalSeitzList == null) JS.SpaceGroup.canonicalSeitzList =  new Array (JS.SpaceGroup.SG.length);
for (var i = 0; i < JS.SpaceGroup.SG.length; i++) {
if (JS.SpaceGroup.canonicalSeitzList[i] == null) JS.SpaceGroup.canonicalSeitzList[i] = JS.SpaceGroup.SG[i].dumpCanonicalSeitzList ();
if (JS.SpaceGroup.canonicalSeitzList[i].indexOf (s) >= 0) return JS.SpaceGroup.SG[i];
}
return null;
}, "~S");
c$.dumpAll = Clazz_defineMethod (c$, "dumpAll", 
 function () {
var sb =  new JU.SB ();
JS.SpaceGroup.getSpaceGroups ();
for (var i = 0; i < JS.SpaceGroup.SG.length; i++) sb.append ("\n----------------------\n" + JS.SpaceGroup.SG[i].dumpInfo (null));

return sb.toString ();
});
c$.dumpAllSeitz = Clazz_defineMethod (c$, "dumpAllSeitz", 
 function () {
JS.SpaceGroup.getSpaceGroups ();
var sb =  new JU.SB ();
for (var i = 0; i < JS.SpaceGroup.SG.length; i++) sb.append ("\n").appendO (JS.SpaceGroup.SG[i].dumpCanonicalSeitzList ());

return sb.toString ();
});
Clazz_defineMethod (c$, "setLattice", 
 function (latticeCode, isCentrosymmetric) {
this.latticeCode = latticeCode;
this.latticeParameter = JS.HallTranslation.getLatticeIndex (latticeCode);
if (!isCentrosymmetric) this.latticeParameter = -this.latticeParameter;
}, "~S,~B");
c$.createSpaceGroupN = Clazz_defineMethod (c$, "createSpaceGroupN", 
 function (name) {
JS.SpaceGroup.getSpaceGroups ();
name = name.trim ();
var sg = JS.SpaceGroup.determineSpaceGroupN (name);
var hallInfo;
if (sg == null) {
hallInfo =  new JS.HallInfo (name);
if (hallInfo.nRotations > 0) {
sg =  new JS.SpaceGroup ("0;--;--;" + name, true);
sg.hallInfo = hallInfo;
} else if (name.indexOf (",") >= 0) {
sg =  new JS.SpaceGroup ("0;--;--;--", true);
sg.doNormalize = false;
sg.generateOperatorsFromXyzInfo (name);
}}if (sg != null) sg.generateAllOperators (null);
return sg;
}, "~S");
Clazz_defineMethod (c$, "addOperation", 
 function (xyz0, opId, allowScaling) {
if (xyz0 == null || xyz0.length < 3) {
this.xyzList =  new java.util.Hashtable ();
this.operationCount = 0;
return -1;
}var isSpecial = (xyz0.charAt (0) == '=');
if (isSpecial) xyz0 = xyz0.substring (1);
if (this.xyzList.containsKey (xyz0)) return this.xyzList.get (xyz0).intValue ();
if (xyz0.startsWith ("x1,x2,x3,x4") && this.modDim == 0) {
this.xyzList.clear ();
this.operationCount = 0;
this.modDim = JU.PT.parseInt (xyz0.substring (xyz0.lastIndexOf ("x") + 1)) - 3;
} else if (xyz0.equals ("x,y,z,m+1")) {
this.xyzList.clear ();
this.operationCount = 0;
}var op =  new JS.SymmetryOperation (null, null, 0, opId, this.doNormalize);
if (!op.setMatrixFromXYZ (xyz0, this.modDim, allowScaling)) {
JU.Logger.error ("couldn't interpret symmetry operation: " + xyz0);
return -1;
}return this.addOp (op, xyz0, isSpecial);
}, "~S,~N,~B");
Clazz_defineMethod (c$, "addOp", 
 function (op, xyz0, isSpecial) {
var ext = "";
var xyz = op.xyz + ext;
var xxx = JU.PT.replaceAllCharacters (xyz, "+123/", "");
if (!isSpecial) {
if (this.xyzList.containsKey (xyz)) return this.xyzList.get (xyz).intValue ();
if (this.latticeOp < 0) {
if (this.xyzList.containsKey (xxx)) this.latticeOp = this.operationCount;
 else this.xyzList.put (xxx, Integer.$valueOf (this.operationCount));
}this.xyzList.put (xyz, Integer.$valueOf (this.operationCount));
}if (!xyz.equals (xyz0 + ext)) this.xyzList.put (xyz0 + ext, Integer.$valueOf (this.operationCount));
if (this.operations == null) this.operations =  new Array (4);
if (this.operationCount == this.operations.length) this.operations = JU.AU.arrayCopyObject (this.operations, this.operationCount * 2);
this.operations[this.operationCount++] = op;
op.index = this.operationCount;
if (JU.Logger.debugging) JU.Logger.debug ("\naddOperation " + this.operationCount + op.dumpInfo ());
return this.operationCount - 1;
}, "JS.SymmetryOperation,~S,~B");
Clazz_defineMethod (c$, "generateOperatorsFromXyzInfo", 
 function (xyzInfo) {
this.init ();
var terms = JU.PT.split (xyzInfo.toLowerCase (), ";");
for (var i = 0; i < terms.length; i++) this.addSymmetry (terms[i], 0, false);

}, "~S");
Clazz_defineMethod (c$, "generateAllOperators", 
 function (h) {
if (h == null) {
h = this.hallInfo;
if (this.operationCount > 0) return;
this.operations =  new Array (4);
if (this.hallInfo == null || this.hallInfo.nRotations == 0) h = this.hallInfo =  new JS.HallInfo (this.hallSymbol);
this.setLattice (this.hallInfo.latticeCode, this.hallInfo.isCentrosymmetric);
this.init ();
}var mat1 =  new JU.M4 ();
var operation =  new JU.M4 ();
var newOps =  new Array (7);
for (var i = 0; i < 7; i++) newOps[i] =  new JU.M4 ();

for (var i = 0; i < h.nRotations; i++) {
mat1.setM4 (h.rotationTerms[i].seitzMatrix12ths);
var nRot = h.rotationTerms[i].order;
newOps[0].setIdentity ();
var nOps = this.operationCount;
for (var j = 1; j <= nRot; j++) {
newOps[j].mul2 (mat1, newOps[0]);
newOps[0].setM4 (newOps[j]);
for (var k = 0; k < nOps; k++) {
operation.mul2 (newOps[j], this.operations[k]);
JS.SymmetryOperation.normalizeTranslation (operation);
var xyz = JS.SymmetryOperation.getXYZFromMatrix (operation, true, true, true);
this.addSymmetrySM (xyz, operation);
}
}
}
}, "JS.HallInfo");
Clazz_defineMethod (c$, "addSymmetrySM", 
function (xyz, operation) {
var iop = this.addOperation (xyz, 0, false);
if (iop >= 0) {
var symmetryOperation = this.operations[iop];
symmetryOperation.setM4 (operation);
}return iop;
}, "~S,JU.M4");
c$.determineSpaceGroupN = Clazz_defineMethod (c$, "determineSpaceGroupN", 
 function (name) {
return JS.SpaceGroup.determineSpaceGroup (name, 0, 0, 0, 0, 0, 0, -1);
}, "~S");
c$.determineSpaceGroupNS = Clazz_defineMethod (c$, "determineSpaceGroupNS", 
 function (name, sg) {
return JS.SpaceGroup.determineSpaceGroup (name, 0, 0, 0, 0, 0, 0, sg.index);
}, "~S,JS.SpaceGroup");
c$.determineSpaceGroupNA = Clazz_defineMethod (c$, "determineSpaceGroupNA", 
 function (name, notionalUnitcell) {
return (notionalUnitcell == null ? JS.SpaceGroup.determineSpaceGroup (name, 0, 0, 0, 0, 0, 0, -1) : JS.SpaceGroup.determineSpaceGroup (name, notionalUnitcell[0], notionalUnitcell[1], notionalUnitcell[2], notionalUnitcell[3], notionalUnitcell[4], notionalUnitcell[5], -1));
}, "~S,~A");
c$.determineSpaceGroup = Clazz_defineMethod (c$, "determineSpaceGroup", 
 function (name, a, b, c, alpha, beta, gamma, lastIndex) {
var i = JS.SpaceGroup.determineSpaceGroupIndex (name, a, b, c, alpha, beta, gamma, lastIndex);
return (i >= 0 ? JS.SpaceGroup.SG[i] : null);
}, "~S,~N,~N,~N,~N,~N,~N,~N");
c$.determineSpaceGroupIndex = Clazz_defineMethod (c$, "determineSpaceGroupIndex", 
 function (name, a, b, c, alpha, beta, gamma, lastIndex) {
JS.SpaceGroup.getSpaceGroups ();
if (lastIndex < 0) lastIndex = JS.SpaceGroup.SG.length;
name = name.trim ().toLowerCase ();
var checkBilbao = false;
if (name.startsWith ("bilbao:")) {
checkBilbao = true;
name = name.substring (7);
}var nameType = (name.startsWith ("hall:") ? 5 : name.startsWith ("hm:") ? 3 : 0);
if (nameType > 0) name = name.substring (nameType);
 else if (name.contains ("[")) {
nameType = 5;
name = name.substring (0, name.indexOf ("[")).trim ();
}var nameExt = name;
var i;
var haveExtension = false;
name = name.$replace ('_', ' ');
if (name.length >= 2) {
i = (name.indexOf ("-") == 0 ? 2 : 1);
if (i < name.length && name.charAt (i) != ' ') name = name.substring (0, i) + " " + name.substring (i);
name = name.substring (0, 2).toUpperCase () + name.substring (2);
}var ext = "";
if ((i = name.indexOf (":")) > 0) {
ext = name.substring (i + 1);
name = name.substring (0, i).trim ();
haveExtension = true;
}if (nameType != 5 && !haveExtension && JU.PT.isOneOf (name, JS.SpaceGroup.ambiguousNames)) {
ext = "?";
haveExtension = true;
}var abbr = JU.PT.replaceAllCharacters (name, " ()", "");
var s;
if (nameType != 3 && !haveExtension) for (i = lastIndex; --i >= 0; ) {
s = JS.SpaceGroup.SG[i];
if (s.hallSymbol.equals (name)) return i;
}
if (nameType != 5) {
if (nameType != 3) for (i = lastIndex; --i >= 0; ) {
s = JS.SpaceGroup.SG[i];
if (s.intlTableNumberFull.equals (nameExt)) return i;
}
for (i = lastIndex; --i >= 0; ) {
s = JS.SpaceGroup.SG[i];
if (s.hmSymbolFull.equals (nameExt)) return i;
}
for (i = lastIndex; --i >= 0; ) {
s = JS.SpaceGroup.SG[i];
if (s.hmSymbolAlternative != null && s.hmSymbolAlternative.equals (nameExt)) return i;
}
if (haveExtension) for (i = lastIndex; --i >= 0; ) {
s = JS.SpaceGroup.SG[i];
if (s.hmSymbolAbbr.equals (abbr) && s.intlTableNumberExt.equals (ext)) return i;
}
if (haveExtension) for (i = lastIndex; --i >= 0; ) {
s = JS.SpaceGroup.SG[i];
if (s.hmSymbolAbbrShort.equals (abbr) && s.intlTableNumberExt.equals (ext)) return i;
}
var uniqueAxis = JS.SpaceGroup.determineUniqueAxis (a, b, c, alpha, beta, gamma);
if (!haveExtension || ext.charAt (0) == '?') for (i = lastIndex; --i >= 0; ) {
s = JS.SpaceGroup.SG[i];
if (s.hmSymbolAbbr.equals (abbr) || s.hmSymbolAbbrShort.equals (abbr)) {
switch (s.ambiguityType) {
case '\0':
return i;
case 'a':
if (s.uniqueAxis == uniqueAxis || uniqueAxis == '\0') return i;
break;
case 'o':
if (ext.length == 0) {
if (s.hmSymbolExt.equals ("2")) return i;
} else if (s.hmSymbolExt.equals (ext)) return i;
break;
case 't':
if (ext.length == 0) {
if (s.axisChoice == 'h') return i;
} else if ((s.axisChoice + "").equals (ext)) return i;
break;
}
}}
}if (ext.length == 0) for (i = 0; i < lastIndex; i++) {
s = JS.SpaceGroup.SG[i];
if (s.intlTableNumber.equals (nameExt) && (!checkBilbao || s.isBilbao)) return i;
}
return -1;
}, "~S,~N,~N,~N,~N,~N,~N,~N");
c$.determineUniqueAxis = Clazz_defineMethod (c$, "determineUniqueAxis", 
 function (a, b, c, alpha, beta, gamma) {
if (a == b) return (b == c ? '\0' : 'c');
if (b == c) return 'a';
if (c == a) return 'b';
if (alpha == beta) return (beta == gamma ? '\0' : 'c');
if (beta == gamma) return 'a';
if (gamma == alpha) return 'b';
return '\0';
}, "~N,~N,~N,~N,~N,~N");
Clazz_defineMethod (c$, "buildSpaceGroup", 
 function (cifLine) {
var terms = JU.PT.split (cifLine.toLowerCase (), ";");
this.isBilbao = (terms.length < 5);
this.intlTableNumberFull = terms[0].trim ();
var parts = JU.PT.split (this.intlTableNumberFull, ":");
this.intlTableNumber = parts[0];
this.intlTableNumberExt = (parts.length == 1 ? "" : parts[1]);
this.ambiguityType = '\0';
if (this.intlTableNumberExt.length > 0) {
var term = this.intlTableNumberExt;
if (term.startsWith ("-")) term = term.substring (1);
if (term.equals ("h") || term.equals ("r")) {
this.ambiguityType = 't';
this.axisChoice = this.intlTableNumberExt.charAt (0);
} else if (this.intlTableNumberExt.startsWith ("1") || this.intlTableNumberExt.startsWith ("2")) {
this.ambiguityType = 'o';
} else if (this.intlTableNumberExt.length <= 2) {
this.ambiguityType = 'a';
this.uniqueAxis = this.intlTableNumberExt.charAt (0);
}}this.hmSymbolFull = Character.toUpperCase (terms[2].charAt (0)) + terms[2].substring (1);
parts = JU.PT.split (this.hmSymbolFull, ":");
this.hmSymbol = parts[0];
this.hmSymbolExt = (parts.length == 1 ? "" : parts[1]);
var pt = this.hmSymbol.indexOf (" -3");
if (pt >= 1) if ("admn".indexOf (this.hmSymbol.charAt (pt - 1)) >= 0) {
this.hmSymbolAlternative = (this.hmSymbol.substring (0, pt) + " 3" + this.hmSymbol.substring (pt + 3)).toLowerCase ();
}this.hmSymbolAbbr = JU.PT.rep (this.hmSymbol, " ", "");
this.hmSymbolAbbrShort = JU.PT.rep (this.hmSymbol, " 1", "");
this.hmSymbolAbbrShort = JU.PT.rep (this.hmSymbolAbbrShort, " ", "");
this.hallSymbol = terms[3];
if (this.hallSymbol.length > 1) this.hallSymbol = this.hallSymbol.substring (0, 2).toUpperCase () + this.hallSymbol.substring (2);
var info = this.intlTableNumber + this.hallSymbol;
if (this.intlTableNumber.charAt (0) != '0' && JS.SpaceGroup.lastInfo.equals (info)) JS.SpaceGroup.ambiguousNames += this.hmSymbol + ";";
JS.SpaceGroup.lastInfo = info;
this.name = this.hallSymbol + " [" + this.hmSymbolFull + "] #" + this.intlTableNumber;
}, "~S");
c$.getSpaceGroups = Clazz_defineMethod (c$, "getSpaceGroups", 
 function () {
return (JS.SpaceGroup.SG == null ? (JS.SpaceGroup.SG = JS.SpaceGroup.createSpaceGroups ()) : JS.SpaceGroup.SG);
});
c$.createSpaceGroups = Clazz_defineMethod (c$, "createSpaceGroups", 
 function () {
var n = JS.SpaceGroup.STR_SG.length;
var defs =  new Array (n);
for (var i = 0; i < n; i++) defs[i] =  new JS.SpaceGroup (JS.SpaceGroup.STR_SG[i], true);

JS.SpaceGroup.STR_SG = null;
return defs;
});
Clazz_defineMethod (c$, "addLatticeVectors", 
function (lattvecs) {
if (this.latticeOp >= 0) return false;
var nOps = this.latticeOp = this.operationCount;
for (var j = 0; j < lattvecs.size (); j++) {
var data = lattvecs.get (j);
var magRev = Clazz_floatToInt (data.length == 5 && Float.isNaN (data[4]) ? data[3] : -2);
if (magRev != -2) data = [data[0], data[1], data[2]];
if (data.length > this.modDim + 3) return false;
for (var i = 0; i < nOps; i++) {
var newOp =  new JS.SymmetryOperation (null, null, 0, 0, this.doNormalize);
newOp.modDim = this.modDim;
var op = this.operations[i];
newOp.linearRotTrans = JU.AU.arrayCopyF (op.linearRotTrans, -1);
newOp.setFromMatrix (data, false);
newOp.xyzOriginal = newOp.xyz;
if (magRev != -2) newOp.setTimeReversal (op.timeReversal * magRev);
this.addOp (newOp, newOp.xyz, true);
}
}
return true;
}, "JU.Lst");
Clazz_defineMethod (c$, "getSiteMultiplicity", 
function (pt, unitCell) {
var n = this.finalOperations.length;
var pts =  new JU.Lst ();
for (var i = n; --i >= 0; ) {
var pt1 = JU.P3.newP (pt);
this.finalOperations[i].rotTrans (pt1);
unitCell.unitize (pt1);
for (var j = pts.size (); --j >= 0; ) {
var pt0 = pts.get (j);
if (pt1.distanceSquared (pt0) < 0.000001) {
pt1 = null;
break;
}}
if (pt1 != null) pts.addLast (pt1);
}
return Clazz_doubleToInt (n / pts.size ());
}, "JU.P3,JS.UnitCell");
Clazz_defineMethod (c$, "getAllLatticeOps", 
function () {
if (this.latticeOp < 0 || this.modDim > 0) return null;
if (this.latticeOps == null) {
this.latticeOps =  Clazz_newIntArray (3, 0);
var nOps = 0;
for (var i = this.latticeOp; i < this.operationCount; i++) {
var o = this.finalOperations[i];
if (o.m00 + o.m01 + o.m02 == 3) {
System.out.println ("spacegroup " + o);
this.latticeOps[nOps++] = i;
}}
}return this.latticeOps;
});
Clazz_defineStatics (c$,
"canonicalSeitzList", null,
"NAME_HALL", 5,
"NAME_HM", 3,
"sgIndex", -1,
"ambiguousNames", "",
"lastInfo", "",
"SG", null,
"STR_SG", ["1;c1^1;p 1;p 1", "2;ci^1;p -1;-p 1", "3:b;c2^1;p 1 2 1;p 2y", "3:b;c2^1;p 2;p 2y", "3:c;c2^1;p 1 1 2;p 2", "3:a;c2^1;p 2 1 1;p 2x", "4:b;c2^2;p 1 21 1;p 2yb", "4:b;c2^2;p 21;p 2yb", "4:b*;c2^2;p 1 21 1*;p 2y1", "4:c;c2^2;p 1 1 21;p 2c", "4:c*;c2^2;p 1 1 21*;p 21", "4:a;c2^2;p 21 1 1;p 2xa", "4:a*;c2^2;p 21 1 1*;p 2x1", "5:b1;c2^3;c 1 2 1;c 2y", "5:b1;c2^3;c 2;c 2y", "5:b2;c2^3;a 1 2 1;a 2y", "5:b3;c2^3;i 1 2 1;i 2y", "5:c1;c2^3;a 1 1 2;a 2", "5:c2;c2^3;b 1 1 2;b 2", "5:c3;c2^3;i 1 1 2;i 2", "5:a1;c2^3;b 2 1 1;b 2x", "5:a2;c2^3;c 2 1 1;c 2x", "5:a3;c2^3;i 2 1 1;i 2x", "6:b;cs^1;p 1 m 1;p -2y", "6:b;cs^1;p m;p -2y", "6:c;cs^1;p 1 1 m;p -2", "6:a;cs^1;p m 1 1;p -2x", "7:b1;cs^2;p 1 c 1;p -2yc", "7:b1;cs^2;p c;p -2yc", "7:b2;cs^2;p 1 n 1;p -2yac", "7:b2;cs^2;p n;p -2yac", "7:b3;cs^2;p 1 a 1;p -2ya", "7:b3;cs^2;p a;p -2ya", "7:c1;cs^2;p 1 1 a;p -2a", "7:c2;cs^2;p 1 1 n;p -2ab", "7:c3;cs^2;p 1 1 b;p -2b", "7:a1;cs^2;p b 1 1;p -2xb", "7:a2;cs^2;p n 1 1;p -2xbc", "7:a3;cs^2;p c 1 1;p -2xc", "8:b1;cs^3;c 1 m 1;c -2y", "8:b1;cs^3;c m;c -2y", "8:b2;cs^3;a 1 m 1;a -2y", "8:b3;cs^3;i 1 m 1;i -2y", "8:b3;cs^3;i m;i -2y", "8:c1;cs^3;a 1 1 m;a -2", "8:c2;cs^3;b 1 1 m;b -2", "8:c3;cs^3;i 1 1 m;i -2", "8:a1;cs^3;b m 1 1;b -2x", "8:a2;cs^3;c m 1 1;c -2x", "8:a3;cs^3;i m 1 1;i -2x", "9:b1;cs^4;c 1 c 1;c -2yc", "9:b1;cs^4;c c;c -2yc", "9:b2;cs^4;a 1 n 1;a -2yab", "9:b3;cs^4;i 1 a 1;i -2ya", "9:-b1;cs^4;a 1 a 1;a -2ya", "9:-b2;cs^4;c 1 n 1;c -2yac", "9:-b3;cs^4;i 1 c 1;i -2yc", "9:c1;cs^4;a 1 1 a;a -2a", "9:c2;cs^4;b 1 1 n;b -2ab", "9:c3;cs^4;i 1 1 b;i -2b", "9:-c1;cs^4;b 1 1 b;b -2b", "9:-c2;cs^4;a 1 1 n;a -2ab", "9:-c3;cs^4;i 1 1 a;i -2a", "9:a1;cs^4;b b 1 1;b -2xb", "9:a2;cs^4;c n 1 1;c -2xac", "9:a3;cs^4;i c 1 1;i -2xc", "9:-a1;cs^4;c c 1 1;c -2xc", "9:-a2;cs^4;b n 1 1;b -2xab", "9:-a3;cs^4;i b 1 1;i -2xb", "10:b;c2h^1;p 1 2/m 1;-p 2y", "10:b;c2h^1;p 2/m;-p 2y", "10:c;c2h^1;p 1 1 2/m;-p 2", "10:a;c2h^1;p 2/m 1 1;-p 2x", "11:b;c2h^2;p 1 21/m 1;-p 2yb", "11:b;c2h^2;p 21/m;-p 2yb", "11:b*;c2h^2;p 1 21/m 1*;-p 2y1", "11:c;c2h^2;p 1 1 21/m;-p 2c", "11:c*;c2h^2;p 1 1 21/m*;-p 21", "11:a;c2h^2;p 21/m 1 1;-p 2xa", "11:a*;c2h^2;p 21/m 1 1*;-p 2x1", "12:b1;c2h^3;c 1 2/m 1;-c 2y", "12:b1;c2h^3;c 2/m;-c 2y", "12:b2;c2h^3;a 1 2/m 1;-a 2y", "12:b3;c2h^3;i 1 2/m 1;-i 2y", "12:b3;c2h^3;i 2/m;-i 2y", "12:c1;c2h^3;a 1 1 2/m;-a 2", "12:c2;c2h^3;b 1 1 2/m;-b 2", "12:c3;c2h^3;i 1 1 2/m;-i 2", "12:a1;c2h^3;b 2/m 1 1;-b 2x", "12:a2;c2h^3;c 2/m 1 1;-c 2x", "12:a3;c2h^3;i 2/m 1 1;-i 2x", "13:b1;c2h^4;p 1 2/c 1;-p 2yc", "13:b1;c2h^4;p 2/c;-p 2yc", "13:b2;c2h^4;p 1 2/n 1;-p 2yac", "13:b2;c2h^4;p 2/n;-p 2yac", "13:b3;c2h^4;p 1 2/a 1;-p 2ya", "13:b3;c2h^4;p 2/a;-p 2ya", "13:c1;c2h^4;p 1 1 2/a;-p 2a", "13:c2;c2h^4;p 1 1 2/n;-p 2ab", "13:c3;c2h^4;p 1 1 2/b;-p 2b", "13:a1;c2h^4;p 2/b 1 1;-p 2xb", "13:a2;c2h^4;p 2/n 1 1;-p 2xbc", "13:a3;c2h^4;p 2/c 1 1;-p 2xc", "14:b1;c2h^5;p 1 21/c 1;-p 2ybc", "14:b1;c2h^5;p 21/c;-p 2ybc", "14:b2;c2h^5;p 1 21/n 1;-p 2yn", "14:b2;c2h^5;p 21/n;-p 2yn", "14:b3;c2h^5;p 1 21/a 1;-p 2yab", "14:b3;c2h^5;p 21/a;-p 2yab", "14:c1;c2h^5;p 1 1 21/a;-p 2ac", "14:c2;c2h^5;p 1 1 21/n;-p 2n", "14:c3;c2h^5;p 1 1 21/b;-p 2bc", "14:a1;c2h^5;p 21/b 1 1;-p 2xab", "14:a2;c2h^5;p 21/n 1 1;-p 2xn", "14:a3;c2h^5;p 21/c 1 1;-p 2xac", "15:b1;c2h^6;c 1 2/c 1;-c 2yc", "15:b1;c2h^6;c 2/c;-c 2yc", "15:b2;c2h^6;a 1 2/n 1;-a 2yab", "15:b3;c2h^6;i 1 2/a 1;-i 2ya", "15:b3;c2h^6;i 2/a;-i 2ya", "15:-b1;c2h^6;a 1 2/a 1;-a 2ya", "15:-b2;c2h^6;c 1 2/n 1;-c 2yac", "15:-b2;c2h^6;c 2/n;-c 2yac", "15:-b3;c2h^6;i 1 2/c 1;-i 2yc", "15:-b3;c2h^6;i 2/c;-i 2yc", "15:c1;c2h^6;a 1 1 2/a;-a 2a", "15:c2;c2h^6;b 1 1 2/n;-b 2ab", "15:c3;c2h^6;i 1 1 2/b;-i 2b", "15:-c1;c2h^6;b 1 1 2/b;-b 2b", "15:-c2;c2h^6;a 1 1 2/n;-a 2ab", "15:-c3;c2h^6;i 1 1 2/a;-i 2a", "15:a1;c2h^6;b 2/b 1 1;-b 2xb", "15:a2;c2h^6;c 2/n 1 1;-c 2xac", "15:a3;c2h^6;i 2/c 1 1;-i 2xc", "15:-a1;c2h^6;c 2/c 1 1;-c 2xc", "15:-a2;c2h^6;b 2/n 1 1;-b 2xab", "15:-a3;c2h^6;i 2/b 1 1;-i 2xb", "16;d2^1;p 2 2 2;p 2 2", "17;d2^2;p 2 2 21;p 2c 2", "17*;d2^2;p 2 2 21*;p 21 2", "17:cab;d2^2;p 21 2 2;p 2a 2a", "17:bca;d2^2;p 2 21 2;p 2 2b", "18;d2^3;p 21 21 2;p 2 2ab", "18:cab;d2^3;p 2 21 21;p 2bc 2", "18:bca;d2^3;p 21 2 21;p 2ac 2ac", "19;d2^4;p 21 21 21;p 2ac 2ab", "20;d2^5;c 2 2 21;c 2c 2", "20*;d2^5;c 2 2 21*;c 21 2", "20:cab;d2^5;a 21 2 2;a 2a 2a", "20:cab*;d2^5;a 21 2 2*;a 2a 21", "20:bca;d2^5;b 2 21 2;b 2 2b", "21;d2^6;c 2 2 2;c 2 2", "21:cab;d2^6;a 2 2 2;a 2 2", "21:bca;d2^6;b 2 2 2;b 2 2", "22;d2^7;f 2 2 2;f 2 2", "23;d2^8;i 2 2 2;i 2 2", "24;d2^9;i 21 21 21;i 2b 2c", "25;c2v^1;p m m 2;p 2 -2", "25:cab;c2v^1;p 2 m m;p -2 2", "25:bca;c2v^1;p m 2 m;p -2 -2", "26;c2v^2;p m c 21;p 2c -2", "26*;c2v^2;p m c 21*;p 21 -2", "26:ba-c;c2v^2;p c m 21;p 2c -2c", "26:ba-c*;c2v^2;p c m 21*;p 21 -2c", "26:cab;c2v^2;p 21 m a;p -2a 2a", "26:-cba;c2v^2;p 21 a m;p -2 2a", "26:bca;c2v^2;p b 21 m;p -2 -2b", "26:a-cb;c2v^2;p m 21 b;p -2b -2", "27;c2v^3;p c c 2;p 2 -2c", "27:cab;c2v^3;p 2 a a;p -2a 2", "27:bca;c2v^3;p b 2 b;p -2b -2b", "28;c2v^4;p m a 2;p 2 -2a", "28*;c2v^4;p m a 2*;p 2 -21", "28:ba-c;c2v^4;p b m 2;p 2 -2b", "28:cab;c2v^4;p 2 m b;p -2b 2", "28:-cba;c2v^4;p 2 c m;p -2c 2", "28:-cba*;c2v^4;p 2 c m*;p -21 2", "28:bca;c2v^4;p c 2 m;p -2c -2c", "28:a-cb;c2v^4;p m 2 a;p -2a -2a", "29;c2v^5;p c a 21;p 2c -2ac", "29:ba-c;c2v^5;p b c 21;p 2c -2b", "29:cab;c2v^5;p 21 a b;p -2b 2a", "29:-cba;c2v^5;p 21 c a;p -2ac 2a", "29:bca;c2v^5;p c 21 b;p -2bc -2c", "29:a-cb;c2v^5;p b 21 a;p -2a -2ab", "30;c2v^6;p n c 2;p 2 -2bc", "30:ba-c;c2v^6;p c n 2;p 2 -2ac", "30:cab;c2v^6;p 2 n a;p -2ac 2", "30:-cba;c2v^6;p 2 a n;p -2ab 2", "30:bca;c2v^6;p b 2 n;p -2ab -2ab", "30:a-cb;c2v^6;p n 2 b;p -2bc -2bc", "31;c2v^7;p m n 21;p 2ac -2", "31:ba-c;c2v^7;p n m 21;p 2bc -2bc", "31:cab;c2v^7;p 21 m n;p -2ab 2ab", "31:-cba;c2v^7;p 21 n m;p -2 2ac", "31:bca;c2v^7;p n 21 m;p -2 -2bc", "31:a-cb;c2v^7;p m 21 n;p -2ab -2", "32;c2v^8;p b a 2;p 2 -2ab", "32:cab;c2v^8;p 2 c b;p -2bc 2", "32:bca;c2v^8;p c 2 a;p -2ac -2ac", "33;c2v^9;p n a 21;p 2c -2n", "33*;c2v^9;p n a 21*;p 21 -2n", "33:ba-c;c2v^9;p b n 21;p 2c -2ab", "33:ba-c*;c2v^9;p b n 21*;p 21 -2ab", "33:cab;c2v^9;p 21 n b;p -2bc 2a", "33:cab*;c2v^9;p 21 n b*;p -2bc 21", "33:-cba;c2v^9;p 21 c n;p -2n 2a", "33:-cba*;c2v^9;p 21 c n*;p -2n 21", "33:bca;c2v^9;p c 21 n;p -2n -2ac", "33:a-cb;c2v^9;p n 21 a;p -2ac -2n", "34;c2v^10;p n n 2;p 2 -2n", "34:cab;c2v^10;p 2 n n;p -2n 2", "34:bca;c2v^10;p n 2 n;p -2n -2n", "35;c2v^11;c m m 2;c 2 -2", "35:cab;c2v^11;a 2 m m;a -2 2", "35:bca;c2v^11;b m 2 m;b -2 -2", "36;c2v^12;c m c 21;c 2c -2", "36*;c2v^12;c m c 21*;c 21 -2", "36:ba-c;c2v^12;c c m 21;c 2c -2c", "36:ba-c*;c2v^12;c c m 21*;c 21 -2c", "36:cab;c2v^12;a 21 m a;a -2a 2a", "36:cab*;c2v^12;a 21 m a*;a -2a 21", "36:-cba;c2v^12;a 21 a m;a -2 2a", "36:-cba*;c2v^12;a 21 a m*;a -2 21", "36:bca;c2v^12;b b 21 m;b -2 -2b", "36:a-cb;c2v^12;b m 21 b;b -2b -2", "37;c2v^13;c c c 2;c 2 -2c", "37:cab;c2v^13;a 2 a a;a -2a 2", "37:bca;c2v^13;b b 2 b;b -2b -2b", "38;c2v^14;a m m 2;a 2 -2", "38:ba-c;c2v^14;b m m 2;b 2 -2", "38:cab;c2v^14;b 2 m m;b -2 2", "38:-cba;c2v^14;c 2 m m;c -2 2", "38:bca;c2v^14;c m 2 m;c -2 -2", "38:a-cb;c2v^14;a m 2 m;a -2 -2", "39;c2v^15;a e m 2;a 2 -2b", "39;c2v^15;a b m 2;a 2 -2b", "39:ba-c;c2v^15;b m a 2;b 2 -2a", "39:cab;c2v^15;b 2 c m;b -2a 2", "39:-cba;c2v^15;c 2 m b;c -2a 2", "39:bca;c2v^15;c m 2 a;c -2a -2a", "39:a-cb;c2v^15;a c 2 m;a -2b -2b", "40;c2v^16;a m a 2;a 2 -2a", "40:ba-c;c2v^16;b b m 2;b 2 -2b", "40:cab;c2v^16;b 2 m b;b -2b 2", "40:-cba;c2v^16;c 2 c m;c -2c 2", "40:bca;c2v^16;c c 2 m;c -2c -2c", "40:a-cb;c2v^16;a m 2 a;a -2a -2a", "41;c2v^17;a e a 2;a 2 -2ab", "41;c2v^17;a b a 2;a 2 -2ab;-b", "41:ba-c;c2v^17;b b a 2;b 2 -2ab", "41:cab;c2v^17;b 2 c b;b -2ab 2", "41:-cba;c2v^17;c 2 c b;c -2ac 2", "41:bca;c2v^17;c c 2 a;c -2ac -2ac", "41:a-cb;c2v^17;a c 2 a;a -2ab -2ab", "42;c2v^18;f m m 2;f 2 -2", "42:cab;c2v^18;f 2 m m;f -2 2", "42:bca;c2v^18;f m 2 m;f -2 -2", "43;c2v^19;f d d 2;f 2 -2d", "43:cab;c2v^19;f 2 d d;f -2d 2", "43:bca;c2v^19;f d 2 d;f -2d -2d", "44;c2v^20;i m m 2;i 2 -2", "44:cab;c2v^20;i 2 m m;i -2 2", "44:bca;c2v^20;i m 2 m;i -2 -2", "45;c2v^21;i b a 2;i 2 -2c", "45:cab;c2v^21;i 2 c b;i -2a 2", "45:bca;c2v^21;i c 2 a;i -2b -2b", "46;c2v^22;i m a 2;i 2 -2a", "46:ba-c;c2v^22;i b m 2;i 2 -2b", "46:cab;c2v^22;i 2 m b;i -2b 2", "46:-cba;c2v^22;i 2 c m;i -2c 2", "46:bca;c2v^22;i c 2 m;i -2c -2c", "46:a-cb;c2v^22;i m 2 a;i -2a -2a", "47;d2h^1;p m m m;-p 2 2", "48:1;d2h^2;p n n n:1;p 2 2 -1n;-b", "48:2;d2h^2;p n n n:2;-p 2ab 2bc", "49;d2h^3;p c c m;-p 2 2c", "49:cab;d2h^3;p m a a;-p 2a 2", "49:bca;d2h^3;p b m b;-p 2b 2b", "50:1;d2h^4;p b a n:1;p 2 2 -1ab;-b", "50:2;d2h^4;p b a n:2;-p 2ab 2b", "50:1cab;d2h^4;p n c b:1;p 2 2 -1bc", "50:2cab;d2h^4;p n c b:2;-p 2b 2bc", "50:1bca;d2h^4;p c n a:1;p 2 2 -1ac", "50:2bca;d2h^4;p c n a:2;-p 2a 2c", "51;d2h^5;p m m a;-p 2a 2a", "51:ba-c;d2h^5;p m m b;-p 2b 2", "51:cab;d2h^5;p b m m;-p 2 2b", "51:-cba;d2h^5;p c m m;-p 2c 2c", "51:bca;d2h^5;p m c m;-p 2c 2", "51:a-cb;d2h^5;p m a m;-p 2 2a", "52;d2h^6;p n n a;-p 2a 2bc", "52:ba-c;d2h^6;p n n b;-p 2b 2n", "52:cab;d2h^6;p b n n;-p 2n 2b", "52:-cba;d2h^6;p c n n;-p 2ab 2c", "52:bca;d2h^6;p n c n;-p 2ab 2n", "52:a-cb;d2h^6;p n a n;-p 2n 2bc", "53;d2h^7;p m n a;-p 2ac 2", "53:ba-c;d2h^7;p n m b;-p 2bc 2bc", "53:cab;d2h^7;p b m n;-p 2ab 2ab", "53:-cba;d2h^7;p c n m;-p 2 2ac", "53:bca;d2h^7;p n c m;-p 2 2bc", "53:a-cb;d2h^7;p m a n;-p 2ab 2", "54;d2h^8;p c c a;-p 2a 2ac", "54:ba-c;d2h^8;p c c b;-p 2b 2c", "54:cab;d2h^8;p b a a;-p 2a 2b", "54:-cba;d2h^8;p c a a;-p 2ac 2c", "54:bca;d2h^8;p b c b;-p 2bc 2b", "54:a-cb;d2h^8;p b a b;-p 2b 2ab", "55;d2h^9;p b a m;-p 2 2ab", "55:cab;d2h^9;p m c b;-p 2bc 2", "55:bca;d2h^9;p c m a;-p 2ac 2ac", "56;d2h^10;p c c n;-p 2ab 2ac", "56:cab;d2h^10;p n a a;-p 2ac 2bc", "56:bca;d2h^10;p b n b;-p 2bc 2ab", "57;d2h^11;p b c m;-p 2c 2b", "57:ba-c;d2h^11;p c a m;-p 2c 2ac", "57:cab;d2h^11;p m c a;-p 2ac 2a", "57:-cba;d2h^11;p m a b;-p 2b 2a", "57:bca;d2h^11;p b m a;-p 2a 2ab", "57:a-cb;d2h^11;p c m b;-p 2bc 2c", "58;d2h^12;p n n m;-p 2 2n", "58:cab;d2h^12;p m n n;-p 2n 2", "58:bca;d2h^12;p n m n;-p 2n 2n", "59:1;d2h^13;p m m n:1;p 2 2ab -1ab;-b", "59:2;d2h^13;p m m n:2;-p 2ab 2a", "59:1cab;d2h^13;p n m m:1;p 2bc 2 -1bc", "59:2cab;d2h^13;p n m m:2;-p 2c 2bc", "59:1bca;d2h^13;p m n m:1;p 2ac 2ac -1ac", "59:2bca;d2h^13;p m n m:2;-p 2c 2a", "60;d2h^14;p b c n;-p 2n 2ab", "60:ba-c;d2h^14;p c a n;-p 2n 2c", "60:cab;d2h^14;p n c a;-p 2a 2n", "60:-cba;d2h^14;p n a b;-p 2bc 2n", "60:bca;d2h^14;p b n a;-p 2ac 2b", "60:a-cb;d2h^14;p c n b;-p 2b 2ac", "61;d2h^15;p b c a;-p 2ac 2ab", "61:ba-c;d2h^15;p c a b;-p 2bc 2ac", "62;d2h^16;p n m a;-p 2ac 2n", "62:ba-c;d2h^16;p m n b;-p 2bc 2a", "62:cab;d2h^16;p b n m;-p 2c 2ab", "62:-cba;d2h^16;p c m n;-p 2n 2ac", "62:bca;d2h^16;p m c n;-p 2n 2a", "62:a-cb;d2h^16;p n a m;-p 2c 2n", "63;d2h^17;c m c m;-c 2c 2", "63:ba-c;d2h^17;c c m m;-c 2c 2c", "63:cab;d2h^17;a m m a;-a 2a 2a", "63:-cba;d2h^17;a m a m;-a 2 2a", "63:bca;d2h^17;b b m m;-b 2 2b", "63:a-cb;d2h^17;b m m b;-b 2b 2", "64;d2h^18;c m c e;-c 2ac 2", "64;d2h^18;c m c a;-c 2ac 2", "64:ba-c;d2h^18;c c m b;-c 2ac 2ac", "64:cab;d2h^18;a b m a;-a 2ab 2ab", "64:-cba;d2h^18;a c a m;-a 2 2ab", "64:bca;d2h^18;b b c m;-b 2 2ab", "64:a-cb;d2h^18;b m a b;-b 2ab 2", "65;d2h^19;c m m m;-c 2 2", "65:cab;d2h^19;a m m m;-a 2 2", "65:bca;d2h^19;b m m m;-b 2 2", "66;d2h^20;c c c m;-c 2 2c", "66:cab;d2h^20;a m a a;-a 2a 2", "66:bca;d2h^20;b b m b;-b 2b 2b", "67;d2h^21;c m m e;-c 2a 2", "67;d2h^21;c m m a;-c 2a 2", "67:ba-c;d2h^21;c m m b;-c 2a 2a", "67:cab;d2h^21;a b m m;-a 2b 2b", "67:-cba;d2h^21;a c m m;-a 2 2b", "67:bca;d2h^21;b m c m;-b 2 2a", "67:a-cb;d2h^21;b m a m;-b 2a 2", "68:1;d2h^22;c c c e:1;c 2 2 -1ac;-b", "68:1;d2h^22;c c c a:1;c 2 2 -1ac;-b", "68:2;d2h^22;c c c e:2;-c 2a 2ac", "68:2;d2h^22;c c c a:2;-c 2a 2ac", "68:1ba-c;d2h^22;c c c b:1;c 2 2 -1ac", "68:2ba-c;d2h^22;c c c b:2;-c 2a 2c", "68:1cab;d2h^22;a b a a:1;a 2 2 -1ab", "68:2cab;d2h^22;a b a a:2;-a 2a 2b", "68:1-cba;d2h^22;a c a a:1;a 2 2 -1ab", "68:2-cba;d2h^22;a c a a:2;-a 2ab 2b", "68:1bca;d2h^22;b b c b:1;b 2 2 -1ab", "68:2bca;d2h^22;b b c b:2;-b 2ab 2b", "68:1a-cb;d2h^22;b b a b:1;b 2 2 -1ab", "68:2a-cb;d2h^22;b b a b:2;-b 2b 2ab", "69;d2h^23;f m m m;-f 2 2", "70:1;d2h^24;f d d d:1;f 2 2 -1d;-b", "70:2;d2h^24;f d d d:2;-f 2uv 2vw", "71;d2h^25;i m m m;-i 2 2", "72;d2h^26;i b a m;-i 2 2c", "72:cab;d2h^26;i m c b;-i 2a 2", "72:bca;d2h^26;i c m a;-i 2b 2b", "73;d2h^27;i b c a;-i 2b 2c", "73:ba-c;d2h^27;i c a b;-i 2a 2b", "74;d2h^28;i m m a;-i 2b 2", "74:ba-c;d2h^28;i m m b;-i 2a 2a", "74:cab;d2h^28;i b m m;-i 2c 2c", "74:-cba;d2h^28;i c m m;-i 2 2b", "74:bca;d2h^28;i m c m;-i 2 2a", "74:a-cb;d2h^28;i m a m;-i 2c 2", "75;c4^1;p 4;p 4", "76;c4^2;p 41;p 4w", "76*;c4^2;p 41*;p 41", "77;c4^3;p 42;p 4c", "77*;c4^3;p 42*;p 42", "78;c4^4;p 43;p 4cw", "78*;c4^4;p 43*;p 43", "79;c4^5;i 4;i 4", "80;c4^6;i 41;i 4bw", "81;s4^1;p -4;p -4", "82;s4^2;i -4;i -4", "83;c4h^1;p 4/m;-p 4", "84;c4h^2;p 42/m;-p 4c", "84*;c4h^2;p 42/m*;-p 42", "85:1;c4h^3;p 4/n:1;p 4ab -1ab;-b", "85:2;c4h^3;p 4/n:2;-p 4a", "86:1;c4h^4;p 42/n:1;p 4n -1n;-b", "86:2;c4h^4;p 42/n:2;-p 4bc", "87;c4h^5;i 4/m;-i 4", "88:1;c4h^6;i 41/a:1;i 4bw -1bw;-b", "88:2;c4h^6;i 41/a:2;-i 4ad", "89;d4^1;p 4 2 2;p 4 2", "90;d4^2;p 4 21 2;p 4ab 2ab", "91;d4^3;p 41 2 2;p 4w 2c", "91*;d4^3;p 41 2 2*;p 41 2c", "92;d4^4;p 41 21 2;p 4abw 2nw", "93;d4^5;p 42 2 2;p 4c 2", "93*;d4^5;p 42 2 2*;p 42 2", "94;d4^6;p 42 21 2;p 4n 2n", "95;d4^7;p 43 2 2;p 4cw 2c", "95*;d4^7;p 43 2 2*;p 43 2c", "96;d4^8;p 43 21 2;p 4nw 2abw", "97;d4^9;i 4 2 2;i 4 2", "98;d4^10;i 41 2 2;i 4bw 2bw", "99;c4v^1;p 4 m m;p 4 -2", "100;c4v^2;p 4 b m;p 4 -2ab", "101;c4v^3;p 42 c m;p 4c -2c", "101*;c4v^3;p 42 c m*;p 42 -2c", "102;c4v^4;p 42 n m;p 4n -2n", "103;c4v^5;p 4 c c;p 4 -2c", "104;c4v^6;p 4 n c;p 4 -2n", "105;c4v^7;p 42 m c;p 4c -2", "105*;c4v^7;p 42 m c*;p 42 -2", "106;c4v^8;p 42 b c;p 4c -2ab", "106*;c4v^8;p 42 b c*;p 42 -2ab", "107;c4v^9;i 4 m m;i 4 -2", "108;c4v^10;i 4 c m;i 4 -2c", "109;c4v^11;i 41 m d;i 4bw -2", "110;c4v^12;i 41 c d;i 4bw -2c", "111;d2d^1;p -4 2 m;p -4 2", "112;d2d^2;p -4 2 c;p -4 2c", "113;d2d^3;p -4 21 m;p -4 2ab", "114;d2d^4;p -4 21 c;p -4 2n", "115;d2d^5;p -4 m 2;p -4 -2", "116;d2d^6;p -4 c 2;p -4 -2c", "117;d2d^7;p -4 b 2;p -4 -2ab", "118;d2d^8;p -4 n 2;p -4 -2n", "119;d2d^9;i -4 m 2;i -4 -2", "120;d2d^10;i -4 c 2;i -4 -2c", "121;d2d^11;i -4 2 m;i -4 2", "122;d2d^12;i -4 2 d;i -4 2bw", "123;d4h^1;p 4/m m m;-p 4 2", "124;d4h^2;p 4/m c c;-p 4 2c", "125:1;d4h^3;p 4/n b m:1;p 4 2 -1ab;-b", "125:2;d4h^3;p 4/n b m:2;-p 4a 2b", "126:1;d4h^4;p 4/n n c:1;p 4 2 -1n;-b", "126:2;d4h^4;p 4/n n c:2;-p 4a 2bc", "127;d4h^5;p 4/m b m;-p 4 2ab", "128;d4h^6;p 4/m n c;-p 4 2n", "129:1;d4h^7;p 4/n m m:1;p 4ab 2ab -1ab;-b", "129:2;d4h^7;p 4/n m m:2;-p 4a 2a", "130:1;d4h^8;p 4/n c c:1;p 4ab 2n -1ab;-b", "130:2;d4h^8;p 4/n c c:2;-p 4a 2ac", "131;d4h^9;p 42/m m c;-p 4c 2", "132;d4h^10;p 42/m c m;-p 4c 2c", "133:1;d4h^11;p 42/n b c:1;p 4n 2c -1n;-b", "133:2;d4h^11;p 42/n b c:2;-p 4ac 2b", "134:1;d4h^12;p 42/n n m:1;p 4n 2 -1n;-b", "134:2;d4h^12;p 42/n n m:2;-p 4ac 2bc", "135;d4h^13;p 42/m b c;-p 4c 2ab", "135*;d4h^13;p 42/m b c*;-p 42 2ab", "136;d4h^14;p 42/m n m;-p 4n 2n", "137:1;d4h^15;p 42/n m c:1;p 4n 2n -1n;-b", "137:2;d4h^15;p 42/n m c:2;-p 4ac 2a", "138:1;d4h^16;p 42/n c m:1;p 4n 2ab -1n;-b", "138:2;d4h^16;p 42/n c m:2;-p 4ac 2ac", "139;d4h^17;i 4/m m m;-i 4 2", "140;d4h^18;i 4/m c m;-i 4 2c", "141:1;d4h^19;i 41/a m d:1;i 4bw 2bw -1bw;-b", "141:2;d4h^19;i 41/a m d:2;-i 4bd 2", "142:1;d4h^20;i 41/a c d:1;i 4bw 2aw -1bw;-b", "142:2;d4h^20;i 41/a c d:2;-i 4bd 2c", "143;c3^1;p 3;p 3", "144;c3^2;p 31;p 31", "145;c3^3;p 32;p 32", "146:h;c3^4;r 3:h;r 3", "146:r;c3^4;r 3:r;p 3*", "147;c3i^1;p -3;-p 3", "148:h;c3i^2;r -3:h;-r 3", "148:r;c3i^2;r -3:r;-p 3*", "149;d3^1;p 3 1 2;p 3 2", "150;d3^2;p 3 2 1;p 3 2\"", "151;d3^3;p 31 1 2;p 31 2 (0 0 4)", "152;d3^4;p 31 2 1;p 31 2\"", "153;d3^5;p 32 1 2;p 32 2 (0 0 2)", "154;d3^6;p 32 2 1;p 32 2\"", "155:h;d3^7;r 3 2:h;r 3 2\"", "155:r;d3^7;r 3 2:r;p 3* 2", "156;c3v^1;p 3 m 1;p 3 -2\"", "157;c3v^2;p 3 1 m;p 3 -2", "158;c3v^3;p 3 c 1;p 3 -2\"c", "159;c3v^4;p 3 1 c;p 3 -2c", "160:h;c3v^5;r 3 m:h;r 3 -2\"", "160:r;c3v^5;r 3 m:r;p 3* -2", "161:h;c3v^6;r 3 c:h;r 3 -2\"c", "161:r;c3v^6;r 3 c:r;p 3* -2n", "162;d3d^1;p -3 1 m;-p 3 2", "163;d3d^2;p -3 1 c;-p 3 2c", "164;d3d^3;p -3 m 1;-p 3 2\"", "165;d3d^4;p -3 c 1;-p 3 2\"c", "166:h;d3d^5;r -3 m:h;-r 3 2\"", "166:r;d3d^5;r -3 m:r;-p 3* 2", "167:h;d3d^6;r -3 c:h;-r 3 2\"c", "167:r;d3d^6;r -3 c:r;-p 3* 2n", "168;c6^1;p 6;p 6", "169;c6^2;p 61;p 61", "170;c6^3;p 65;p 65", "171;c6^4;p 62;p 62", "172;c6^5;p 64;p 64", "173;c6^6;p 63;p 6c", "173*;c6^6;p 63*;p 63 ", "174;c3h^1;p -6;p -6", "175;c6h^1;p 6/m;-p 6", "176;c6h^2;p 63/m;-p 6c", "176*;c6h^2;p 63/m*;-p 63", "177;d6^1;p 6 2 2;p 6 2", "178;d6^2;p 61 2 2;p 61 2 (0 0 5)", "179;d6^3;p 65 2 2;p 65 2 (0 0 1)", "180;d6^4;p 62 2 2;p 62 2 (0 0 4)", "181;d6^5;p 64 2 2;p 64 2 (0 0 2)", "182;d6^6;p 63 2 2;p 6c 2c", "182*;d6^6;p 63 2 2*;p 63 2c", "183;c6v^1;p 6 m m;p 6 -2", "184;c6v^2;p 6 c c;p 6 -2c", "185;c6v^3;p 63 c m;p 6c -2", "185*;c6v^3;p 63 c m*;p 63 -2", "186;c6v^4;p 63 m c;p 6c -2c", "186*;c6v^4;p 63 m c*;p 63 -2c", "187;d3h^1;p -6 m 2;p -6 2", "188;d3h^2;p -6 c 2;p -6c 2", "189;d3h^3;p -6 2 m;p -6 -2", "190;d3h^4;p -6 2 c;p -6c -2c", "191;d6h^1;p 6/m m m;-p 6 2", "192;d6h^2;p 6/m c c;-p 6 2c", "193;d6h^3;p 63/m c m;-p 6c 2", "193*;d6h^3;p 63/m c m*;-p 63 2", "194;d6h^4;p 63/m m c;-p 6c 2c", "194*;d6h^4;p 63/m m c*;-p 63 2c", "195;t^1;p 2 3;p 2 2 3", "196;t^2;f 2 3;f 2 2 3", "197;t^3;i 2 3;i 2 2 3", "198;t^4;p 21 3;p 2ac 2ab 3", "199;t^5;i 21 3;i 2b 2c 3", "200;th^1;p m -3;-p 2 2 3", "201:1;th^2;p n -3:1;p 2 2 3 -1n;-b", "201:2;th^2;p n -3:2;-p 2ab 2bc 3", "202;th^3;f m -3;-f 2 2 3", "203:1;th^4;f d -3:1;f 2 2 3 -1d;-b", "203:2;th^4;f d -3:2;-f 2uv 2vw 3", "204;th^5;i m -3;-i 2 2 3", "205;th^6;p a -3;-p 2ac 2ab 3", "206;th^7;i a -3;-i 2b 2c 3", "207;o^1;p 4 3 2;p 4 2 3", "208;o^2;p 42 3 2;p 4n 2 3", "209;o^3;f 4 3 2;f 4 2 3", "210;o^4;f 41 3 2;f 4d 2 3", "211;o^5;i 4 3 2;i 4 2 3", "212;o^6;p 43 3 2;p 4acd 2ab 3", "213;o^7;p 41 3 2;p 4bd 2ab 3", "214;o^8;i 41 3 2;i 4bd 2c 3", "215;td^1;p -4 3 m;p -4 2 3", "216;td^2;f -4 3 m;f -4 2 3", "217;td^3;i -4 3 m;i -4 2 3", "218;td^4;p -4 3 n;p -4n 2 3", "219;td^5;f -4 3 c;f -4a 2 3", "220;td^6;i -4 3 d;i -4bd 2c 3", "221;oh^1;p m -3 m;-p 4 2 3", "222:1;oh^2;p n -3 n:1;p 4 2 3 -1n;-b", "222:2;oh^2;p n -3 n:2;-p 4a 2bc 3", "223;oh^3;p m -3 n;-p 4n 2 3", "224:1;oh^4;p n -3 m:1;p 4n 2 3 -1n;-b", "224:2;oh^4;p n -3 m:2;-p 4bc 2bc 3", "225;oh^5;f m -3 m;-f 4 2 3", "226;oh^6;f m -3 c;-f 4a 2 3", "227:1;oh^7;f d -3 m:1;f 4d 2 3 -1d;-b", "227:2;oh^7;f d -3 m:2;-f 4vw 2vw 3", "228:1;oh^8;f d -3 c:1;f 4d 2 3 -1ad;-b", "228:2;oh^8;f d -3 c:2;-f 4ud 2vw 3", "229;oh^9;i m -3 m;-i 4 2 3", "230;oh^10;i a -3 d;-i 4bd 2c 3"]);
});
Clazz_declarePackage ("JS");
Clazz_load (null, "JS.HallInfo", ["JU.P3i", "$.SB", "JS.HallRotationTerm", "$.HallTranslation", "JU.Logger"], function () {
c$ = Clazz_decorateAsClass (function () {
this.hallSymbol = null;
this.primitiveHallSymbol = null;
this.latticeCode = '\0';
this.latticeExtension = null;
this.isCentrosymmetric = false;
this.nRotations = 0;
this.rotationTerms = null;
this.vector12ths = null;
this.vectorCode = null;
Clazz_instantialize (this, arguments);
}, JS, "HallInfo");
Clazz_prepareFields (c$, function () {
this.rotationTerms =  new Array (16);
});
Clazz_makeConstructor (c$, 
function (hallSymbol) {
try {
var str = this.hallSymbol = hallSymbol.trim ();
str = this.extractLatticeInfo (str);
if (JS.HallTranslation.getLatticeIndex (this.latticeCode) == 0) return;
this.latticeExtension = JS.HallTranslation.getLatticeExtension (this.latticeCode, this.isCentrosymmetric);
str = this.extractVectorInfo (str) + this.latticeExtension;
if (JU.Logger.debugging) JU.Logger.debug ("Hallinfo: " + hallSymbol + " " + str);
var prevOrder = 0;
var prevAxisType = '\u0000';
this.primitiveHallSymbol = "P";
while (str.length > 0 && this.nRotations < 16) {
str = this.extractRotationInfo (str, prevOrder, prevAxisType);
var r = this.rotationTerms[this.nRotations - 1];
prevOrder = r.order;
prevAxisType = r.axisType;
this.primitiveHallSymbol += " " + r.primitiveCode;
}
this.primitiveHallSymbol += this.vectorCode;
} catch (e) {
if (Clazz_exceptionOf (e, Exception)) {
JU.Logger.error ("Invalid Hall symbol " + e);
this.nRotations = 0;
} else {
throw e;
}
}
}, "~S");
Clazz_defineMethod (c$, "dumpInfo", 
function () {
var sb =  new JU.SB ();
sb.append ("\nHall symbol: ").append (this.hallSymbol).append ("\nprimitive Hall symbol: ").append (this.primitiveHallSymbol).append ("\nlattice type: ").append (this.getLatticeDesignation ());
for (var i = 0; i < this.nRotations; i++) {
sb.append ("\n\nrotation term ").appendI (i + 1).append (this.rotationTerms[i].dumpInfo (this.vectorCode));
}
return sb.toString ();
});
Clazz_defineMethod (c$, "getLatticeDesignation", 
 function () {
return JS.HallTranslation.getLatticeDesignation2 (this.latticeCode, this.isCentrosymmetric);
});
Clazz_defineMethod (c$, "extractLatticeInfo", 
 function (name) {
var i = name.indexOf (" ");
if (i < 0) return "";
var term = name.substring (0, i).toUpperCase ();
this.latticeCode = term.charAt (0);
if (this.latticeCode == '-') {
this.isCentrosymmetric = true;
this.latticeCode = term.charAt (1);
}return name.substring (i + 1).trim ();
}, "~S");
Clazz_defineMethod (c$, "extractVectorInfo", 
 function (name) {
this.vector12ths =  new JU.P3i ();
this.vectorCode = "";
var i = name.indexOf ("(");
var j = name.indexOf (")", i);
if (i > 0 && j > i) {
var term = name.substring (i + 1, j);
this.vectorCode = " (" + term + ")";
name = name.substring (0, i).trim ();
i = term.indexOf (" ");
if (i >= 0) {
this.vector12ths.x = Integer.parseInt (term.substring (0, i));
term = term.substring (i + 1).trim ();
i = term.indexOf (" ");
if (i >= 0) {
this.vector12ths.y = Integer.parseInt (term.substring (0, i));
term = term.substring (i + 1).trim ();
}}this.vector12ths.z = Integer.parseInt (term);
}return name;
}, "~S");
Clazz_defineMethod (c$, "extractRotationInfo", 
 function (name, prevOrder, prevAxisType) {
var i = name.indexOf (" ");
var code;
if (i >= 0) {
code = name.substring (0, i);
name = name.substring (i + 1).trim ();
} else {
code = name;
name = "";
}this.rotationTerms[this.nRotations] =  new JS.HallRotationTerm (this, code, prevOrder, prevAxisType);
this.nRotations++;
return name;
}, "~S,~N,~S");
});
Clazz_declarePackage ("JS");
Clazz_load (["JU.M4"], "JS.HallRotationTerm", ["JU.SB", "JS.HallRotation", "$.HallTranslation", "$.SymmetryOperation", "JU.Logger"], function () {
c$ = Clazz_decorateAsClass (function () {
this.inputCode = null;
this.primitiveCode = null;
this.lookupCode = null;
this.translationString = null;
this.rotation = null;
this.translation = null;
this.seitzMatrix12ths = null;
this.isImproper = false;
this.order = 0;
this.axisType = '\0';
this.diagonalReferenceAxis = '\0';
this.allPositive = true;
Clazz_instantialize (this, arguments);
}, JS, "HallRotationTerm");
Clazz_prepareFields (c$, function () {
this.seitzMatrix12ths =  new JU.M4 ();
});
Clazz_makeConstructor (c$, 
function (hallInfo, code, prevOrder, prevAxisType) {
this.inputCode = code;
code += "   ";
if (code.charAt (0) == '-') {
this.isImproper = true;
code = code.substring (1);
}this.primitiveCode = "";
this.order = code.charCodeAt (0) - 48;
this.diagonalReferenceAxis = '\0';
this.axisType = '\0';
var ptr = 2;
var c;
switch (c = code.charAt (1)) {
case 'x':
case 'y':
case 'z':
switch (code.charAt (2)) {
case '\'':
case '"':
this.diagonalReferenceAxis = c;
c = code.charAt (2);
ptr++;
}
case '*':
this.axisType = c;
break;
case '\'':
case '"':
this.axisType = c;
switch (code.charAt (2)) {
case 'x':
case 'y':
case 'z':
this.diagonalReferenceAxis = code.charAt (2);
ptr++;
break;
default:
this.diagonalReferenceAxis = prevAxisType;
}
break;
default:
this.axisType = (this.order == 1 ? '_' : hallInfo.nRotations == 0 ? 'z' : hallInfo.nRotations == 2 ? '*' : prevOrder == 2 || prevOrder == 4 ? 'x' : '\'');
code = code.substring (0, 1) + this.axisType + code.substring (1);
}
this.primitiveCode += (this.axisType == '_' ? "1" : code.substring (0, 2));
if (this.diagonalReferenceAxis != '\0') {
code = code.substring (0, 1) + this.diagonalReferenceAxis + this.axisType + code.substring (ptr);
this.primitiveCode += this.diagonalReferenceAxis;
ptr = 3;
}this.lookupCode = code.substring (0, ptr);
this.rotation = JS.HallRotation.lookup (this.lookupCode);
if (this.rotation == null) {
JU.Logger.error ("Rotation lookup could not find " + this.inputCode + " ? " + this.lookupCode);
return;
}this.translation =  new JS.HallTranslation ('\0', null);
this.translationString = "";
var len = code.length;
for (var i = ptr; i < len; i++) {
var translationCode = code.charAt (i);
var t = JS.HallTranslation.getHallTranslation (translationCode, this.order);
if (t != null) {
this.translationString += "" + t.translationCode;
this.translation.rotationShift12ths += t.rotationShift12ths;
this.translation.vectorShift12ths.add (t.vectorShift12ths);
}}
this.primitiveCode = (this.isImproper ? "-" : "") + this.primitiveCode + this.translationString;
if (this.isImproper) {
this.seitzMatrix12ths.setM4 (this.rotation.seitzMatrixInv);
} else {
this.seitzMatrix12ths.setM4 (this.rotation.seitzMatrix);
}this.seitzMatrix12ths.m03 = this.translation.vectorShift12ths.x;
this.seitzMatrix12ths.m13 = this.translation.vectorShift12ths.y;
this.seitzMatrix12ths.m23 = this.translation.vectorShift12ths.z;
switch (this.axisType) {
case 'x':
this.seitzMatrix12ths.m03 += this.translation.rotationShift12ths;
break;
case 'y':
this.seitzMatrix12ths.m13 += this.translation.rotationShift12ths;
break;
case 'z':
this.seitzMatrix12ths.m23 += this.translation.rotationShift12ths;
break;
}
if (hallInfo.vectorCode.length > 0) {
var m1 = JU.M4.newM4 (null);
var m2 = JU.M4.newM4 (null);
var v = hallInfo.vector12ths;
m1.m03 = v.x;
m1.m13 = v.y;
m1.m23 = v.z;
m2.m03 = -v.x;
m2.m13 = -v.y;
m2.m23 = -v.z;
this.seitzMatrix12ths.mul2 (m1, this.seitzMatrix12ths);
this.seitzMatrix12ths.mul (m2);
}if (JU.Logger.debugging) {
JU.Logger.debug ("code = " + code + "; primitive code =" + this.primitiveCode + "\n Seitz Matrix(12ths):" + this.seitzMatrix12ths);
}}, "JS.HallInfo,~S,~N,~S");
Clazz_defineMethod (c$, "dumpInfo", 
function (vectorCode) {
var sb =  new JU.SB ();
sb.append ("\ninput code: ").append (this.inputCode).append ("; primitive code: ").append (this.primitiveCode).append ("\norder: ").appendI (this.order).append (this.isImproper ? " (improper axis)" : "");
if (this.axisType != '_') {
sb.append ("; axisType: ").appendC (this.axisType);
if (this.diagonalReferenceAxis != '\0') sb.appendC (this.diagonalReferenceAxis);
}if (this.translationString.length > 0) sb.append ("; translation: ").append (this.translationString);
if (vectorCode.length > 0) sb.append ("; vector offset: ").append (vectorCode);
if (this.rotation != null) sb.append ("\noperator: ").append (this.getXYZ (this.allPositive)).append ("\nSeitz matrix:\n").append (JS.SymmetryOperation.dumpSeitz (this.seitzMatrix12ths, false));
return sb.toString ();
}, "~S");
Clazz_defineMethod (c$, "getXYZ", 
function (allPositive) {
return JS.SymmetryOperation.getXYZFromMatrix (this.seitzMatrix12ths, true, allPositive, true);
}, "~B");
});
Clazz_declarePackage ("JS");
Clazz_load (["JU.M4"], "JS.HallRotation", null, function () {
c$ = Clazz_decorateAsClass (function () {
this.rotCode = null;
this.seitzMatrix = null;
this.seitzMatrixInv = null;
Clazz_instantialize (this, arguments);
}, JS, "HallRotation");
Clazz_prepareFields (c$, function () {
this.seitzMatrix =  new JU.M4 ();
this.seitzMatrixInv =  new JU.M4 ();
});
Clazz_makeConstructor (c$, 
 function (code, matrixData) {
this.rotCode = code;
var data =  Clazz_newFloatArray (16, 0);
var dataInv =  Clazz_newFloatArray (16, 0);
data[15] = dataInv[15] = 1;
for (var i = 0, ipt = 0; ipt < 11; i++) {
var value = 0;
switch (matrixData.charAt (i)) {
case ' ':
ipt++;
continue;
case '+':
case '1':
value = 1;
break;
case '-':
value = -1;
break;
}
data[ipt] = value;
dataInv[ipt] = -value;
ipt++;
}
this.seitzMatrix.setA (data);
this.seitzMatrixInv.setA (dataInv);
}, "~S,~S");
c$.lookup = Clazz_defineMethod (c$, "lookup", 
function (code) {
for (var i = JS.HallRotation.getHallTerms ().length; --i >= 0; ) if (JS.HallRotation.hallRotationTerms[i].rotCode.equals (code)) return JS.HallRotation.hallRotationTerms[i];

return null;
}, "~S");
c$.getHallTerms = Clazz_defineMethod (c$, "getHallTerms", 
 function () {
return (JS.HallRotation.hallRotationTerms == null ? JS.HallRotation.hallRotationTerms = [ new JS.HallRotation ("1_", "+00 0+0 00+"),  new JS.HallRotation ("2x", "+00 0-0 00-"),  new JS.HallRotation ("2y", "-00 0+0 00-"),  new JS.HallRotation ("2z", "-00 0-0 00+"),  new JS.HallRotation ("2\'", "0-0 -00 00-"),  new JS.HallRotation ("2\"", "0+0 +00 00-"),  new JS.HallRotation ("2x\'", "-00 00- 0-0"),  new JS.HallRotation ("2x\"", "-00 00+ 0+0"),  new JS.HallRotation ("2y\'", "00- 0-0 -00"),  new JS.HallRotation ("2y\"", "00+ 0-0 +00"),  new JS.HallRotation ("2z\'", "0-0 -00 00-"),  new JS.HallRotation ("2z\"", "0+0 +00 00-"),  new JS.HallRotation ("3x", "+00 00- 0+-"),  new JS.HallRotation ("3y", "-0+ 0+0 -00"),  new JS.HallRotation ("3z", "0-0 +-0 00+"),  new JS.HallRotation ("3*", "00+ +00 0+0"),  new JS.HallRotation ("4x", "+00 00- 0+0"),  new JS.HallRotation ("4y", "00+ 0+0 -00"),  new JS.HallRotation ("4z", "0-0 +00 00+"),  new JS.HallRotation ("6x", "+00 0+- 0+0"),  new JS.HallRotation ("6y", "00+ 0+0 -0+"),  new JS.HallRotation ("6z", "+-0 +00 00+")] : JS.HallRotation.hallRotationTerms);
});
Clazz_defineStatics (c$,
"hallRotationTerms", null);
});
Clazz_declarePackage ("JS");
Clazz_load (null, "JS.HallTranslation", ["JU.P3i"], function () {
c$ = Clazz_decorateAsClass (function () {
this.translationCode = '\0';
this.rotationOrder = 0;
this.rotationShift12ths = 0;
this.vectorShift12ths = null;
Clazz_instantialize (this, arguments);
}, JS, "HallTranslation");
Clazz_makeConstructor (c$, 
function (translationCode, params) {
this.translationCode = translationCode;
if (params != null) {
if (params.z >= 0) {
this.vectorShift12ths = params;
return;
}this.rotationOrder = params.x;
this.rotationShift12ths = params.y;
}this.vectorShift12ths =  new JU.P3i ();
}, "~S,JU.P3i");
c$.getHallLatticeEquivalent = Clazz_defineMethod (c$, "getHallLatticeEquivalent", 
function (latticeParameter) {
var latticeCode = JS.HallTranslation.getLatticeCode (latticeParameter);
var isCentrosymmetric = (latticeParameter > 0);
return (isCentrosymmetric ? "-" : "") + latticeCode + " 1";
}, "~N");
c$.getLatticeIndex = Clazz_defineMethod (c$, "getLatticeIndex", 
function (latt) {
for (var i = 1, ipt = 3; i <= JS.HallTranslation.nLatticeTypes; i++, ipt += 3) if (JS.HallTranslation.latticeTranslationData[ipt].charAt (0) == latt) return i;

return 0;
}, "~S");
c$.getLatticeCode = Clazz_defineMethod (c$, "getLatticeCode", 
function (latt) {
if (latt < 0) latt = -latt;
return (latt == 0 ? '\0' : latt > JS.HallTranslation.nLatticeTypes ? JS.HallTranslation.getLatticeCode (JS.HallTranslation.getLatticeIndex (String.fromCharCode (latt))) : JS.HallTranslation.latticeTranslationData[latt * 3].charAt (0));
}, "~N");
c$.getLatticeDesignation = Clazz_defineMethod (c$, "getLatticeDesignation", 
function (latt) {
var isCentrosymmetric = (latt > 0);
var str = (isCentrosymmetric ? "-" : "");
if (latt < 0) latt = -latt;
if (latt == 0 || latt > JS.HallTranslation.nLatticeTypes) return "";
return str + JS.HallTranslation.getLatticeCode (latt) + ": " + (isCentrosymmetric ? "centrosymmetric " : "") + JS.HallTranslation.latticeTranslationData[latt * 3 + 1];
}, "~N");
c$.getLatticeDesignation2 = Clazz_defineMethod (c$, "getLatticeDesignation2", 
function (latticeCode, isCentrosymmetric) {
var latt = JS.HallTranslation.getLatticeIndex (latticeCode);
if (!isCentrosymmetric) latt = -latt;
return JS.HallTranslation.getLatticeDesignation (latt);
}, "~S,~B");
c$.getLatticeExtension = Clazz_defineMethod (c$, "getLatticeExtension", 
function (latt, isCentrosymmetric) {
for (var i = 1, ipt = 3; i <= JS.HallTranslation.nLatticeTypes; i++, ipt += 3) if (JS.HallTranslation.latticeTranslationData[ipt].charAt (0) == latt) return JS.HallTranslation.latticeTranslationData[ipt + 2] + (isCentrosymmetric ? " -1" : "");

return "";
}, "~S,~B");
c$.getHallTerms = Clazz_defineMethod (c$, "getHallTerms", 
 function () {
return (JS.HallTranslation.hallTranslationTerms == null ? JS.HallTranslation.hallTranslationTerms = [ new JS.HallTranslation ('a', JU.P3i.new3 (6, 0, 0)),  new JS.HallTranslation ('b', JU.P3i.new3 (0, 6, 0)),  new JS.HallTranslation ('c', JU.P3i.new3 (0, 0, 6)),  new JS.HallTranslation ('n', JU.P3i.new3 (6, 6, 6)),  new JS.HallTranslation ('u', JU.P3i.new3 (3, 0, 0)),  new JS.HallTranslation ('v', JU.P3i.new3 (0, 3, 0)),  new JS.HallTranslation ('w', JU.P3i.new3 (0, 0, 3)),  new JS.HallTranslation ('d', JU.P3i.new3 (3, 3, 3)),  new JS.HallTranslation ('1', JU.P3i.new3 (2, 6, -1)),  new JS.HallTranslation ('1', JU.P3i.new3 (3, 4, -1)),  new JS.HallTranslation ('2', JU.P3i.new3 (3, 8, -1)),  new JS.HallTranslation ('1', JU.P3i.new3 (4, 3, -1)),  new JS.HallTranslation ('3', JU.P3i.new3 (4, 9, -1)),  new JS.HallTranslation ('1', JU.P3i.new3 (6, 2, -1)),  new JS.HallTranslation ('2', JU.P3i.new3 (6, 4, -1)),  new JS.HallTranslation ('4', JU.P3i.new3 (6, 8, -1)),  new JS.HallTranslation ('5', JU.P3i.new3 (6, 10, -1)),  new JS.HallTranslation ('r', JU.P3i.new3 (4, 8, 8)),  new JS.HallTranslation ('s', JU.P3i.new3 (8, 8, 4)),  new JS.HallTranslation ('t', JU.P3i.new3 (8, 4, 8))] : JS.HallTranslation.hallTranslationTerms);
});
c$.getHallTranslation = Clazz_defineMethod (c$, "getHallTranslation", 
function (translationCode, order) {
var ht = null;
for (var i = JS.HallTranslation.getHallTerms ().length; --i >= 0; ) {
var h = JS.HallTranslation.hallTranslationTerms[i];
if (h.translationCode == translationCode) {
if (h.rotationOrder == 0 || h.rotationOrder == order) {
ht =  new JS.HallTranslation (translationCode, null);
ht.translationCode = translationCode;
ht.rotationShift12ths = h.rotationShift12ths;
ht.vectorShift12ths = h.vectorShift12ths;
return ht;
}}}
return ht;
}, "~S,~N");
Clazz_defineStatics (c$,
"latticeTranslationData", ["\0", "unknown", "", "P", "primitive", "", "I", "body-centered", " 1n", "R", "rhombohedral", " 1r 1r", "F", "face-centered", " 1ab 1bc 1ac", "A", "A-centered", " 1bc", "B", "B-centered", " 1ac", "C", "C-centered", " 1ab", "S", "rhombohedral(S)", " 1s 1s", "T", "rhombohedral(T)", " 1t 1t"]);
c$.nLatticeTypes = c$.prototype.nLatticeTypes = Clazz_doubleToInt (JS.HallTranslation.latticeTranslationData.length / 3) - 1;
Clazz_defineStatics (c$,
"hallTranslationTerms", null);
});
Clazz_declarePackage ("JS");
Clazz_load (["JU.M4"], "JS.SymmetryOperation", ["java.lang.Float", "JU.Lst", "$.Matrix", "$.P3", "$.P4", "$.PT", "$.Quat", "$.SB", "$.V3", "JU.Escape", "$.Logger", "$.Measure", "$.Parser"], function () {
c$ = Clazz_decorateAsClass (function () {
this.xyzOriginal = null;
this.xyz = null;
this.doNormalize = true;
this.isFinalized = false;
this.opId = 0;
this.centering = null;
this.atomTest = null;
this.myLabels = null;
this.modDim = 0;
this.linearRotTrans = null;
this.rsvs = null;
this.isBio = false;
this.sigma = null;
this.index = 0;
this.subsystemCode = null;
this.timeReversal = 0;
this.magOp = 3.4028235E38;
this.isCenteringOp = false;
this.unCentered = false;
Clazz_instantialize (this, arguments);
}, JS, "SymmetryOperation", JU.M4);
Clazz_defineMethod (c$, "setSigma", 
function (subsystemCode, sigma) {
this.subsystemCode = subsystemCode;
this.sigma = sigma;
}, "~S,JU.Matrix");
Clazz_overrideConstructor (c$, 
function (op, atoms, atomIndex, countOrId, doNormalize) {
this.doNormalize = doNormalize;
if (op == null) {
this.opId = countOrId;
return;
}this.xyzOriginal = op.xyzOriginal;
this.xyz = op.xyz;
this.opId = op.opId;
this.modDim = op.modDim;
this.myLabels = op.myLabels;
this.index = op.index;
this.linearRotTrans = op.linearRotTrans;
this.sigma = op.sigma;
this.subsystemCode = op.subsystemCode;
this.timeReversal = op.timeReversal;
this.setMatrix (false);
if (!op.isFinalized) this.doFinalize ();
if (doNormalize && this.sigma == null) this.setOffset (atoms, atomIndex, countOrId);
}, "JS.SymmetryOperation,~A,~N,~N,~B");
Clazz_defineMethod (c$, "setGamma", 
 function (isReverse) {
var n = 3 + this.modDim;
var a = (this.rsvs =  new JU.Matrix (null, n + 1, n + 1)).getArray ();
var t =  Clazz_newDoubleArray (n, 0);
var pt = 0;
for (var i = 0; i < n; i++) {
for (var j = 0; j < n; j++) a[i][j] = this.linearRotTrans[pt++];

t[i] = (isReverse ? -1 : 1) * this.linearRotTrans[pt++];
}
a[n][n] = 1;
if (isReverse) this.rsvs = this.rsvs.inverse ();
for (var i = 0; i < n; i++) a[i][n] = t[i];

a = this.rsvs.getSubmatrix (0, 0, 3, 3).getArray ();
for (var i = 0; i < 3; i++) for (var j = 0; j < 4; j++) this.setElement (i, j, (j < 3 ? a[i][j] : t[i]));


this.setElement (3, 3, 1);
}, "~B");
Clazz_defineMethod (c$, "doFinalize", 
function () {
this.m03 /= 12;
this.m13 /= 12;
this.m23 /= 12;
if (this.modDim > 0) {
var a = this.rsvs.getArray ();
for (var i = a.length - 1; --i >= 0; ) a[i][3 + this.modDim] /= 12;

}this.isFinalized = true;
});
Clazz_defineMethod (c$, "getXyz", 
function (normalized) {
return (normalized && this.modDim == 0 || this.xyzOriginal == null ? this.xyz : this.xyzOriginal);
}, "~B");
Clazz_defineMethod (c$, "newPoint", 
function (atom1, atom2, x, y, z) {
this.rotTrans2 (atom1, atom2);
atom2.add3 (x, y, z);
}, "JU.P3,JU.P3,~N,~N,~N");
Clazz_defineMethod (c$, "dumpInfo", 
function () {
return "\n" + this.xyz + "\ninternal matrix representation:\n" + this.toString ();
});
c$.dumpSeitz = Clazz_defineMethod (c$, "dumpSeitz", 
function (s, isCanonical) {
var sb =  new JU.SB ();
var r =  Clazz_newFloatArray (4, 0);
for (var i = 0; i < 3; i++) {
s.getRow (i, r);
sb.append ("[\t");
for (var j = 0; j < 3; j++) sb.appendI (Clazz_floatToInt (r[j])).append ("\t");

sb.append (JS.SymmetryOperation.twelfthsOf (isCanonical ? (Clazz_floatToInt (r[3]) + 12) % 12 : Clazz_floatToInt (r[3]))).append ("\t]\n");
}
return sb.toString ();
}, "JU.M4,~B");
Clazz_defineMethod (c$, "setMatrixFromXYZ", 
function (xyz, modDim, allowScaling) {
if (xyz == null) return false;
this.xyzOriginal = xyz;
xyz = xyz.toLowerCase ();
var n = (modDim + 4) * (modDim + 4);
this.modDim = modDim;
if (modDim > 0) this.myLabels = JS.SymmetryOperation.labelsXn;
this.linearRotTrans =  Clazz_newFloatArray (n, 0);
var isReverse = (xyz.startsWith ("!"));
if (isReverse) xyz = xyz.substring (1);
if (xyz.indexOf ("xyz matrix:") == 0) {
this.xyz = xyz;
JU.Parser.parseStringInfestedFloatArray (xyz, null, this.linearRotTrans);
return this.setFromMatrix (null, isReverse);
}if (xyz.indexOf ("[[") == 0) {
xyz = xyz.$replace ('[', ' ').$replace (']', ' ').$replace (',', ' ');
JU.Parser.parseStringInfestedFloatArray (xyz, null, this.linearRotTrans);
for (var i = 0; i < n; i++) {
var v = this.linearRotTrans[i];
if (Float.isNaN (v)) return false;
}
this.setMatrix (isReverse);
this.isFinalized = true;
this.isBio = (xyz.indexOf ("bio") >= 0);
this.xyz = (this.isBio ? this.toString () : JS.SymmetryOperation.getXYZFromMatrix (this, false, false, false));
return true;
}if (xyz.indexOf (",m") >= 0) {
this.timeReversal = (xyz.indexOf (",m-1") >= 0 ? -1 : 1);
allowScaling = true;
}var strOut = JS.SymmetryOperation.getMatrixFromString (this, xyz, this.linearRotTrans, allowScaling);
if (strOut == null) return false;
this.setMatrix (isReverse);
this.xyz = (isReverse ? JS.SymmetryOperation.getXYZFromMatrix (this, true, false, false) : strOut);
if (this.timeReversal != 0) this.xyz += (this.timeReversal == 1 ? ",m+1" : ",m-1");
if (JU.Logger.debugging) JU.Logger.debug ("" + this);
return true;
}, "~S,~N,~B");
Clazz_defineMethod (c$, "setMatrix", 
 function (isReverse) {
if (this.linearRotTrans.length > 16) {
this.setGamma (isReverse);
} else {
this.setA (this.linearRotTrans);
if (isReverse) {
var p3 = JU.P3.new3 (this.m03, this.m13, this.m23);
this.invertM (this);
this.rotate (p3);
p3.scale (-1);
this.setTranslation (p3);
}}}, "~B");
Clazz_defineMethod (c$, "setFromMatrix", 
function (offset, isReverse) {
var v = 0;
var pt = 0;
this.myLabels = (this.modDim == 0 ? JS.SymmetryOperation.labelsXYZ : JS.SymmetryOperation.labelsXn);
var rowPt = 0;
var n = 3 + this.modDim;
for (var i = 0; rowPt < n; i++) {
if (Float.isNaN (this.linearRotTrans[i])) return false;
v = this.linearRotTrans[i];
if (Math.abs (v) < 0.00001) v = 0;
var isTrans = ((i + 1) % (n + 1) == 0);
if (isTrans) {
if (offset != null) {
v /= 12;
if (pt < offset.length) v += offset[pt++];
}v = JS.SymmetryOperation.normalizeTwelfths ((v < 0 ? -1 : 1) * Math.round (Math.abs (v * 12)) / 12, this.doNormalize);
rowPt++;
}this.linearRotTrans[i] = v;
}
this.linearRotTrans[this.linearRotTrans.length - 1] = 1;
this.setMatrix (isReverse);
this.isFinalized = (offset == null);
this.xyz = JS.SymmetryOperation.getXYZFromMatrix (this, true, false, false);
return true;
}, "~A,~B");
c$.getMatrixFromString = Clazz_defineMethod (c$, "getMatrixFromString", 
function (op, xyz, linearRotTrans, allowScaling) {
var isDenominator = false;
var isDecimal = false;
var isNegative = false;
var modDim = (op == null ? 0 : op.modDim);
var nRows = 4 + modDim;
var doNormalize = (op != null && op.doNormalize);
var dimOffset = (modDim > 0 ? 3 : 0);
linearRotTrans[linearRotTrans.length - 1] = 1;
var transPt = xyz.indexOf (';') + 1;
if (transPt != 0) {
allowScaling = true;
if (transPt == xyz.length) xyz += "0,0,0";
}var rotPt = -1;
var myLabels = (op == null || modDim == 0 ? null : op.myLabels);
if (myLabels == null) myLabels = JS.SymmetryOperation.labelsXYZ;
xyz = xyz.toLowerCase () + ",";
if (modDim > 0) for (var i = modDim + 3; --i >= 0; ) xyz = JU.PT.rep (xyz, JS.SymmetryOperation.labelsXn[i], JS.SymmetryOperation.labelsXnSub[i]);

var xpt = 0;
var tpt0 = 0;
var rowPt = 0;
var ch;
var iValue = 0;
var decimalMultiplier = 1;
var strT = "";
var strOut = "";
for (var i = 0; i < xyz.length; i++) {
switch (ch = xyz.charAt (i)) {
case ';':
break;
case '\'':
case ' ':
case '{':
case '}':
case '!':
continue;
case '-':
isNegative = true;
continue;
case '+':
isNegative = false;
continue;
case '/':
isDenominator = true;
continue;
case 'x':
case 'y':
case 'z':
case 'a':
case 'b':
case 'c':
case 'd':
case 'e':
case 'f':
case 'g':
case 'h':
tpt0 = rowPt * nRows;
var ipt = (ch >= 'x' ? ch.charCodeAt (0) - 120 : ch.charCodeAt (0) - 97 + dimOffset);
xpt = tpt0 + ipt;
var val = (isNegative ? -1 : 1);
if (allowScaling && iValue != 0) {
linearRotTrans[xpt] = iValue;
val = Clazz_floatToInt (iValue);
iValue = 0;
} else {
linearRotTrans[xpt] = val;
}strT += JS.SymmetryOperation.plusMinus (strT, val, myLabels[ipt]);
break;
case ',':
if (transPt != 0) {
if (transPt > 0) {
rotPt = i;
i = transPt - 1;
transPt = -i;
iValue = 0;
continue;
}transPt = i + 1;
i = rotPt;
}iValue = JS.SymmetryOperation.normalizeTwelfths (iValue, doNormalize);
linearRotTrans[tpt0 + nRows - 1] = iValue;
strT += JS.SymmetryOperation.xyzFraction (iValue, false, true);
strOut += (strOut === "" ? "" : ",") + strT;
if (rowPt == nRows - 2) return strOut;
iValue = 0;
strT = "";
if (rowPt++ > 2 && modDim == 0) {
JU.Logger.warn ("Symmetry Operation? " + xyz);
return null;
}break;
case '.':
isDecimal = true;
decimalMultiplier = 1;
continue;
case '0':
if (!isDecimal && (isDenominator || !allowScaling)) continue;
default:
var ich = ch.charCodeAt (0) - 48;
if (isDecimal && ich >= 0 && ich <= 9) {
decimalMultiplier /= 10;
if (iValue < 0) isNegative = true;
iValue += decimalMultiplier * ich * (isNegative ? -1 : 1);
continue;
}if (ich >= 0 && ich <= 9) {
if (isDenominator) {
if (iValue == 0) {
linearRotTrans[xpt] /= ich;
} else {
iValue /= ich;
}} else {
iValue = iValue * 10 + (isNegative ? -1 : 1) * ich;
isNegative = false;
}} else {
JU.Logger.warn ("symmetry character?" + ch);
}}
isDecimal = isDenominator = isNegative = false;
}
return null;
}, "JS.SymmetryOperation,~S,~A,~B");
c$.xyzFraction = Clazz_defineMethod (c$, "xyzFraction", 
 function (n12ths, allPositive, halfOrLess) {
n12ths = Math.round (n12ths);
if (allPositive) {
while (n12ths < 0) n12ths += 12;

} else if (halfOrLess) {
while (n12ths > 6) n12ths -= 12;

while (n12ths < -6.0) n12ths += 12;

}var s = JS.SymmetryOperation.twelfthsOf (n12ths);
return (s.charAt (0) == '0' ? "" : n12ths > 0 ? "+" + s : s);
}, "~N,~B,~B");
c$.twelfthsOf = Clazz_defineMethod (c$, "twelfthsOf", 
 function (n12ths) {
var str = "";
var i12ths = Math.round (n12ths);
if (i12ths == 12) return "1";
if (i12ths == -12) return "-1";
if (i12ths < 0) {
i12ths = -i12ths;
if (i12ths % 12 != 0) str = "-";
}var n = Clazz_doubleToInt (i12ths / 12);
if (n < 1) return str + JS.SymmetryOperation.twelfths[i12ths % 12];
var m = 0;
switch (i12ths % 12) {
case 0:
return str + n;
case 1:
case 5:
case 7:
case 11:
m = 12;
break;
case 2:
case 10:
m = 6;
break;
case 3:
case 9:
m = 4;
break;
case 4:
case 8:
m = 3;
break;
case 6:
m = 2;
break;
}
return str + (Clazz_doubleToInt (i12ths * m / 12)) + "/" + m;
}, "~N");
c$.plusMinus = Clazz_defineMethod (c$, "plusMinus", 
 function (strT, x, sx) {
return (x == 0 ? "" : (x < 0 ? "-" : strT.length == 0 ? "" : "+") + (x == 1 || x == -1 ? "" : "" + Clazz_floatToInt (Math.abs (x))) + sx);
}, "~S,~N,~S");
c$.normalizeTwelfths = Clazz_defineMethod (c$, "normalizeTwelfths", 
 function (iValue, doNormalize) {
iValue *= 12;
if (doNormalize) {
while (iValue > 6) iValue -= 12;

while (iValue <= -6) iValue += 12;

}return iValue;
}, "~N,~B");
c$.getXYZFromMatrix = Clazz_defineMethod (c$, "getXYZFromMatrix", 
function (mat, is12ths, allPositive, halfOrLess) {
var str = "";
var op = (Clazz_instanceOf (mat, JS.SymmetryOperation) ? mat : null);
if (op != null && op.modDim > 0) return JS.SymmetryOperation.getXYZFromRsVs (op.rsvs.getRotation (), op.rsvs.getTranslation (), is12ths);
var row =  Clazz_newFloatArray (4, 0);
for (var i = 0; i < 3; i++) {
var lpt = (i < 3 ? 0 : 3);
mat.getRow (i, row);
var term = "";
for (var j = 0; j < 3; j++) if (row[j] != 0) term += JS.SymmetryOperation.plusMinus (term, row[j], JS.SymmetryOperation.labelsXYZ[j + lpt]);

term += JS.SymmetryOperation.xyzFraction ((is12ths ? row[3] : row[3] * 12), allPositive, halfOrLess);
str += "," + term;
}
return str.substring (1);
}, "JU.M4,~B,~B,~B");
Clazz_defineMethod (c$, "setOffset", 
 function (atoms, atomIndex, count) {
var i1 = atomIndex;
var i2 = i1 + count;
var x = 0;
var y = 0;
var z = 0;
if (this.atomTest == null) this.atomTest =  new JU.P3 ();
for (var i = i1; i < i2; i++) {
this.newPoint (atoms[i], this.atomTest, 0, 0, 0);
x += this.atomTest.x;
y += this.atomTest.y;
z += this.atomTest.z;
}
while (x < -0.001 || x >= count + 0.001) {
this.m03 += (x < 0 ? 1 : -1);
x += (x < 0 ? count : -count);
}
while (y < -0.001 || y >= count + 0.001) {
this.m13 += (y < 0 ? 1 : -1);
y += (y < 0 ? count : -count);
}
while (z < -0.001 || z >= count + 0.001) {
this.m23 += (z < 0 ? 1 : -1);
z += (z < 0 ? count : -count);
}
}, "~A,~N,~N");
Clazz_defineMethod (c$, "rotateAxes", 
function (vectors, unitcell, ptTemp, mTemp) {
var vRot =  new Array (3);
this.getRotationScale (mTemp);
for (var i = vectors.length; --i >= 0; ) {
ptTemp.setT (vectors[i]);
unitcell.toFractional (ptTemp, true);
mTemp.rotate (ptTemp);
unitcell.toCartesian (ptTemp, true);
vRot[i] = JU.V3.newV (ptTemp);
}
return vRot;
}, "~A,JS.UnitCell,JU.P3,JU.M3");
Clazz_defineMethod (c$, "getDescription", 
function (modelSet, uc, pta00, ptTarget, id) {
if (!this.isFinalized) this.doFinalize ();
var vtemp =  new JU.V3 ();
var ptemp =  new JU.P3 ();
var ptemp2 =  new JU.P3 ();
var pta01 =  new JU.P3 ();
var pta02 =  new JU.P3 ();
var ftrans =  new JU.V3 ();
var vtrans =  new JU.V3 ();
var vcentering = null;
if (this.centering != null) {
vcentering = JU.V3.newV (this.centering);
uc.toCartesian (vcentering, false);
}var haveCentering = false;
if (pta00 == null || Float.isNaN (pta00.x)) pta00 =  new JU.P3 ();
if (ptTarget != null) {
JS.SymmetryOperation.setFractional (uc, pta01, pta00, ptemp);
this.rotTrans (pta01);
uc.toCartesian (pta01, false);
uc.toUnitCell (pta01, ptemp);
pta02.setT (ptTarget);
uc.toUnitCell (pta02, ptemp);
if (pta01.distance (pta02) > 0.1) return null;
JS.SymmetryOperation.setFractional (uc, pta01, pta00, null);
this.rotTrans (pta01);
JS.SymmetryOperation.setFractional (uc, pta02, ptTarget, null);
vtrans.sub2 (pta02, pta01);
}pta01.set (1, 0, 0);
pta02.set (0, 1, 0);
var pta03 = JU.P3.new3 (0, 0, 1);
pta01.add (pta00);
pta02.add (pta00);
pta03.add (pta00);
if (haveCentering) vtrans.sub (this.centering);
var pt0 = this.rotTransCart (pta00, uc, vtrans);
var pt1 = this.rotTransCart (pta01, uc, vtrans);
var pt2 = this.rotTransCart (pta02, uc, vtrans);
var pt3 = this.rotTransCart (pta03, uc, vtrans);
var vt1 = JU.V3.newVsub (pt1, pt0);
var vt2 = JU.V3.newVsub (pt2, pt0);
var vt3 = JU.V3.newVsub (pt3, pt0);
JS.SymmetryOperation.approx (vtrans);
vtemp.cross (vt1, vt2);
var haveInversion = (vtemp.dot (vt3) < 0);
if (haveInversion) {
pt1.sub2 (pt0, vt1);
pt2.sub2 (pt0, vt2);
pt3.sub2 (pt0, vt3);
}var info;
info = JU.Measure.computeHelicalAxis (null, 135266306, pta00, pt0, JU.Quat.getQuaternionFrame (pt0, pt1, pt2).div (JU.Quat.getQuaternionFrame (pta00, pta01, pta02)));
var pa1 = info[0];
var ax1 = info[1];
var ang1 = Clazz_floatToInt (Math.abs (JU.PT.approx ((info[3]).x, 1)));
var pitch1 = JS.SymmetryOperation.approxF ((info[3]).y);
if (haveInversion) {
pt1.add2 (pt0, vt1);
pt2.add2 (pt0, vt2);
pt3.add2 (pt0, vt3);
}var trans = JU.V3.newVsub (pt0, pta00);
if (trans.length () < 0.1) trans = null;
var ptinv = null;
var ipt = null;
var ptref = null;
var isTranslation = (ang1 == 0);
var isRotation = !isTranslation;
var isInversionOnly = false;
var isMirrorPlane = false;
if (isRotation || haveInversion) trans = null;
if (haveInversion && isTranslation) {
ipt = JU.P3.newP (pta00);
ipt.add (pt0);
ipt.scale (0.5);
ptinv = pt0;
isInversionOnly = true;
} else if (haveInversion) {
var d = (pitch1 == 0 ?  new JU.V3 () : ax1);
var f = 0;
switch (ang1) {
case 60:
f = 0.6666667;
break;
case 120:
f = 2;
break;
case 90:
f = 1;
break;
case 180:
ptref = JU.P3.newP (pta00);
ptref.add (d);
pa1.scaleAdd2 (0.5, d, pta00);
if (ptref.distance (pt0) > 0.1) {
trans = JU.V3.newVsub (pt0, ptref);
JS.SymmetryOperation.setFractional (uc, ptemp, trans, null);
ftrans.setT (ptemp);
} else {
trans = null;
}isRotation = false;
haveInversion = false;
isMirrorPlane = true;
}
if (f != 0) {
vtemp.sub2 (pta00, pa1);
vtemp.add (pt0);
vtemp.sub (pa1);
vtemp.sub (d);
vtemp.scale (f);
pa1.add (vtemp);
ipt =  new JU.P3 ();
ipt.scaleAdd2 (0.5, d, pa1);
ptinv =  new JU.P3 ();
ptinv.scaleAdd2 (-2, ipt, pta00);
ptinv.scale (-1);
}} else if (trans != null) {
ptemp.setT (trans);
uc.toFractional (ptemp, false);
if (JS.SymmetryOperation.approxF (ptemp.x) == 1) {
ptemp.x = 0;
}if (JS.SymmetryOperation.approxF (ptemp.y) == 1) {
ptemp.y = 0;
}if (JS.SymmetryOperation.approxF (ptemp.z) == 1) {
ptemp.z = 0;
}ftrans.setT (ptemp);
uc.toCartesian (ptemp, false);
trans.setT (ptemp);
}var ang = ang1;
JS.SymmetryOperation.approx0 (ax1);
if (isRotation) {
var ptr =  new JU.P3 ();
vtemp.setT (ax1);
var ang2 = ang1;
if (haveInversion) {
ptr.add2 (pa1, vtemp);
ang2 = Math.round (JU.Measure.computeTorsion (ptinv, pa1, ptr, pt0, true));
} else if (pitch1 == 0) {
ptr.setT (pa1);
ptemp.scaleAdd2 (1, ptr, vtemp);
ang2 = Math.round (JU.Measure.computeTorsion (pta00, pa1, ptemp, pt0, true));
} else {
ptemp.add2 (pa1, vtemp);
ptr.scaleAdd2 (0.5, vtemp, pa1);
ang2 = Math.round (JU.Measure.computeTorsion (pta00, pa1, ptemp, pt0, true));
}if (ang2 != 0) ang1 = ang2;
}if (isRotation && !haveInversion && pitch1 == 0) {
if (ax1.z < 0 || ax1.z == 0 && (ax1.y < 0 || ax1.y == 0 && ax1.x < 0)) {
ax1.scale (-1);
ang1 = -ang1;
}}var info1 = "identity";
if (isInversionOnly) {
ptemp.setT (ipt);
uc.toFractional (ptemp, false);
info1 = "Ci: " + this.strCoord (ptemp);
} else if (isRotation) {
if (haveInversion) {
info1 = "" + (Clazz_doubleToInt (360 / ang)) + "-bar axis";
} else if (pitch1 != 0) {
info1 = "" + (Clazz_doubleToInt (360 / ang)) + "-fold screw axis";
ptemp.setT (ax1);
uc.toFractional (ptemp, false);
info1 += "|translation: " + this.strCoord (ptemp);
} else {
info1 = "C" + (Clazz_doubleToInt (360 / ang)) + " axis";
}} else if (trans != null) {
var s = " " + this.strCoord (ftrans);
if (isTranslation) {
info1 = "translation:" + s;
} else if (isMirrorPlane) {
var fx = JS.SymmetryOperation.approxF (ftrans.x);
var fy = JS.SymmetryOperation.approxF (ftrans.y);
var fz = JS.SymmetryOperation.approxF (ftrans.z);
s = " " + this.strCoord (ftrans);
if (fx != 0 && fy != 0 && fz != 0) {
if (Math.abs (fx) == Math.abs (fy) && Math.abs (fy) == Math.abs (fz)) info1 = "d-";
 else info1 = "g-";
} else if (fx != 0 && fy != 0 || fy != 0 && fz != 0 || fz != 0 && fx != 0) info1 = "n-";
 else if (fx != 0) info1 = "a-";
 else if (fy != 0) info1 = "b-";
 else info1 = "c-";
info1 += "glide plane|translation:" + s;
}} else if (isMirrorPlane) {
info1 = "mirror plane";
}if (haveInversion && !isInversionOnly) {
ptemp.setT (ipt);
uc.toFractional (ptemp, false);
info1 += "|inversion center at " + this.strCoord (ptemp);
}if (haveCentering) info1 += "|centering " + this.strCoord (this.centering);
if (this.timeReversal != 0) info1 += "|spin timeReversal " + (this.timeReversal == 1 ? "+1" : "-1");
var cmds = null;
var xyzNew = (this.isBio ? this.xyzOriginal : JS.SymmetryOperation.getXYZFromMatrix (this, false, false, false));
var draw1 =  new JU.SB ();
if (id != null) {
var drawid;
var opType = null;
drawid = "\ndraw ID " + id + "_";
draw1 =  new JU.SB ();
draw1.append (("// " + this.xyzOriginal + "|" + xyzNew + "|" + info1).$replace ('\n', ' ')).append ("\n").append (drawid).append ("* delete");
JS.SymmetryOperation.drawLine (draw1, drawid + "frame1X", 0.15, pta00, pta01, "red");
JS.SymmetryOperation.drawLine (draw1, drawid + "frame1Y", 0.15, pta00, pta02, "green");
JS.SymmetryOperation.drawLine (draw1, drawid + "frame1Z", 0.15, pta00, pta03, "blue");
var color;
if (isRotation) {
var ptr =  new JU.P3 ();
color = "red";
ang = ang1;
var scale = 1.0;
vtemp.setT (ax1);
if (haveInversion) {
opType = drawid + "rotinv";
ptr.add2 (pa1, vtemp);
if (pitch1 == 0) {
ptr.setT (ipt);
vtemp.scale (3);
ptemp.scaleAdd2 (-1, vtemp, pa1);
draw1.append (drawid).append ("rotVector2 diameter 0.1 ").append (JU.Escape.eP (pa1)).append (JU.Escape.eP (ptemp)).append (" color red");
}scale = pt0.distance (ptr);
draw1.append (drawid).append ("rotLine1 ").append (JU.Escape.eP (ptr)).append (JU.Escape.eP (ptinv)).append (" color red");
draw1.append (drawid).append ("rotLine2 ").append (JU.Escape.eP (ptr)).append (JU.Escape.eP (pt0)).append (" color red");
} else if (pitch1 == 0) {
opType = drawid + "rot";
var isSpecial = (pta00.distance (pt0) < 0.2);
if (!isSpecial) {
draw1.append (drawid).append ("rotLine1 ").append (JU.Escape.eP (pta00)).append (JU.Escape.eP (pa1)).append (" color red");
draw1.append (drawid).append ("rotLine2 ").append (JU.Escape.eP (pt0)).append (JU.Escape.eP (pa1)).append (" color red");
}vtemp.scale (3);
ptemp.scaleAdd2 (-1, vtemp, pa1);
draw1.append (drawid).append ("rotVector2 diameter 0.1 ").append (JU.Escape.eP (pa1)).append (JU.Escape.eP (ptemp)).append (" color red");
ptr.setT (pa1);
if (pitch1 == 0 && pta00.distance (pt0) < 0.2) ptr.scaleAdd2 (0.5, ptr, vtemp);
} else {
opType = drawid + "screw";
color = "orange";
draw1.append (drawid).append ("rotLine1 ").append (JU.Escape.eP (pta00)).append (JU.Escape.eP (pa1)).append (" color red");
ptemp.add2 (pa1, vtemp);
draw1.append (drawid).append ("rotLine2 ").append (JU.Escape.eP (pt0)).append (JU.Escape.eP (ptemp)).append (" color red");
ptr.scaleAdd2 (0.5, vtemp, pa1);
}ptemp.add2 (ptr, vtemp);
if (haveInversion && pitch1 != 0) {
draw1.append (drawid).append ("rotRotLine1").append (JU.Escape.eP (ptr)).append (JU.Escape.eP (ptinv)).append (" color red");
draw1.append (drawid).append ("rotRotLine2").append (JU.Escape.eP (ptr)).append (JU.Escape.eP (pt0)).append (" color red");
}draw1.append (drawid).append ("rotRotArrow arrow width 0.10 scale " + scale + " arc ").append (JU.Escape.eP (ptr)).append (JU.Escape.eP (ptemp));
ptemp.setT (haveInversion ? ptinv : pta00);
if (ptemp.distance (pt0) < 0.1) ptemp.set (Math.random (), Math.random (), Math.random ());
draw1.append (JU.Escape.eP (ptemp));
ptemp.set (0, ang, 0);
draw1.append (JU.Escape.eP (ptemp)).append (" color red");
draw1.append (drawid).append ("rotVector1 vector diameter 0.1 ").append (JU.Escape.eP (pa1)).append (JU.Escape.eP (vtemp)).append ("color ").append (color);
}if (isMirrorPlane) {
if (pta00.distance (ptref) > 0.2) draw1.append (drawid).append ("planeVector arrow ").append (JU.Escape.eP (pta00)).append (JU.Escape.eP (ptref)).append (" color indigo");
opType = drawid + "plane";
if (trans != null) {
this.drawFrameLine ("X", ptref, vt1, 0.15, ptemp, draw1, opType, "red");
this.drawFrameLine ("Y", ptref, vt2, 0.15, ptemp, draw1, opType, "green");
this.drawFrameLine ("Z", ptref, vt3, 0.15, ptemp, draw1, opType, "blue");
opType = drawid + "glide";
}color = (trans == null ? "green" : "blue");
vtemp.setT (ax1);
vtemp.normalize ();
var w = -vtemp.x * pa1.x - vtemp.y * pa1.y - vtemp.z * pa1.z;
var plane = JU.P4.new4 (vtemp.x, vtemp.y, vtemp.z, w);
var v =  new JU.Lst ();
var margin = 1.05;
v.addLast (uc.getCanonicalCopy (margin, false));
modelSet.intersectPlane (plane, v, 3);
for (var i = v.size (); --i >= 0; ) {
var pts = v.get (i);
draw1.append (drawid).append ("planep").appendI (i).append (" ").append (JU.Escape.eP (pts[0])).append (JU.Escape.eP (pts[1]));
if (pts.length == 3) draw1.append (JU.Escape.eP (pts[2]));
draw1.append (" color translucent ").append (color);
}
if (v.size () == 0) {
ptemp.add2 (pa1, ax1);
draw1.append (drawid).append ("planeCircle scale 2.0 circle ").append (JU.Escape.eP (pa1)).append (JU.Escape.eP (ptemp)).append (" color translucent ").append (color).append (" mesh fill");
}}if (haveInversion) {
opType = drawid + "inv";
draw1.append (drawid).append ("invPoint diameter 0.4 ").append (JU.Escape.eP (ipt));
draw1.append (drawid).append ("invArrow arrow ").append (JU.Escape.eP (pta00)).append (JU.Escape.eP (ptinv)).append (" color indigo");
if (!isInversionOnly && !haveCentering) {
this.drawFrameLine ("X", ptinv, vt1, 0.15, ptemp, draw1, opType, "red");
this.drawFrameLine ("Y", ptinv, vt2, 0.15, ptemp, draw1, opType, "green");
this.drawFrameLine ("Z", ptinv, vt3, 0.15, ptemp, draw1, opType, "blue");
}}if (trans != null) {
if (ptref == null) ptref = JU.P3.newP (pta00);
draw1.append (drawid).append ("transVector vector ").append (JU.Escape.eP (ptref)).append (JU.Escape.eP (trans));
}if (haveCentering) {
if (opType != null) {
this.drawFrameLine ("X", pt0, vt1, 0.15, ptemp, draw1, opType, "red");
this.drawFrameLine ("Y", pt0, vt2, 0.15, ptemp, draw1, opType, "green");
this.drawFrameLine ("Z", pt0, vt3, 0.15, ptemp, draw1, opType, "blue");
}if (ptTarget == null) {
ptTarget = ptemp;
ptemp.add2 (pt0, vcentering);
}draw1.append (drawid).append ("centeringVector arrow ").append (JU.Escape.eP (pt0)).append (JU.Escape.eP (ptTarget)).append (" color cyan");
}ptemp2.setT (pt0);
if (haveCentering) ptemp2.add (vcentering);
ptemp.sub2 (pt1, pt0);
ptemp.scaleAdd2 (0.9, ptemp, ptemp2);
JS.SymmetryOperation.drawLine (draw1, drawid + "frame2X", 0.2, ptemp2, ptemp, "red");
ptemp.sub2 (pt2, pt0);
ptemp.scaleAdd2 (0.9, ptemp, ptemp2);
JS.SymmetryOperation.drawLine (draw1, drawid + "frame2Y", 0.2, ptemp2, ptemp, "green");
ptemp.sub2 (pt3, pt0);
ptemp.scaleAdd2 (0.9, ptemp, ptemp2);
JS.SymmetryOperation.drawLine (draw1, drawid + "frame2Z", 0.2, ptemp2, ptemp, "purple");
draw1.append ("\nvar pt00 = " + JU.Escape.eP (pta00));
draw1.append ("\nsym_point = pt00");
draw1.append ("\nvar p0 = " + JU.Escape.eP (ptemp2));
draw1.append ("\nvar set2 = within(0.2,p0);if(!set2){set2 = within(0.2,p0.uxyz.xyz)}");
draw1.append ("\nsym_target = set2;if (set2) {");
draw1.append (drawid).append ("offsetFrameX diameter 0.20 @{set2.xyz} @{set2.xyz + ").append (JU.Escape.eP (vt1)).append ("*0.9} color red");
draw1.append (drawid).append ("offsetFrameY diameter 0.20 @{set2.xyz} @{set2.xyz + ").append (JU.Escape.eP (vt2)).append ("*0.9} color green");
draw1.append (drawid).append ("offsetFrameZ diameter 0.20 @{set2.xyz} @{set2.xyz + ").append (JU.Escape.eP (vt3)).append ("*0.9} color purple");
draw1.append ("\n}}\n");
cmds = draw1.toString ();
draw1 = null;
drawid = null;
}if (trans == null) ftrans = null;
if (isRotation) {
if (haveInversion) {
} else if (pitch1 == 0) {
} else {
trans = JU.V3.newV (ax1);
ptemp.setT (trans);
uc.toFractional (ptemp, false);
ftrans = JU.V3.newV (ptemp);
}}if (isMirrorPlane) {
ang1 = 0;
}if (haveInversion) {
if (isInversionOnly) {
pa1 = null;
ax1 = null;
trans = null;
ftrans = null;
}} else if (isTranslation) {
pa1 = null;
ax1 = null;
}if (ax1 != null) ax1.normalize ();
var m2 = null;
m2 = JU.M4.newM4 (this);
if (haveCentering) vtrans.add (this.centering);
if (vtrans.length () != 0) {
m2.m03 += vtrans.x;
m2.m13 += vtrans.y;
m2.m23 += vtrans.z;
}xyzNew = (this.isBio ? m2.toString () : this.modDim > 0 ? this.xyzOriginal : JS.SymmetryOperation.getXYZFromMatrix (m2, false, false, false));
if (this.timeReversal != 0) xyzNew += (this.timeReversal == 1 ? ",m+1" : ",m-1");
return [xyzNew, this.xyzOriginal, info1, cmds, JS.SymmetryOperation.approx0 (ftrans), JS.SymmetryOperation.approx0 (trans), JS.SymmetryOperation.approx0 (ipt), JS.SymmetryOperation.approx0 (pa1), JS.SymmetryOperation.approx0 (ax1), Integer.$valueOf (ang1), m2, vtrans, this.centering];
}, "JM.ModelSet,J.api.SymmetryInterface,JU.P3,JU.P3,~S");
c$.setFractional = Clazz_defineMethod (c$, "setFractional", 
 function (uc, pt01, pt00, offset) {
pt01.setT (pt00);
if (offset != null) uc.toUnitCell (pt01, offset);
uc.toFractional (pt01, false);
}, "J.api.SymmetryInterface,JU.P3,JU.T3,JU.P3");
Clazz_defineMethod (c$, "drawFrameLine", 
 function (xyz, pt, v, width, ptemp, draw1, key, color) {
ptemp.setT (pt);
ptemp.add (v);
JS.SymmetryOperation.drawLine (draw1, key + "Frame" + xyz, width, pt, ptemp, "translucent " + color);
}, "~S,JU.P3,JU.V3,~N,JU.P3,JU.SB,~S,~S");
Clazz_defineMethod (c$, "rotTransCart", 
 function (pt00, uc, vtrans) {
var p0 = JU.P3.newP (pt00);
uc.toFractional (p0, false);
this.rotTrans2 (p0, p0);
p0.add (vtrans);
uc.toCartesian (p0, false);
return p0;
}, "JU.P3,J.api.SymmetryInterface,JU.V3");
Clazz_defineMethod (c$, "strCoord", 
 function (p) {
JS.SymmetryOperation.approx0 (p);
return (this.isBio ? p.x + " " + p.y + " " + p.z : JS.SymmetryOperation.fcoord (p));
}, "JU.T3");
c$.drawLine = Clazz_defineMethod (c$, "drawLine", 
 function (s, id, diameter, pt0, pt1, color) {
s.append (id).append (" diameter ").appendF (diameter).append (JU.Escape.eP (pt0)).append (JU.Escape.eP (pt1)).append (" color ").append (color);
}, "JU.SB,~S,~N,JU.P3,JU.P3,~S");
c$.fcoord = Clazz_defineMethod (c$, "fcoord", 
function (p) {
return JS.SymmetryOperation.fc (p.x) + " " + JS.SymmetryOperation.fc (p.y) + " " + JS.SymmetryOperation.fc (p.z);
}, "JU.T3");
c$.fc = Clazz_defineMethod (c$, "fc", 
 function (x) {
var xabs = Math.abs (x);
var x24 = Clazz_floatToInt (JS.SymmetryOperation.approxF (xabs * 24));
var m = (x < 0 ? "-" : "");
if (x24 % 8 != 0) return m + JS.SymmetryOperation.twelfthsOf (x24 >> 1);
return (x24 == 0 ? "0" : x24 == 24 ? m + "1" : m + (Clazz_doubleToInt (x24 / 8)) + "/3");
}, "~N");
c$.approx0 = Clazz_defineMethod (c$, "approx0", 
 function (pt) {
if (pt != null) {
if (Math.abs (pt.x) < 0.0001) pt.x = 0;
if (Math.abs (pt.y) < 0.0001) pt.y = 0;
if (Math.abs (pt.z) < 0.0001) pt.z = 0;
}return pt;
}, "JU.T3");
c$.approx = Clazz_defineMethod (c$, "approx", 
 function (pt) {
if (pt != null) {
pt.x = JS.SymmetryOperation.approxF (pt.x);
pt.y = JS.SymmetryOperation.approxF (pt.y);
pt.z = JS.SymmetryOperation.approxF (pt.z);
}return pt;
}, "JU.T3");
c$.approxF = Clazz_defineMethod (c$, "approxF", 
 function (f) {
return JU.PT.approx (f, 100);
}, "~N");
c$.normalizeTranslation = Clazz_defineMethod (c$, "normalizeTranslation", 
function (operation) {
operation.m03 = (Clazz_floatToInt (operation.m03) + 12) % 12;
operation.m13 = (Clazz_floatToInt (operation.m13) + 12) % 12;
operation.m23 = (Clazz_floatToInt (operation.m23) + 12) % 12;
}, "JU.M4");
c$.getXYZFromRsVs = Clazz_defineMethod (c$, "getXYZFromRsVs", 
function (rs, vs, is12ths) {
var ra = rs.getArray ();
var va = vs.getArray ();
var d = ra.length;
var s = "";
for (var i = 0; i < d; i++) {
s += ",";
for (var j = 0; j < d; j++) {
var r = ra[i][j];
if (r != 0) {
s += (r < 0 ? "-" : s.endsWith (",") ? "" : "+") + (Math.abs (r) == 1 ? "" : "" + Clazz_doubleToInt (Math.abs (r))) + "x" + (j + 1);
}}
s += JS.SymmetryOperation.xyzFraction (Clazz_doubleToInt (va[i][0] * (is12ths ? 1 : 12)), false, true);
}
return JU.PT.rep (s.substring (1), ",+", ",");
}, "JU.Matrix,JU.Matrix,~B");
Clazz_defineMethod (c$, "toString", 
function () {
return (this.rsvs == null ? Clazz_superCall (this, JS.SymmetryOperation, "toString", []) : Clazz_superCall (this, JS.SymmetryOperation, "toString", []) + " " + this.rsvs.toString ());
});
Clazz_defineMethod (c$, "getSpinOp", 
function () {
if (this.magOp == 3.4028235E38) this.magOp = this.determinant3 () * this.timeReversal;
return this.magOp;
});
Clazz_defineMethod (c$, "setTimeReversal", 
function (magRev) {
this.timeReversal = magRev;
if (this.xyz.indexOf (",m") >= 0) this.xyz = this.xyz.substring (0, this.xyz.indexOf (",m"));
this.xyz += (magRev == 1 ? ",m+1" : ",m-1");
}, "~N");
c$.cleanMatrix = Clazz_defineMethod (c$, "cleanMatrix", 
function (m4) {
var sb =  new JU.SB ();
sb.append ("[ ");
var row =  Clazz_newFloatArray (4, 0);
for (var i = 0; i < 3; i++) {
m4.getRow (i, row);
sb.append ("[ ").appendI (Clazz_floatToInt (row[0])).append (" ").appendI (Clazz_floatToInt (row[1])).append (" ").appendI (Clazz_floatToInt (row[2])).append (" ");
sb.append (JS.SymmetryOperation.twelfthsOf (row[3] * 12)).append (" ]");
}
return sb.append (" ]").toString ();
}, "JU.M4");
Clazz_defineMethod (c$, "setCentering", 
function (c, isFinal) {
if (this.centering == null && !this.unCentered) {
if (this.modDim == 0 && this.index > 1 && this.m00 == 1 && this.m11 == 1 && this.m22 == 1 && this.m01 == 0 && this.m02 == 0 && this.m10 == 0 && this.m12 == 0 && this.m20 == 0 && this.m21 == 0) {
this.centering = JU.V3.new3 (this.m03, this.m13, this.m23);
if (this.centering.lengthSquared () == 0) {
this.unCentered = true;
this.centering = null;
} else if (!isFinal) this.centering.scale (0.083333336);
this.isCenteringOp = true;
} else {
this.centering = c;
}}return this.centering;
}, "JU.V3,~B");
Clazz_defineStatics (c$,
"twelfths", ["0", "1/12", "1/6", "1/4", "1/3", "5/12", "1/2", "7/12", "2/3", "3/4", "5/6", "11/12"]);
c$.labelsXYZ = c$.prototype.labelsXYZ = ["x", "y", "z"];
c$.labelsXn = c$.prototype.labelsXn = ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13"];
c$.labelsXnSub = c$.prototype.labelsXnSub = ["x", "y", "z", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j"];
});
Clazz_declarePackage ("JS");
Clazz_load (null, "JS.SymmetryInfo", ["JU.Lst", "$.PT", "JU.SimpleUnitCell"], function () {
c$ = Clazz_decorateAsClass (function () {
this.coordinatesAreFractional = false;
this.isMultiCell = false;
this.sgName = null;
this.symmetryOperations = null;
this.infoStr = null;
this.cellRange = null;
this.periodicOriginXyz = null;
this.centerings = null;
Clazz_instantialize (this, arguments);
}, JS, "SymmetryInfo");
Clazz_defineMethod (c$, "isPeriodic", 
function () {
return this.periodicOriginXyz != null;
});
Clazz_makeConstructor (c$, 
function () {
});
Clazz_defineMethod (c$, "setSymmetryInfo", 
function (info, notionalUnitCell) {
this.cellRange = info.get ("unitCellRange");
this.periodicOriginXyz = info.get ("periodicOriginXyz");
this.sgName = info.get ("spaceGroup");
if (this.sgName == null || this.sgName === "") this.sgName = "spacegroup unspecified";
var symmetryCount = info.containsKey ("symmetryCount") ? (info.get ("symmetryCount")).intValue () : 0;
this.symmetryOperations = info.remove ("symmetryOps");
this.infoStr = "Spacegroup: " + this.sgName;
if (this.symmetryOperations == null) {
this.infoStr += "\nNumber of symmetry operations: ?\nSymmetry Operations: unspecified\n";
} else {
this.centerings =  new JU.Lst ();
var c = "";
var s = "\nNumber of symmetry operations: " + (symmetryCount == 0 ? 1 : symmetryCount) + "\nSymmetry Operations:";
for (var i = 0; i < symmetryCount; i++) {
var op = this.symmetryOperations[i];
s += "\n" + op.xyz;
if (op.isCenteringOp) {
this.centerings.addLast (op.centering);
var oc = JU.PT.replaceAllCharacters (op.xyz, "xyz", "0");
c += " (" + JU.PT.rep (oc, "0+", "") + ")";
}}
if (c.length > 0) this.infoStr += "\nCentering: " + c;
this.infoStr += s;
}this.infoStr += "\n";
if (notionalUnitCell == null) notionalUnitCell = info.get ("notionalUnitcell");
if (!JU.SimpleUnitCell.isValid (notionalUnitCell)) return null;
this.coordinatesAreFractional = info.containsKey ("coordinatesAreFractional") ? (info.get ("coordinatesAreFractional")).booleanValue () : false;
this.isMultiCell = (this.coordinatesAreFractional && this.symmetryOperations != null);
return notionalUnitCell;
}, "java.util.Map,~A");
});
Clazz_declarePackage ("JS");
Clazz_load (["JU.SimpleUnitCell", "JU.P3", "JV.JC"], "JS.UnitCell", ["java.lang.Float", "JU.M4", "$.Quat", "$.T4", "$.V3", "J.api.Interface", "JU.BoxInfo", "$.Escape"], function () {
c$ = Clazz_decorateAsClass (function () {
this.vertices = null;
this.cartesianOffset = null;
this.fractionalOffset = null;
this.allFractionalRelative = false;
this.unitCellMultiplier = null;
this.moreInfo = null;
this.name = "";
Clazz_instantialize (this, arguments);
}, JS, "UnitCell", JU.SimpleUnitCell);
Clazz_prepareFields (c$, function () {
this.cartesianOffset =  new JU.P3 ();
});
Clazz_makeConstructor (c$, 
function () {
Clazz_superConstructor (this, JS.UnitCell, []);
});
c$.newP = Clazz_defineMethod (c$, "newP", 
function (points, setRelative) {
var c =  new JS.UnitCell ();
var parameters = [-1, 0, 0, 0, 0, 0, points[1].x, points[1].y, points[1].z, points[2].x, points[2].y, points[2].z, points[3].x, points[3].y, points[3].z];
c.init (parameters);
c.allFractionalRelative = setRelative;
c.initUnitcellVertices ();
c.setCartesianOffset (points[0]);
return c;
}, "~A,~B");
c$.newA = Clazz_defineMethod (c$, "newA", 
function (notionalUnitcell, setRelative) {
var c =  new JS.UnitCell ();
c.init (notionalUnitcell);
c.initUnitcellVertices ();
c.allFractionalRelative = setRelative;
return c;
}, "~A,~B");
Clazz_defineMethod (c$, "initOrientation", 
function (mat) {
if (mat == null) return;
var m =  new JU.M4 ();
m.setToM3 (mat);
this.matrixFractionalToCartesian.mul2 (m, this.matrixFractionalToCartesian);
this.matrixCartesianToFractional.invertM (this.matrixFractionalToCartesian);
this.initUnitcellVertices ();
}, "JU.M3");
Clazz_defineMethod (c$, "toUnitCell", 
function (pt, offset) {
if (this.matrixCartesianToFractional == null) return;
if (offset == null) {
this.matrixCartesianToFractional.rotTrans (pt);
this.unitize (pt);
this.matrixFractionalToCartesian.rotTrans (pt);
} else {
this.matrixCtoFANoOffset.rotTrans (pt);
this.unitize (pt);
pt.add (offset);
this.matrixFtoCNoOffset.rotTrans (pt);
}}, "JU.P3,JU.P3");
Clazz_defineMethod (c$, "unitize", 
function (pt) {
switch (this.dimension) {
case 3:
pt.z = JS.UnitCell.toFractionalX (pt.z);
case 2:
pt.y = JS.UnitCell.toFractionalX (pt.y);
case 1:
pt.x = JS.UnitCell.toFractionalX (pt.x);
}
}, "JU.P3");
Clazz_defineMethod (c$, "reset", 
function () {
this.unitCellMultiplier = null;
this.setOffset (JU.P3.new3 (0, 0, 0));
});
Clazz_defineMethod (c$, "setOffset", 
function (pt) {
if (pt == null) return;
var pt4 = (Clazz_instanceOf (pt, JU.T4) ? pt : null);
if (pt4 != null ? pt4.w <= 0 : pt.x >= 100 || pt.y >= 100) {
this.unitCellMultiplier = (pt.z == 0 && pt.x == pt.y ? null : JU.P3.newP (pt));
if (pt4 == null || pt4.w == 0) return;
}if (this.hasOffset () || pt.lengthSquared () > 0) {
this.fractionalOffset =  new JU.P3 ();
this.fractionalOffset.setT (pt);
}this.matrixCartesianToFractional.m03 = -pt.x;
this.matrixCartesianToFractional.m13 = -pt.y;
this.matrixCartesianToFractional.m23 = -pt.z;
this.cartesianOffset.setT (pt);
this.matrixFractionalToCartesian.m03 = 0;
this.matrixFractionalToCartesian.m13 = 0;
this.matrixFractionalToCartesian.m23 = 0;
this.matrixFractionalToCartesian.rotTrans (this.cartesianOffset);
this.matrixFractionalToCartesian.m03 = this.cartesianOffset.x;
this.matrixFractionalToCartesian.m13 = this.cartesianOffset.y;
this.matrixFractionalToCartesian.m23 = this.cartesianOffset.z;
if (this.allFractionalRelative) {
this.matrixCtoFANoOffset.setM4 (this.matrixCartesianToFractional);
this.matrixFtoCNoOffset.setM4 (this.matrixFractionalToCartesian);
}}, "JU.T3");
Clazz_defineMethod (c$, "setCartesianOffset", 
 function (origin) {
this.cartesianOffset.setT (origin);
this.matrixFractionalToCartesian.m03 = this.cartesianOffset.x;
this.matrixFractionalToCartesian.m13 = this.cartesianOffset.y;
this.matrixFractionalToCartesian.m23 = this.cartesianOffset.z;
var wasOffset = this.hasOffset ();
this.fractionalOffset =  new JU.P3 ();
this.fractionalOffset.setT (this.cartesianOffset);
this.matrixCartesianToFractional.m03 = 0;
this.matrixCartesianToFractional.m13 = 0;
this.matrixCartesianToFractional.m23 = 0;
this.matrixCartesianToFractional.rotTrans (this.fractionalOffset);
this.matrixCartesianToFractional.m03 = -this.fractionalOffset.x;
this.matrixCartesianToFractional.m13 = -this.fractionalOffset.y;
this.matrixCartesianToFractional.m23 = -this.fractionalOffset.z;
if (this.allFractionalRelative) {
this.matrixCtoFANoOffset.setM4 (this.matrixCartesianToFractional);
this.matrixFtoCNoOffset.setM4 (this.matrixFractionalToCartesian);
}if (!wasOffset && this.fractionalOffset.lengthSquared () == 0) this.fractionalOffset = null;
}, "JU.T3");
Clazz_defineMethod (c$, "setMinMaxLatticeParameters", 
function (minXYZ, maxXYZ) {
if (maxXYZ.x <= maxXYZ.y && maxXYZ.y >= 555) {
var pt =  new JU.P3 ();
JU.SimpleUnitCell.ijkToPoint3f (maxXYZ.x, pt, 0);
minXYZ.x = Clazz_floatToInt (pt.x);
minXYZ.y = Clazz_floatToInt (pt.y);
minXYZ.z = Clazz_floatToInt (pt.z);
JU.SimpleUnitCell.ijkToPoint3f (maxXYZ.y, pt, 1);
maxXYZ.x = Clazz_floatToInt (pt.x);
maxXYZ.y = Clazz_floatToInt (pt.y);
maxXYZ.z = Clazz_floatToInt (pt.z);
}switch (this.dimension) {
case 1:
minXYZ.y = 0;
maxXYZ.y = 1;
case 2:
minXYZ.z = 0;
maxXYZ.z = 1;
}
}, "JU.P3i,JU.P3i");
Clazz_defineMethod (c$, "dumpInfo", 
function (isFull) {
return "a=" + this.a + ", b=" + this.b + ", c=" + this.c + ", alpha=" + this.alpha + ", beta=" + this.beta + ", gamma=" + this.gamma + "\n" + JU.Escape.eAP (this.getUnitCellVectors ()) + "\nvolume=" + this.volume + (isFull ? "\nfractional to cartesian: " + this.matrixFractionalToCartesian + "\ncartesian to fractional: " + this.matrixCartesianToFractional : "");
}, "~B");
Clazz_defineMethod (c$, "getVertices", 
function () {
return this.vertices;
});
Clazz_defineMethod (c$, "getCartesianOffset", 
function () {
return this.cartesianOffset;
});
Clazz_defineMethod (c$, "getFractionalOffset", 
function () {
return this.fractionalOffset;
});
Clazz_defineMethod (c$, "getTensor", 
function (parBorU) {
var t = (J.api.Interface.getUtil ("Tensor"));
if (parBorU[0] == 0 && parBorU[1] == 0 && parBorU[2] == 0) {
var f = parBorU[7];
var eigenValues = [f, f, f];
return t.setFromEigenVectors (JS.UnitCell.unitVectors, eigenValues, "iso", "Uiso=" + f, null);
}t.parBorU = parBorU;
var Bcart =  Clazz_newDoubleArray (6, 0);
var ortepType = Clazz_floatToInt (parBorU[6]);
if (ortepType == 12) {
Bcart[0] = parBorU[0] * 19.739208802178716;
Bcart[1] = parBorU[1] * 19.739208802178716;
Bcart[2] = parBorU[2] * 19.739208802178716;
Bcart[3] = parBorU[3] * 19.739208802178716 * 2;
Bcart[4] = parBorU[4] * 19.739208802178716 * 2;
Bcart[5] = parBorU[5] * 19.739208802178716 * 2;
parBorU[7] = (parBorU[0] + parBorU[1] + parBorU[3]) / 3;
} else {
var isFractional = (ortepType == 4 || ortepType == 5 || ortepType == 8 || ortepType == 9);
var cc = 2 - (ortepType % 2);
var dd = (ortepType == 8 || ortepType == 9 || ortepType == 10 ? 19.739208802178716 : ortepType == 4 || ortepType == 5 ? 0.25 : ortepType == 2 || ortepType == 3 ? Math.log (2) : 1);
var B11 = parBorU[0] * dd * (isFractional ? this.a_ * this.a_ : 1);
var B22 = parBorU[1] * dd * (isFractional ? this.b_ * this.b_ : 1);
var B33 = parBorU[2] * dd * (isFractional ? this.c_ * this.c_ : 1);
var B12 = parBorU[3] * dd * (isFractional ? this.a_ * this.b_ : 1) * cc;
var B13 = parBorU[4] * dd * (isFractional ? this.a_ * this.c_ : 1) * cc;
var B23 = parBorU[5] * dd * (isFractional ? this.b_ * this.c_ : 1) * cc;
parBorU[7] = Math.pow (B11 / 19.739208802178716 / this.a_ / this.a_ * B22 / 19.739208802178716 / this.b_ / this.b_ * B33 / 19.739208802178716 / this.c_ / this.c_, 0.3333);
Bcart[0] = this.a * this.a * B11 + this.b * this.b * this.cosGamma * this.cosGamma * B22 + this.c * this.c * this.cosBeta * this.cosBeta * B33 + this.a * this.b * this.cosGamma * B12 + this.b * this.c * this.cosGamma * this.cosBeta * B23 + this.a * this.c * this.cosBeta * B13;
Bcart[1] = this.b * this.b * this.sinGamma * this.sinGamma * B22 + this.c * this.c * this.cA_ * this.cA_ * B33 + this.b * this.c * this.cA_ * this.sinGamma * B23;
Bcart[2] = this.c * this.c * this.cB_ * this.cB_ * B33;
Bcart[3] = 2 * this.b * this.b * this.cosGamma * this.sinGamma * B22 + 2 * this.c * this.c * this.cA_ * this.cosBeta * B33 + this.a * this.b * this.sinGamma * B12 + this.b * this.c * (this.cA_ * this.cosGamma + this.sinGamma * this.cosBeta) * B23 + this.a * this.c * this.cA_ * B13;
Bcart[4] = 2 * this.c * this.c * this.cB_ * this.cosBeta * B33 + this.b * this.c * this.cosGamma * B23 + this.a * this.c * this.cB_ * B13;
Bcart[5] = 2 * this.c * this.c * this.cA_ * this.cB_ * B33 + this.b * this.c * this.cB_ * this.sinGamma * B23;
}return t.setFromThermalEquation (Bcart, JU.Escape.eAF (parBorU));
}, "~A");
Clazz_defineMethod (c$, "getCanonicalCopy", 
function (scale, withOffset) {
var pts =  new Array (8);
var cell0 = null;
var cell1 = null;
if (withOffset && this.unitCellMultiplier != null) {
cell0 =  new JU.P3 ();
cell1 =  new JU.P3 ();
JU.SimpleUnitCell.ijkToPoint3f (Clazz_floatToInt (this.unitCellMultiplier.x), cell0, 0);
JU.SimpleUnitCell.ijkToPoint3f (Clazz_floatToInt (this.unitCellMultiplier.y), cell1, 0);
cell1.sub (cell0);
}for (var i = 0; i < 8; i++) {
var pt = pts[i] = JU.P3.newP (JU.BoxInfo.unitCubePoints[i]);
if (cell0 != null) {
scale *= (this.unitCellMultiplier.z == 0 ? 1 : this.unitCellMultiplier.z);
pts[i].add3 (cell0.x + cell1.x * pt.x, cell0.y + cell1.y * pt.y, cell0.z + cell1.z * pt.z);
}this.matrixFractionalToCartesian.rotTrans (pt);
if (!withOffset) pt.sub (this.cartesianOffset);
}
return JU.BoxInfo.getCanonicalCopy (pts, scale);
}, "~N,~B");
c$.toFractionalX = Clazz_defineMethod (c$, "toFractionalX", 
 function (x) {
x = (x - Math.floor (x));
if (x > 0.9999 || x < 0.0001) x = 0;
return x;
}, "~N");
Clazz_defineMethod (c$, "initUnitcellVertices", 
 function () {
if (this.matrixFractionalToCartesian == null) return;
this.matrixCtoFANoOffset = JU.M4.newM4 (this.matrixCartesianToFractional);
this.matrixFtoCNoOffset = JU.M4.newM4 (this.matrixFractionalToCartesian);
this.vertices =  new Array (8);
for (var i = 8; --i >= 0; ) {
this.vertices[i] =  new JU.P3 ();
this.matrixFractionalToCartesian.rotTrans2 (JU.BoxInfo.unitCubePoints[i], this.vertices[i]);
}
});
Clazz_defineMethod (c$, "checkDistance", 
function (f1, f2, distance, dx, iRange, jRange, kRange, ptOffset) {
var p1 = JU.P3.newP (f1);
this.toCartesian (p1, true);
for (var i = -iRange; i <= iRange; i++) for (var j = -jRange; j <= jRange; j++) for (var k = -kRange; k <= kRange; k++) {
ptOffset.set (f2.x + i, f2.y + j, f2.z + k);
this.toCartesian (ptOffset, true);
var d = p1.distance (ptOffset);
if (dx > 0 ? Math.abs (d - distance) <= dx : d <= distance && d > 0.1) {
ptOffset.set (i, j, k);
return true;
}}


return false;
}, "JU.P3,JU.P3,~N,~N,~N,~N,~N,JU.P3");
Clazz_defineMethod (c$, "getUnitCellMultiplier", 
function () {
return this.unitCellMultiplier;
});
Clazz_defineMethod (c$, "getUnitCellVectors", 
function () {
var m = this.matrixFractionalToCartesian;
return [JU.V3.newV (this.cartesianOffset), JU.V3.new3 (this.fix (m.m00), this.fix (m.m10), this.fix (m.m20)), JU.V3.new3 (this.fix (m.m01), this.fix (m.m11), this.fix (m.m21)), JU.V3.new3 (this.fix (m.m02), this.fix (m.m12), this.fix (m.m22))];
});
Clazz_defineMethod (c$, "fix", 
 function (x) {
return (Math.abs (x) < 0.001 ? 0 : x);
}, "~N");
Clazz_defineMethod (c$, "isSameAs", 
function (uc) {
if (uc.notionalUnitcell.length != this.notionalUnitcell.length) return false;
for (var i = this.notionalUnitcell.length; --i >= 0; ) if (this.notionalUnitcell[i] != uc.notionalUnitcell[i] && !(Float.isNaN (this.notionalUnitcell[i]) && Float.isNaN (uc.notionalUnitcell[i]))) return false;

return (this.fractionalOffset == null ? !uc.hasOffset () : uc.fractionalOffset == null ? !this.hasOffset () : this.fractionalOffset.distanceSquared (uc.fractionalOffset) == 0);
}, "JS.UnitCell");
Clazz_defineMethod (c$, "hasOffset", 
function () {
return (this.fractionalOffset != null && this.fractionalOffset.lengthSquared () != 0);
});
Clazz_defineMethod (c$, "getState", 
function () {
var s = "";
if (this.fractionalOffset != null && this.fractionalOffset.lengthSquared () != 0) s += "  unitcell offset " + JU.Escape.eP (this.fractionalOffset) + ";\n";
if (this.unitCellMultiplier != null) s += "  unitcell range " + JU.Escape.eP (this.unitCellMultiplier) + ";\n";
return s;
});
Clazz_defineMethod (c$, "getQuaternionRotation", 
function (abc) {
var a = this.vertices[4];
var b = this.vertices[2];
var c = this.vertices[1];
var x =  new JU.V3 ();
var v =  new JU.V3 ();
switch ("abc".indexOf (abc)) {
case 0:
x.cross (a, c);
v.cross (x, a);
break;
case 1:
x.cross (b, a);
v.cross (x, b);
break;
case 2:
x.cross (c, b);
v.cross (x, c);
break;
default:
return null;
}
return JU.Quat.getQuaternionFrame (null, v, x).inv ();
}, "~S");
Clazz_defineStatics (c$,
"twoP2", 19.739208802178716);
c$.unitVectors = c$.prototype.unitVectors = [JV.JC.axisX, JV.JC.axisY, JV.JC.axisZ];
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
