Clazz.declarePackage ("JM");
Clazz.load (["J.api.JmolBioResolver"], "JM.Resolver", ["java.lang.Boolean", "$.NullPointerException", "java.util.Arrays", "$.Hashtable", "JU.BS", "$.P3", "$.P4", "$.PT", "$.SB", "$.V3", "J.c.STR", "JM.Group", "JM.AlphaMonomer", "$.AlphaPolymer", "$.AminoMonomer", "$.AminoPolymer", "$.BioModel", "$.CarbohydrateMonomer", "$.CarbohydratePolymer", "$.Monomer", "$.NucleicMonomer", "$.NucleicPolymer", "$.PhosphorusMonomer", "$.PhosphorusPolymer", "JU.BSUtil", "$.Logger", "$.Measure", "JV.JC"], function () {
c$ = Clazz.decorateAsClass (function () {
this.ml = null;
this.ms = null;
this.bsAddedHydrogens = null;
this.bsAtomsForHs = null;
this.htBondMap = null;
this.htGroupBonds = null;
this.hNames = null;
this.lastSetH = -2147483648;
this.maxSerial = 0;
this.baseBondIndex = 0;
this.haveHsAlready = false;
this.vAB = null;
this.vAC = null;
this.vNorm = null;
this.plane = null;
if (!Clazz.isClassDefined ("JM.Resolver.BondSorter")) {
JM.Resolver.$Resolver$BondSorter$ ();
}
this.bsAssigned = null;
Clazz.instantialize (this, arguments);
}, JM, "Resolver", null, J.api.JmolBioResolver);
Clazz.makeConstructor (c$, 
function () {
});
Clazz.overrideMethod (c$, "getBioModel", 
function (modelIndex, trajectoryBaseIndex, jmolData, modelProperties, modelAuxiliaryInfo) {
return  new JM.BioModel (this.ms, modelIndex, trajectoryBaseIndex, jmolData, modelProperties, modelAuxiliaryInfo);
}, "~N,~N,~S,java.util.Properties,java.util.Map");
Clazz.overrideMethod (c$, "distinguishAndPropagateGroup", 
function (chain, group3, seqcode, firstAtomIndex, maxAtomIndex, modelIndex, specialAtomIndexes, atoms) {
var lastAtomIndex = maxAtomIndex - 1;
var distinguishingBits = 0;
for (var i = JV.JC.ATOMID_MAX; --i >= 0; ) specialAtomIndexes[i] = -2147483648;

for (var i = maxAtomIndex; --i >= firstAtomIndex; ) {
var specialAtomID = atoms[i].getAtomID ();
if (specialAtomID <= 0) continue;
if (specialAtomID < 14) {
distinguishingBits |= (1 << specialAtomID);
}specialAtomIndexes[specialAtomID] = i;
}
if (lastAtomIndex < firstAtomIndex) throw  new NullPointerException ();
var m = null;
if ((distinguishingBits & 14) == 14) m = JM.AminoMonomer.validateAndAllocate (chain, group3, seqcode, firstAtomIndex, lastAtomIndex, specialAtomIndexes, atoms);
 else if (distinguishingBits == 4) m = JM.AlphaMonomer.validateAndAllocateA (chain, group3, seqcode, firstAtomIndex, lastAtomIndex, specialAtomIndexes);
 else if (((distinguishingBits & 8128) == 8128)) m = JM.NucleicMonomer.validateAndAllocate (chain, group3, seqcode, firstAtomIndex, lastAtomIndex, specialAtomIndexes);
 else if (distinguishingBits == 8192) m = JM.PhosphorusMonomer.validateAndAllocateP (chain, group3, seqcode, firstAtomIndex, lastAtomIndex, specialAtomIndexes);
 else if (JV.JC.checkCarbohydrate (group3)) m = JM.CarbohydrateMonomer.validateAndAllocate (chain, group3, seqcode, firstAtomIndex, lastAtomIndex);
return (m != null && m.leadAtomIndex >= 0 ? m : null);
}, "JM.Chain,~S,~N,~N,~N,~N,~A,~A");
Clazz.overrideMethod (c$, "setHaveHsAlready", 
function (b) {
this.haveHsAlready = b;
}, "~B");
Clazz.overrideMethod (c$, "initialize", 
function (modelLoader) {
this.ml = modelLoader;
this.ms = modelLoader.ms;
}, "JM.ModelLoader");
Clazz.overrideMethod (c$, "initializeHydrogenAddition", 
function () {
this.baseBondIndex = this.ml.ms.bondCount;
this.bsAddedHydrogens =  new JU.BS ();
this.bsAtomsForHs =  new JU.BS ();
this.htBondMap =  new java.util.Hashtable ();
this.htGroupBonds =  new java.util.Hashtable ();
this.hNames =  new Array (3);
this.vAB =  new JU.V3 ();
this.vAC =  new JU.V3 ();
this.vNorm =  new JU.V3 ();
this.plane =  new JU.P4 ();
});
Clazz.overrideMethod (c$, "addImplicitHydrogenAtoms", 
function (adapter, iGroup, nH) {
var group3 = this.ml.getGroup3 (iGroup);
var nH1;
if (this.haveHsAlready || group3 == null || (nH1 = JV.JC.getStandardPdbHydrogenCount (JM.Group.lookupGroupID (group3))) == 0) return;
nH = (nH1 < 0 ? -1 : nH1 + nH);
var model = null;
var iFirst = this.ml.getFirstAtomIndex (iGroup);
var ac = this.ms.getAtomCount ();
if (nH < 0) {
if (ac - iFirst == 1) return;
model = this.ms.vwr.getLigandModel (group3, "ligand_", "_data", null);
if (model == null) return;
nH = adapter.getHydrogenAtomCount (model);
if (nH < 1) return;
}this.getBondInfo (adapter, group3, model);
this.ms.am[this.ms.at[iFirst].mi].isPdbWithMultipleBonds = true;
this.bsAtomsForHs.setBits (iFirst, ac);
this.bsAddedHydrogens.setBits (ac, ac + nH);
var isHetero = this.ms.at[iFirst].isHetero ();
var xyz = JU.P3.new3 (NaN, NaN, NaN);
for (var i = 0; i < nH; i++) this.ms.addAtom (this.ms.at[iFirst].mi, this.ms.at[iFirst].getGroup (), 1, "H", 0, 0, xyz, NaN, null, 0, 0, 1, 0, null, isHetero, 0, null).deleteBonds (null);

}, "J.api.JmolAdapter,~N,~N");
Clazz.defineMethod (c$, "getBondInfo", 
function (adapter, group3, model) {
if (this.htGroupBonds.get (group3) != null) return;
var bondInfo;
if (model == null) {
bondInfo = this.ms.vwr.getPdbBondInfo (group3);
} else {
bondInfo = this.getLigandBondInfo (adapter, model, group3);
}if (bondInfo == null) return;
this.htGroupBonds.put (group3, Boolean.TRUE);
for (var i = 0; i < bondInfo.length; i++) {
if (bondInfo[i] == null) continue;
if (bondInfo[i][1].charAt (0) == 'H') this.htBondMap.put (group3 + "." + bondInfo[i][0], bondInfo[i][1]);
 else this.htBondMap.put (group3 + ":" + bondInfo[i][0] + ":" + bondInfo[i][1], bondInfo[i][2]);
}
}, "J.api.JmolAdapter,~S,~O");
Clazz.defineMethod (c$, "getLigandBondInfo", 
 function (adapter, model, group3) {
var dataIn = adapter.getBondList (model);
var htAtoms =  new java.util.Hashtable ();
var iterAtom = adapter.getAtomIterator (model);
while (iterAtom.hasNext ()) htAtoms.put (iterAtom.getAtomName (), iterAtom.getXYZ ());

var bondInfo =  new Array (dataIn.length * 2);
var n = 0;
for (var i = 0; i < dataIn.length; i++) {
var b = dataIn[i];
if (b[0].charAt (0) != 'H') bondInfo[n++] = [b[0], b[1], b[2], b[1].startsWith ("H") ? "0" : "1"];
if (b[1].charAt (0) != 'H') bondInfo[n++] = [b[1], b[0], b[2], b[0].startsWith ("H") ? "0" : "1"];
}
java.util.Arrays.sort (bondInfo, Clazz.innerTypeInstance (JM.Resolver.BondSorter, this, null));
var t;
for (var i = 0; i < n; ) {
t = bondInfo[i];
var a1 = t[0];
var nH = 0;
var nC = 0;
for (; i < n && (t = bondInfo[i])[0].equals (a1); i++) {
if (t[3].equals ("0")) {
nH++;
continue;
}if (t[3].equals ("1")) nC++;
}
var pt = i - nH - nC;
if (nH == 1) continue;
switch (nC) {
case 1:
var sep = (nH == 2 ? '@' : '|');
for (var j = 1; j < nH; j++) {
bondInfo[pt][1] += sep + bondInfo[pt + j][1];
bondInfo[pt + j] = null;
}
continue;
case 2:
if (nH != 2) continue;
var name = bondInfo[pt][0];
var name1 = bondInfo[pt + nH][1];
var name2 = bondInfo[pt + nH + 1][1];
var factor = name1.compareTo (name2);
JU.Measure.getPlaneThroughPoints (htAtoms.get (name1), htAtoms.get (name), htAtoms.get (name2), this.vNorm, this.vAB, this.vAC, this.plane);
var d = JU.Measure.distanceToPlane (this.plane, htAtoms.get (bondInfo[pt][1])) * factor;
bondInfo[pt][1] = (d > 0 ? bondInfo[pt][1] + "@" + bondInfo[pt + 1][1] : bondInfo[pt + 1][1] + "@" + bondInfo[pt][1]);
bondInfo[pt + 1] = null;
}
}
for (var i = 0; i < n; i++) {
if ((t = bondInfo[i]) != null && t[1].charAt (0) != 'H' && t[0].compareTo (t[1]) > 0) {
bondInfo[i] = null;
continue;
}if (t != null) JU.Logger.info (" ligand " + group3 + ": " + bondInfo[i][0] + " - " + bondInfo[i][1] + " order " + bondInfo[i][2]);
}
return bondInfo;
}, "J.api.JmolAdapter,~O,~S");
Clazz.overrideMethod (c$, "finalizeHydrogens", 
function () {
this.ms.vwr.getLigandModel (null, null, null, null);
this.finalizePdbMultipleBonds ();
this.addHydrogens ();
});
Clazz.defineMethod (c$, "addHydrogens", 
 function () {
if (this.bsAddedHydrogens.nextSetBit (0) < 0) return;
this.finalizePdbCharges ();
var nTotal =  Clazz.newIntArray (1, 0);
var pts = this.ms.calculateHydrogens (this.bsAtomsForHs, nTotal, true, false, null);
var groupLast = null;
var ipt = 0;
for (var i = 0; i < pts.length; i++) {
if (pts[i] == null) continue;
var atom = this.ms.at[i];
var g = atom.getGroup ();
if (g !== groupLast) {
groupLast = g;
ipt = g.lastAtomIndex;
while (this.bsAddedHydrogens.get (ipt)) ipt--;

}var gName = atom.getGroup3 (false);
var aName = atom.getAtomName ();
var hName = this.htBondMap.get (gName + "." + aName);
if (hName == null) continue;
var isChiral = hName.contains ("@");
var isMethyl = (hName.endsWith ("?") || hName.indexOf ("|") >= 0);
var n = pts[i].length;
if (n == 3 && !isMethyl && hName.equals ("H@H2")) {
hName = "H|H2|H3";
isMethyl = true;
isChiral = false;
}if (isChiral && n == 3 || isMethyl != (n == 3)) {
JU.Logger.info ("Error adding H atoms to " + gName + g.getResno () + ": " + pts[i].length + " atoms should not be added to " + aName);
continue;
}var pt = hName.indexOf ("@");
switch (pts[i].length) {
case 1:
if (pt > 0) hName = hName.substring (0, pt);
this.setHydrogen (i, ++ipt, hName, pts[i][0]);
break;
case 2:
var hName1;
var hName2;
var d = -1;
var bonds = atom.getBonds ();
if (bonds != null) switch (bonds.length) {
case 2:
var atom1 = bonds[0].getOtherAtom (atom);
var atom2 = bonds[1].getOtherAtom (atom);
var factor = atom1.getAtomName ().compareTo (atom2.getAtomName ());
JU.Measure.getPlaneThroughPoints (atom1, atom, atom2, this.vNorm, this.vAB, this.vAC, this.plane);
d = JU.Measure.distanceToPlane (this.plane, pts[i][0]) * factor;
break;
}
if (pt < 0) {
JU.Logger.info ("Error adding H atoms to " + gName + g.getResno () + ": expected to only need 1 H but needed 2");
hName1 = hName2 = "H";
} else if (d < 0) {
hName2 = hName.substring (0, pt);
hName1 = hName.substring (pt + 1);
} else {
hName1 = hName.substring (0, pt);
hName2 = hName.substring (pt + 1);
}this.setHydrogen (i, ++ipt, hName1, pts[i][0]);
this.setHydrogen (i, ++ipt, hName2, pts[i][1]);
break;
case 3:
var pt1 = hName.indexOf ('|');
if (pt1 >= 0) {
var pt2 = hName.lastIndexOf ('|');
this.hNames[0] = hName.substring (0, pt1);
this.hNames[1] = hName.substring (pt1 + 1, pt2);
this.hNames[2] = hName.substring (pt2 + 1);
} else {
this.hNames[0] = hName.$replace ('?', '1');
this.hNames[1] = hName.$replace ('?', '2');
this.hNames[2] = hName.$replace ('?', '3');
}this.setHydrogen (i, ++ipt, this.hNames[0], pts[i][0]);
this.setHydrogen (i, ++ipt, this.hNames[1], pts[i][2]);
this.setHydrogen (i, ++ipt, this.hNames[2], pts[i][1]);
break;
}
}
this.deleteUnneededAtoms ();
this.ms.fixFormalCharges (JU.BSUtil.newBitSet2 (this.ml.baseAtomIndex, this.ml.ms.ac));
});
Clazz.defineMethod (c$, "deleteUnneededAtoms", 
 function () {
var bsBondsDeleted =  new JU.BS ();
for (var i = this.bsAtomsForHs.nextSetBit (0); i >= 0; i = this.bsAtomsForHs.nextSetBit (i + 1)) {
var atom = this.ms.at[i];
if (!atom.isHetero () || atom.getElementNumber () != 8 || atom.getFormalCharge () != 0 || atom.getCovalentBondCount () != 2) continue;
var bonds = atom.getBonds ();
var atom1 = bonds[0].getOtherAtom (atom);
var atomH = bonds[1].getOtherAtom (atom);
if (atom1.getElementNumber () == 1) {
var a = atom1;
atom1 = atomH;
atomH = a;
}if (atomH.getElementNumber () != 1) continue;
var bonds1 = atom1.getBonds ();
for (var j = 0; j < bonds1.length; j++) {
if (bonds1[j].order == 2) {
var atomO = bonds1[j].getOtherAtom (atom1);
if (atomO.getElementNumber () == 8) {
this.bsAddedHydrogens.set (atomH.i);
atomH.deleteBonds (bsBondsDeleted);
break;
}}}
}
this.ms.deleteBonds (bsBondsDeleted, true);
this.ml.deleteAtoms (this.bsAddedHydrogens);
});
Clazz.defineMethod (c$, "finalizePdbCharges", 
 function () {
var atoms = this.ms.at;
for (var i = this.bsAtomsForHs.nextSetBit (0); i >= 0; i = this.bsAtomsForHs.nextSetBit (i + 1)) {
var a = atoms[i];
if (a.getGroup ().getNitrogenAtom () === a && a.getCovalentBondCount () == 1) a.setFormalCharge (1);
if ((i = this.bsAtomsForHs.nextClearBit (i + 1)) < 0) break;
}
});
Clazz.defineMethod (c$, "finalizePdbMultipleBonds", 
 function () {
var htKeysUsed =  new java.util.Hashtable ();
var bondCount = this.ms.bondCount;
var bonds = this.ms.bo;
for (var i = this.baseBondIndex; i < bondCount; i++) {
var a1 = bonds[i].getAtom1 ();
var a2 = bonds[i].getAtom2 ();
var g = a1.getGroup ();
if (g !== a2.getGroup ()) continue;
var key =  new JU.SB ().append (g.getGroup3 ());
key.append (":");
var n1 = a1.getAtomName ();
var n2 = a2.getAtomName ();
if (n1.compareTo (n2) > 0) key.append (n2).append (":").append (n1);
 else key.append (n1).append (":").append (n2);
var skey = key.toString ();
var type = this.htBondMap.get (skey);
if (type == null) continue;
htKeysUsed.put (skey, Boolean.TRUE);
bonds[i].setOrder (JU.PT.parseInt (type));
}
for (var key, $key = this.htBondMap.keySet ().iterator (); $key.hasNext () && ((key = $key.next ()) || true);) {
if (htKeysUsed.get (key) != null) continue;
if (key.indexOf (":") < 0) {
htKeysUsed.put (key, Boolean.TRUE);
continue;
}var value = this.htBondMap.get (key);
JU.Logger.info ("bond " + key + " was not used; order=" + value);
if (this.htBondMap.get (key).equals ("1")) {
htKeysUsed.put (key, Boolean.TRUE);
continue;
}}
var htKeysBad =  new java.util.Hashtable ();
for (var key, $key = this.htBondMap.keySet ().iterator (); $key.hasNext () && ((key = $key.next ()) || true);) {
if (htKeysUsed.get (key) != null) continue;
htKeysBad.put (key.substring (0, key.lastIndexOf (":")), this.htBondMap.get (key));
}
if (htKeysBad.isEmpty ()) return;
for (var i = 0; i < bondCount; i++) {
var a1 = bonds[i].getAtom1 ();
var a2 = bonds[i].getAtom2 ();
if (a1.getGroup () === a2.getGroup ()) continue;
var value;
if ((value = htKeysBad.get (a1.getGroup3 (false) + ":" + a1.getAtomName ())) == null && ((value = htKeysBad.get (a2.getGroup3 (false) + ":" + a2.getAtomName ())) == null)) continue;
bonds[i].setOrder (JU.PT.parseInt (value));
JU.Logger.info ("assigning order " + bonds[i].order + " to bond " + bonds[i]);
}
});
Clazz.defineMethod (c$, "setHydrogen", 
 function (iTo, iAtom, name, pt) {
if (!this.bsAddedHydrogens.get (iAtom)) return;
var atoms = this.ms.at;
if (this.lastSetH == -2147483648 || atoms[iAtom].mi != atoms[this.lastSetH].mi) this.maxSerial = (this.ms.getInfo (atoms[this.lastSetH = iAtom].mi, "PDB_CONECT_firstAtom_count_max"))[2];
this.bsAddedHydrogens.clear (iAtom);
this.ms.setAtomName (iAtom, name);
atoms[iAtom].setT (pt);
this.ms.setAtomNumber (iAtom, ++this.maxSerial);
atoms[iAtom].setAtomSymmetry (atoms[iTo].getAtomSymmetry ());
this.ml.undeleteAtom (iAtom);
this.ms.bondAtoms (atoms[iTo], atoms[iAtom], 1, this.ms.getDefaultMadFromOrder (1), null, 0, true, false);
}, "~N,~N,~S,JU.P3");
Clazz.overrideMethod (c$, "fixPropertyValue", 
function (bsAtoms, data) {
var aData = JU.PT.split (data, "\n");
var atoms = this.ms.at;
var newData =  new Array (bsAtoms.cardinality ());
var lastData = "";
for (var pt = 0, iAtom = 0, i = bsAtoms.nextSetBit (0); i >= 0; i = bsAtoms.nextSetBit (i), iAtom++) {
if (atoms[i].getElementNumber () != 1) lastData = aData[pt++];
newData[iAtom] = lastData;
}
return JU.PT.join (newData, '\n', 0);
}, "JU.BS,~S");
c$.allocateBioPolymer = Clazz.defineMethod (c$, "allocateBioPolymer", 
function (groups, firstGroupIndex, checkConnections) {
var previous = null;
var count = 0;
for (var i = firstGroupIndex; i < groups.length; ++i) {
var group = groups[i];
var current;
if (!(Clazz.instanceOf (group, JM.Monomer)) || (current = group).bioPolymer != null || previous != null && previous.getClass () !== current.getClass () || checkConnections && !current.isConnectedAfter (previous)) break;
previous = current;
count++;
}
if (count == 0) return null;
var monomers =  new Array (count);
for (var j = 0; j < count; ++j) monomers[j] = groups[firstGroupIndex + j];

if (Clazz.instanceOf (previous, JM.AminoMonomer)) return  new JM.AminoPolymer (monomers);
if (Clazz.instanceOf (previous, JM.AlphaMonomer)) return  new JM.AlphaPolymer (monomers);
if (Clazz.instanceOf (previous, JM.NucleicMonomer)) return  new JM.NucleicPolymer (monomers);
if (Clazz.instanceOf (previous, JM.PhosphorusMonomer)) return  new JM.PhosphorusPolymer (monomers);
if (Clazz.instanceOf (previous, JM.CarbohydrateMonomer)) return  new JM.CarbohydratePolymer (monomers);
JU.Logger.error ("Polymer.allocatePolymer() ... no matching polymer for monomor " + previous);
throw  new NullPointerException ();
}, "~A,~N,~B");
Clazz.overrideMethod (c$, "iterateOverAllNewStructures", 
function (adapter, atomSetCollection) {
var iterStructure = adapter.getStructureIterator (atomSetCollection);
if (iterStructure == null) return;
var bs = iterStructure.getStructuredModels ();
if (bs != null) for (var i = bs.nextSetBit (0); i >= 0; i = bs.nextSetBit (i + 1)) this.ml.structuresDefinedInFile.set (this.ml.baseModelIndex + i);

while (iterStructure.hasNext ()) if (iterStructure.getStructureType () !== J.c.STR.TURN) this.setStructure (iterStructure);

iterStructure = adapter.getStructureIterator (atomSetCollection);
while (iterStructure.hasNext ()) if (iterStructure.getStructureType () === J.c.STR.TURN) this.setStructure (iterStructure);

}, "J.api.JmolAdapter,~O");
Clazz.defineMethod (c$, "setStructure", 
 function (iterStructure) {
var t = iterStructure.getSubstructureType ();
var id = iterStructure.getStructureID ();
var serID = iterStructure.getSerialID ();
var count = iterStructure.getStrandCount ();
var atomRange = iterStructure.getAtomIndices ();
var modelRange = iterStructure.getModelIndices ();
if (this.bsAssigned == null) this.bsAssigned =  new JU.BS ();
this.defineStructure (t, id, serID, count, iterStructure.getStartChainID (), iterStructure.getStartSequenceNumber (), iterStructure.getStartInsertionCode (), iterStructure.getEndChainID (), iterStructure.getEndSequenceNumber (), iterStructure.getEndInsertionCode (), atomRange, modelRange, this.bsAssigned);
}, "J.api.JmolAdapterStructureIterator");
Clazz.defineMethod (c$, "defineStructure", 
 function (subType, structureID, serialID, strandCount, startChainID, startSequenceNumber, startInsertionCode, endChainID, endSequenceNumber, endInsertionCode, atomRange, modelRange, bsAssigned) {
var type = (subType === J.c.STR.NOT ? J.c.STR.NONE : subType);
var startSeqCode = JM.Group.getSeqcodeFor (startSequenceNumber, startInsertionCode);
var endSeqCode = JM.Group.getSeqcodeFor (endSequenceNumber, endInsertionCode);
var models = this.ms.am;
if (this.ml.isTrajectory) {
modelRange[1] = modelRange[0];
} else {
modelRange[0] += this.ml.baseModelIndex;
modelRange[1] += this.ml.baseModelIndex;
}this.ml.structuresDefinedInFile.setBits (modelRange[0], modelRange[1] + 1);
for (var i = modelRange[0]; i <= modelRange[1]; i++) {
var i0 = models[i].firstAtomIndex;
if (Clazz.instanceOf (models[i], JM.BioModel)) (models[i]).addSecondaryStructure (type, structureID, serialID, strandCount, startChainID, startSeqCode, endChainID, endSeqCode, i0 + atomRange[0], i0 + atomRange[1], bsAssigned);
}
}, "J.c.STR,~S,~N,~N,~N,~N,~S,~N,~N,~S,~A,~A,JU.BS");
c$.$Resolver$BondSorter$ = function () {
Clazz.pu$h(self.c$);
c$ = Clazz.decorateAsClass (function () {
Clazz.prepareCallback (this, arguments);
Clazz.instantialize (this, arguments);
}, JM.Resolver, "BondSorter", null, java.util.Comparator);
Clazz.overrideMethod (c$, "compare", 
function (a, b) {
return (b == null ? (a == null ? 0 : -1) : a == null ? 1 : a[0].compareTo (b[0]) < 0 ? -1 : a[0].compareTo (b[0]) > 0 ? 1 : a[3].compareTo (b[3]) < 0 ? -1 : a[3].compareTo (b[3]) > 0 ? 1 : a[1].compareTo (b[1]) < 0 ? -1 : a[1].compareTo (b[1]) > 0 ? 1 : 0);
}, "~A,~A");
c$ = Clazz.p0p ();
};
});
