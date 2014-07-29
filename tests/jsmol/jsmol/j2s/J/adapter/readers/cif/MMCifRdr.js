Clazz.declarePackage ("J.adapter.readers.cif");
Clazz.load (null, "J.adapter.readers.cif.MMCifRdr", ["java.util.Hashtable", "JU.BS", "$.Lst", "$.M4", "$.P3", "$.PT", "$.SB", "J.adapter.smarter.Atom", "$.Structure", "J.api.JmolAdapter", "J.c.STR", "JU.Logger"], function () {
c$ = Clazz.decorateAsClass (function () {
this.cr = null;
this.isBiomolecule = false;
this.byChain = false;
this.bySymop = false;
this.isCourseGrained = false;
this.chainAtomMap = null;
this.chainAtomCounts = null;
this.vBiomolecules = null;
this.thisBiomolecule = null;
this.htBiomts = null;
this.htSites = null;
this.assemblyIdAtoms = null;
this.thisChain = -1;
this.chainSum = null;
this.chainAtomCount = null;
this.assem = null;
this.hetatmData = null;
this.field = null;
this.firstChar = '\0';
this.htHetero = null;
this.propertyCount = 0;
this.fieldOf = null;
Clazz.instantialize (this, arguments);
}, J.adapter.readers.cif, "MMCifRdr");
Clazz.makeConstructor (c$, 
function () {
});
Clazz.defineMethod (c$, "initialize", 
function (r) {
this.cr = r;
this.byChain = r.checkFilterKey ("BYCHAIN");
this.bySymop = r.checkFilterKey ("BYSYMOP");
this.isCourseGrained = this.byChain || this.bySymop;
if (this.byChain) {
this.chainAtomMap =  new java.util.Hashtable ();
this.chainAtomCounts =  new java.util.Hashtable ();
}if (this.cr.checkFilterKey ("BIOMOLECULE")) this.cr.filter = JU.PT.rep (this.cr.filter, "BIOMOLECULE", "ASSEMBLY");
this.isBiomolecule = this.cr.checkFilterKey ("ASSEMBLY");
return this.isCourseGrained;
}, "J.adapter.readers.cif.CifReader");
Clazz.defineMethod (c$, "finalizeReader", 
function (nAtoms) {
if (this.byChain && !this.isBiomolecule) for (var id, $id = this.chainAtomMap.keySet ().iterator (); $id.hasNext () && ((id = $id.next ()) || true);) this.createParticle (id);

var asc = this.cr.asc;
if (!this.isCourseGrained && asc.ac == nAtoms) asc.removeCurrentAtomSet ();
 else this.cr.applySymmetryAndSetTrajectory ();
if (this.htSites != null) this.cr.addSites (this.htSites);
if (this.vBiomolecules != null && this.vBiomolecules.size () == 1 && (this.isCourseGrained || asc.ac > 0)) {
asc.setAtomSetAuxiliaryInfo ("biomolecules", this.vBiomolecules);
var ht = this.vBiomolecules.get (0);
this.cr.appendLoadNote ("Constructing " + ht.get ("name"));
this.setBiomolecules (ht);
if (this.thisBiomolecule != null) {
asc.getXSymmetry ().applySymmetryBio (this.thisBiomolecule, this.cr.notionalUnitCell, this.cr.applySymmetryToBonds, this.cr.filter);
asc.xtalSymmetry = null;
}}}, "~N");
Clazz.defineMethod (c$, "processEntry", 
function () {
if (this.cr.key.startsWith ("_pdbx_entity_nonpoly")) this.processDataNonpoly ();
 else if (this.cr.key.startsWith ("_pdbx_struct_assembly_gen")) this.processDataAssemblyGen ();
});
Clazz.defineMethod (c$, "processSequence", 
 function () {
this.parseLoopParameters (J.adapter.readers.cif.MMCifRdr.structRefFields);
while (this.cr.parser.getData ()) {
var g1 = null;
var g3 = null;
var n = this.cr.parser.getFieldCount ();
for (var i = 0; i < n; ++i) {
switch (this.fieldProperty (i)) {
case 0:
g3 = this.field;
break;
case 1:
if (this.field.length == 1) g1 = this.field.toLowerCase ();
}
}
if (g1 != null && g3 != null) {
if (this.cr.htGroup1 == null) this.cr.asc.setInfo ("htGroup1", this.cr.htGroup1 =  new java.util.Hashtable ());
this.cr.htGroup1.put (g3, g1);
}}
return true;
});
Clazz.defineMethod (c$, "processDataNonpoly", 
 function () {
if (this.hetatmData == null) this.hetatmData =  new Array (3);
for (var i = J.adapter.readers.cif.MMCifRdr.nonpolyFields.length; --i >= 0; ) if (this.cr.key.equals (J.adapter.readers.cif.MMCifRdr.nonpolyFields[i])) {
this.hetatmData[i] = this.cr.data;
break;
}
if (this.hetatmData[1] == null || this.hetatmData[2] == null) return;
this.addHetero (this.hetatmData[2], this.hetatmData[1]);
this.hetatmData = null;
});
Clazz.defineMethod (c$, "processDataAssemblyGen", 
 function () {
if (this.assem == null) this.assem =  new Array (3);
if (this.cr.key.indexOf ("assembly_id") >= 0) this.assem[0] = this.cr.parser.fullTrim (this.cr.data);
 else if (this.cr.key.indexOf ("oper_expression") >= 0) this.assem[1] = this.cr.parser.fullTrim (this.cr.data);
 else if (this.cr.key.indexOf ("asym_id_list") >= 0) this.assem[2] = this.cr.parser.fullTrim (this.cr.data);
if (this.assem[0] != null && this.assem[1] != null && this.assem[2] != null) this.addAssembly ();
});
Clazz.defineMethod (c$, "processAssemblyGenBlock", 
 function () {
this.parseLoopParameters (J.adapter.readers.cif.MMCifRdr.assemblyFields);
while (this.cr.parser.getData ()) {
this.assem =  new Array (3);
var count = 0;
var p;
var n = this.cr.parser.getFieldCount ();
for (var i = 0; i < n; ++i) {
switch (p = this.fieldProperty (i)) {
case 0:
case 1:
case 2:
count++;
this.assem[p] = this.field;
break;
}
}
if (count == 3) this.addAssembly ();
}
this.assem = null;
return true;
});
Clazz.defineMethod (c$, "addAssembly", 
 function () {
var id = this.assem[0];
var iMolecule = this.cr.parseIntStr (id);
var list = this.assem[2];
this.cr.appendLoadNote ("found biomolecule " + id + ": " + list);
if (!this.cr.checkFilterKey ("ASSEMBLY " + id + ";")) return;
if (this.vBiomolecules == null) {
this.vBiomolecules =  new JU.Lst ();
}var info =  new java.util.Hashtable ();
info.put ("name", "biomolecule " + id);
info.put ("molecule", iMolecule == -2147483648 ? id : Integer.$valueOf (iMolecule));
info.put ("assemblies", "$" + list.$replace (',', '$'));
info.put ("operators", this.decodeAssemblyOperators (this.assem[1]));
info.put ("biomts",  new JU.Lst ());
this.thisBiomolecule = info;
JU.Logger.info ("assembly " + id + " operators " + this.assem[1] + " ASYM_IDs " + this.assem[2]);
this.vBiomolecules.addLast (info);
this.assem = null;
});
Clazz.defineMethod (c$, "decodeAssemblyOperators", 
 function (ops) {
var pt = ops.indexOf (")(");
if (pt >= 0) return this.crossBinary (this.decodeAssemblyOperators (ops.substring (0, pt + 1)), this.decodeAssemblyOperators (ops.substring (pt + 1)));
if (ops.startsWith ("(")) {
if (ops.indexOf ("-") >= 0) ops = JU.BS.unescape ("({" + ops.substring (1, ops.length - 1).$replace ('-', ':') + "})").toString ();
ops = JU.PT.rep (ops, " ", "");
ops = ops.substring (1, ops.length - 1);
}return ops;
}, "~S");
Clazz.defineMethod (c$, "crossBinary", 
 function (ops1, ops2) {
var sb =  new JU.SB ();
var opsLeft = JU.PT.split (ops1, ",");
var opsRight = JU.PT.split (ops2, ",");
for (var i = 0; i < opsLeft.length; i++) for (var j = 0; j < opsRight.length; j++) sb.append (",").append (opsLeft[i]).append ("|").append (opsRight[j]);


return sb.toString ().substring (1);
}, "~S,~S");
Clazz.defineMethod (c$, "processStructOperListBlock", 
 function () {
this.parseLoopParameters (J.adapter.readers.cif.MMCifRdr.operFields);
var m =  Clazz.newFloatArray (16, 0);
m[15] = 1;
while (this.cr.parser.getData ()) {
var count = 0;
var id = null;
var xyz = null;
var n = this.cr.parser.getFieldCount ();
for (var i = 0; i < n; ++i) {
var p = this.fieldProperty (i);
switch (p) {
case -1:
break;
case 12:
id = this.field;
break;
case 13:
xyz = this.field;
break;
default:
m[p] = this.cr.parseFloatStr (this.field);
++count;
}
}
if (id != null && (count == 12 || xyz != null && this.cr.symmetry != null)) {
JU.Logger.info ("assembly operator " + id + " " + xyz);
var m4 =  new JU.M4 ();
if (count != 12) {
this.cr.symmetry.getMatrixFromString (xyz, m, false, 0);
m[3] *= this.cr.symmetry.getUnitCellInfoType (0) / 12;
m[7] *= this.cr.symmetry.getUnitCellInfoType (1) / 12;
m[11] *= this.cr.symmetry.getUnitCellInfoType (2) / 12;
}m4.setA (m);
if (this.htBiomts == null) this.htBiomts =  new java.util.Hashtable ();
this.htBiomts.put (id, m4);
}}
return true;
});
Clazz.defineMethod (c$, "processChemCompLoopBlock", 
 function () {
this.parseLoopParameters (J.adapter.readers.cif.MMCifRdr.chemCompFields);
while (this.cr.parser.getData ()) {
var groupName = null;
var hetName = null;
var n = this.cr.parser.getFieldCount ();
for (var i = 0; i < n; ++i) {
switch (this.fieldProperty (i)) {
case -1:
break;
case 0:
groupName = this.field;
break;
case 1:
hetName = this.field;
break;
}
}
if (groupName != null && hetName != null) this.addHetero (groupName, hetName);
}
return true;
});
Clazz.defineMethod (c$, "processNonpolyLoopBlock", 
 function () {
this.parseLoopParameters (J.adapter.readers.cif.MMCifRdr.nonpolyFields);
while (this.cr.parser.getData ()) {
var groupName = null;
var hetName = null;
var n = this.cr.parser.getFieldCount ();
for (var i = 0; i < n; ++i) {
switch (this.fieldProperty (i)) {
case -1:
case 0:
break;
case 2:
groupName = this.field;
break;
case 1:
hetName = this.field;
break;
}
}
if (groupName == null || hetName == null) return false;
this.addHetero (groupName, hetName);
}
return true;
});
Clazz.defineMethod (c$, "addHetero", 
 function (groupName, hetName) {
if (!J.api.JmolAdapter.isHetero (groupName)) return;
if (this.htHetero == null) this.htHetero =  new java.util.Hashtable ();
this.htHetero.put (groupName, hetName);
if (JU.Logger.debugging) {
JU.Logger.debug ("hetero: " + groupName + " = " + hetName);
}}, "~S,~S");
Clazz.defineMethod (c$, "processStructConfLoopBlock", 
 function () {
this.parseLoopParameters (J.adapter.readers.cif.MMCifRdr.structConfFields);
for (var i = this.propertyCount; --i >= 0; ) if (this.fieldOf[i] == -1) {
JU.Logger.warn ("?que? missing property: " + J.adapter.readers.cif.MMCifRdr.structConfFields[i]);
return false;
}
while (this.cr.parser.getData ()) {
var structure =  new J.adapter.smarter.Structure (-1, J.c.STR.HELIX, J.c.STR.HELIX, null, 0, 0);
var n = this.cr.parser.getFieldCount ();
for (var i = 0; i < n; ++i) {
switch (this.fieldProperty (i)) {
case -1:
break;
case 0:
if (this.field.startsWith ("TURN")) structure.structureType = structure.substructureType = J.c.STR.TURN;
 else if (!this.field.startsWith ("HELX")) structure.structureType = structure.substructureType = J.c.STR.NONE;
break;
case 1:
structure.startChainStr = this.field;
structure.startChainID = this.cr.vwr.getChainID (this.field);
break;
case 2:
structure.startSequenceNumber = this.cr.parseIntStr (this.field);
break;
case 3:
structure.startInsertionCode = this.firstChar;
break;
case 4:
structure.endChainStr = this.field;
structure.endChainID = this.cr.vwr.getChainID (this.field);
break;
case 5:
structure.endSequenceNumber = this.cr.parseIntStr (this.field);
break;
case 9:
structure.substructureType = J.adapter.smarter.Structure.getHelixType (this.cr.parseIntStr (this.field));
break;
case 6:
structure.endInsertionCode = this.firstChar;
break;
case 7:
structure.structureID = this.field;
break;
case 8:
structure.serialID = this.cr.parseIntStr (this.field);
break;
}
}
this.cr.asc.addStructure (structure);
}
return true;
});
Clazz.defineMethod (c$, "processStructSheetRangeLoopBlock", 
 function () {
this.parseLoopParameters (J.adapter.readers.cif.MMCifRdr.structSheetRangeFields);
for (var i = this.propertyCount; --i >= 0; ) if (this.fieldOf[i] == -1) {
JU.Logger.warn ("?que? missing property:" + J.adapter.readers.cif.MMCifRdr.structSheetRangeFields[i]);
return false;
}
while (this.cr.parser.getData ()) {
var structure =  new J.adapter.smarter.Structure (-1, J.c.STR.SHEET, J.c.STR.SHEET, null, 0, 0);
var n = this.cr.parser.getFieldCount ();
for (var i = 0; i < n; ++i) {
switch (this.fieldProperty (i)) {
case 1:
structure.startChainID = this.cr.vwr.getChainID (this.field);
break;
case 2:
structure.startSequenceNumber = this.cr.parseIntStr (this.field);
break;
case 3:
structure.startInsertionCode = this.firstChar;
break;
case 4:
structure.endChainID = this.cr.vwr.getChainID (this.field);
break;
case 5:
structure.endSequenceNumber = this.cr.parseIntStr (this.field);
break;
case 6:
structure.endInsertionCode = this.firstChar;
break;
case 0:
structure.strandCount = 1;
structure.structureID = this.field;
break;
case 7:
structure.serialID = this.cr.parseIntStr (this.field);
break;
}
}
this.cr.asc.addStructure (structure);
}
return true;
});
Clazz.defineMethod (c$, "parseLoopParameters", 
 function (fields) {
this.cr.parseLoopParameters (fields);
this.propertyCount = fields.length;
this.fieldOf = this.cr.fieldOf;
}, "~A");
Clazz.defineMethod (c$, "processStructSiteBlock", 
 function () {
this.parseLoopParameters (J.adapter.readers.cif.MMCifRdr.structSiteRangeFields);
for (var i = 3; --i >= 0; ) if (this.fieldOf[i] == -1) {
JU.Logger.warn ("?que? missing property: " + J.adapter.readers.cif.MMCifRdr.structSiteRangeFields[i]);
return false;
}
var siteID = "";
var seqNum = "";
var insCode = "";
var chainID = "";
var resID = "";
var group = "";
var htSite = null;
this.htSites =  new java.util.Hashtable ();
while (this.cr.parser.getData ()) {
var n = this.cr.parser.getFieldCount ();
for (var i = 0; i < n; ++i) {
switch (this.fieldProperty (i)) {
case 0:
if (group !== "") {
var groups = htSite.get ("groups");
groups += (groups.length == 0 ? "" : ",") + group;
group = "";
htSite.put ("groups", groups);
}siteID = this.field;
htSite = this.htSites.get (siteID);
if (htSite == null) {
htSite =  new java.util.Hashtable ();
htSite.put ("groups", "");
this.htSites.put (siteID, htSite);
}seqNum = "";
insCode = "";
chainID = "";
resID = "";
break;
case 1:
resID = this.field;
break;
case 2:
chainID = this.field;
break;
case 3:
seqNum = this.field;
break;
case 4:
insCode = this.field;
break;
}
if (seqNum !== "" && resID !== "") group = "[" + resID + "]" + seqNum + (insCode.length > 0 ? "^" + insCode : "") + (chainID.length > 0 ? ":" + chainID : "");
}
}
if (group !== "") {
var groups = htSite.get ("groups");
groups += (groups.length == 0 ? "" : ",") + group;
group = "";
htSite.put ("groups", groups);
}return true;
});
Clazz.defineMethod (c$, "fieldProperty", 
 function (i) {
return ((this.field = this.cr.parser.getLoopData (i)).length > 0 && (this.firstChar = this.field.charAt (0)) != '\0' ? this.cr.propertyOf[i] : -1);
}, "~N");
Clazz.defineMethod (c$, "setBiomolecules", 
 function (biomolecule) {
if (!this.isBiomolecule || this.assemblyIdAtoms == null && this.chainAtomCounts == null) return;
var mident = JU.M4.newM4 (null);
var ops = JU.PT.split (biomolecule.get ("operators"), ",");
var assemblies = biomolecule.get ("assemblies");
var biomts =  new JU.Lst ();
biomolecule.put ("biomts", biomts);
biomts.addLast (mident);
for (var j = 0; j < ops.length; j++) {
var m = this.getOpMatrix (ops[j]);
if (m != null && !m.equals (mident)) biomts.addLast (m);
}
var bsAll =  new JU.BS ();
var sum =  new JU.P3 ();
var count = 0;
var nAtoms = 0;
var ids = JU.PT.split (assemblies, "$");
for (var j = 1; j < ids.length; j++) {
var id = ids[j];
if (this.assemblyIdAtoms != null) {
var bs = this.assemblyIdAtoms.get (id);
if (bs != null) {
bsAll.or (bs);
}} else if (this.isCourseGrained) {
var asum = this.chainAtomMap.get (id);
var c = this.chainAtomCounts.get (id)[0];
if (asum != null) {
if (this.bySymop) {
sum.add (asum);
count += c;
} else {
this.createParticle (id);
nAtoms++;
}}}}
if (this.isCourseGrained) {
if (this.bySymop) {
nAtoms = 1;
var a1 =  new J.adapter.smarter.Atom ();
a1.setT (sum);
a1.scale (1 / count);
a1.radius = 16;
}} else {
nAtoms = bsAll.cardinality ();
if (nAtoms < this.cr.asc.ac) this.cr.asc.bsAtoms = bsAll;
}biomolecule.put ("atomCount", Integer.$valueOf (nAtoms * ops.length));
}, "java.util.Map");
Clazz.defineMethod (c$, "createParticle", 
 function (id) {
var asum = this.chainAtomMap.get (id);
var c = this.chainAtomCounts.get (id)[0];
var a =  new J.adapter.smarter.Atom ();
a.setT (asum);
a.scale (1 / c);
a.elementSymbol = "Pt";
a.chainID = this.cr.vwr.getChainID (id);
a.radius = 16;
this.cr.asc.addAtom (a);
}, "~S");
Clazz.defineMethod (c$, "getOpMatrix", 
 function (ops) {
if (this.htBiomts == null) return JU.M4.newM4 (null);
var pt = ops.indexOf ("|");
if (pt >= 0) {
var m = JU.M4.newM4 (this.htBiomts.get (ops.substring (0, pt)));
m.mul (this.htBiomts.get (ops.substring (pt + 1)));
return m;
}return this.htBiomts.get (ops);
}, "~S");
Clazz.defineMethod (c$, "processLigandBondLoopBlock", 
 function () {
this.parseLoopParameters (J.adapter.readers.cif.MMCifRdr.chemCompBondFields);
for (var i = this.propertyCount; --i >= 0; ) if (this.fieldOf[i] == -1) {
JU.Logger.warn ("?que? missing property: " + J.adapter.readers.cif.MMCifRdr.chemCompBondFields[i]);
return false;
}
var order = 0;
var isAromatic = false;
while (this.cr.parser.getData ()) {
var atom1 = null;
var atom2 = null;
order = 0;
isAromatic = false;
var n = this.cr.parser.getFieldCount ();
for (var i = 0; i < n; ++i) {
switch (this.fieldProperty (i)) {
case 0:
atom1 = this.cr.asc.getAtomFromName (this.field);
break;
case 1:
atom2 = this.cr.asc.getAtomFromName (this.field);
break;
case 3:
isAromatic = (this.field.charAt (0) == 'Y');
break;
case 2:
order = this.cr.getBondOrder (this.field);
break;
}
}
if (isAromatic) switch (order) {
case 1:
order = 513;
break;
case 2:
order = 514;
break;
}
this.cr.asc.addNewBondWithOrderA (atom1, atom2, order);
}
return true;
});
Clazz.defineMethod (c$, "checkAtom", 
function (atom, assemblyId, index) {
if (this.byChain && !this.isBiomolecule) {
if (this.thisChain != atom.chainID) {
this.thisChain = atom.chainID;
var id = "" + atom.chainID;
this.chainSum = this.chainAtomMap.get (id);
if (this.chainSum == null) {
this.chainAtomMap.put (id, this.chainSum =  new JU.P3 ());
this.chainAtomCounts.put (id, this.chainAtomCount =  Clazz.newIntArray (1, 0));
}}this.chainSum.add (atom);
this.chainAtomCount[0]++;
return false;
}if (this.isBiomolecule && this.isCourseGrained) {
var sum = this.chainAtomMap.get (assemblyId);
if (sum == null) {
this.chainAtomMap.put (assemblyId, sum =  new JU.P3 ());
this.chainAtomCounts.put (assemblyId,  Clazz.newIntArray (1, 0));
}this.chainAtomCounts.get (assemblyId)[0]++;
sum.add (atom);
return false;
}if (assemblyId != null) {
if (this.assemblyIdAtoms == null) this.assemblyIdAtoms =  new java.util.Hashtable ();
var bs = this.assemblyIdAtoms.get (assemblyId);
if (bs == null) this.assemblyIdAtoms.put (assemblyId, bs =  new JU.BS ());
bs.set (index);
}if (atom.isHetero && this.htHetero != null) {
this.cr.asc.setAtomSetAuxiliaryInfo ("hetNames", this.htHetero);
this.cr.asc.setInfo ("hetNames", this.htHetero);
this.htHetero = null;
}return true;
}, "J.adapter.smarter.Atom,~S,~N");
Clazz.defineMethod (c$, "processLoopBlock", 
function () {
var key = this.cr.key;
if (key.startsWith ("_pdbx_struct_oper_list")) return this.processStructOperListBlock ();
if (key.startsWith ("_pdbx_struct_assembly_gen")) return this.processAssemblyGenBlock ();
if (key.startsWith ("_struct_ref_seq_dif")) return this.processSequence ();
if (this.isCourseGrained) return false;
if (key.startsWith ("_struct_site_gen")) return this.processStructSiteBlock ();
if (key.startsWith ("_chem_comp_bond")) return this.processLigandBondLoopBlock ();
if (key.startsWith ("_chem_comp")) return this.processChemCompLoopBlock ();
if (key.startsWith ("_pdbx_entity_nonpoly")) return this.processNonpolyLoopBlock ();
if (key.startsWith ("_struct_conf") && !key.startsWith ("_struct_conf_type")) return this.processStructConfLoopBlock ();
if (key.startsWith ("_struct_sheet_range")) return this.processStructSheetRangeLoopBlock ();
return false;
});
Clazz.defineStatics (c$,
"NONE", -1,
"OPER_ID", 12,
"OPER_XYZ", 13,
"operFields", ["_pdbx_struct_oper_list_matrix[1][1]", "_pdbx_struct_oper_list_matrix[1][2]", "_pdbx_struct_oper_list_matrix[1][3]", "_pdbx_struct_oper_list_vector[1]", "_pdbx_struct_oper_list_matrix[2][1]", "_pdbx_struct_oper_list_matrix[2][2]", "_pdbx_struct_oper_list_matrix[2][3]", "_pdbx_struct_oper_list_vector[2]", "_pdbx_struct_oper_list_matrix[3][1]", "_pdbx_struct_oper_list_matrix[3][2]", "_pdbx_struct_oper_list_matrix[3][3]", "_pdbx_struct_oper_list_vector[3]", "_pdbx_struct_oper_list_id", "_pdbx_struct_oper_list_symmetry_operation"],
"ASSEM_ID", 0,
"ASSEM_OPERS", 1,
"ASSEM_LIST", 2,
"assemblyFields", ["_pdbx_struct_assembly_gen_assembly_id", "_pdbx_struct_assembly_gen_oper_expression", "_pdbx_struct_assembly_gen_asym_id_list"],
"STRUCT_REF_G3", 0,
"STRUCT_REF_G1", 1,
"structRefFields", ["_struct_ref_seq_dif_mon_id", "_struct_ref_seq_dif.db_mon_id"],
"NONPOLY_ENTITY_ID", 0,
"NONPOLY_NAME", 1,
"NONPOLY_COMP_ID", 2,
"nonpolyFields", ["_pdbx_entity_nonpoly_entity_id", "_pdbx_entity_nonpoly_name", "_pdbx_entity_nonpoly_comp_id"],
"CHEM_COMP_ID", 0,
"CHEM_COMP_NAME", 1,
"chemCompFields", ["_chem_comp_id", "_chem_comp_name"],
"CONF_TYPE_ID", 0,
"BEG_ASYM_ID", 1,
"BEG_SEQ_ID", 2,
"BEG_INS_CODE", 3,
"END_ASYM_ID", 4,
"END_SEQ_ID", 5,
"END_INS_CODE", 6,
"STRUCT_ID", 7,
"SERIAL_NO", 8,
"HELIX_CLASS", 9,
"structConfFields", ["_struct_conf_conf_type_id", "_struct_conf_beg_auth_asym_id", "_struct_conf_beg_auth_seq_id", "_struct_conf_pdbx_beg_pdb_ins_code", "_struct_conf_end_auth_asym_id", "_struct_conf_end_auth_seq_id", "_struct_conf_pdbx_end_pdb_ins_code", "_struct_conf_id", "_struct_conf_pdbx_pdb_helix_id", "_struct_conf_pdbx_pdb_helix_class"],
"SHEET_ID", 0,
"STRAND_ID", 7,
"structSheetRangeFields", ["_struct_sheet_range_sheet_id", "_struct_sheet_range_beg_auth_asym_id", "_struct_sheet_range_beg_auth_seq_id", "_struct_sheet_range_pdbx_beg_pdb_ins_code", "_struct_sheet_range_end_auth_asym_id", "_struct_sheet_range_end_auth_seq_id", "_struct_sheet_range_pdbx_end_pdb_ins_code", "_struct_sheet_range_id"],
"SITE_ID", 0,
"SITE_COMP_ID", 1,
"SITE_ASYM_ID", 2,
"SITE_SEQ_ID", 3,
"SITE_INS_CODE", 4,
"structSiteRangeFields", ["_struct_site_gen_site_id", "_struct_site_gen_auth_comp_id", "_struct_site_gen_auth_asym_id", "_struct_site_gen_auth_seq_id", "_struct_site_gen_label_alt_id"],
"CHEM_COMP_BOND_ATOM_ID_1", 0,
"CHEM_COMP_BOND_ATOM_ID_2", 1,
"CHEM_COMP_BOND_VALUE_ORDER", 2,
"CHEM_COMP_BOND_AROMATIC_FLAG", 3,
"chemCompBondFields", ["_chem_comp_bond_atom_id_1", "_chem_comp_bond_atom_id_2", "_chem_comp_bond_value_order", "_chem_comp_bond_pdbx_aromatic_flag"]);
});
