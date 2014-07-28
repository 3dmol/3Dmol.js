Clazz.declarePackage ("JV");
Clazz.load (["J.api.JmolDataManager", "java.util.Hashtable"], "JV.DataManager", ["java.lang.Boolean", "JU.AU", "$.BS", "$.PT", "$.SB", "J.c.VDW", "JS.T", "JU.BSUtil", "$.Elements", "$.Escape", "$.Logger", "$.Parser"], function () {
c$ = Clazz.decorateAsClass (function () {
this.dataValues = null;
this.vwr = null;
Clazz.instantialize (this, arguments);
}, JV, "DataManager", null, J.api.JmolDataManager);
Clazz.prepareFields (c$, function () {
this.dataValues =  new java.util.Hashtable ();
});
Clazz.makeConstructor (c$, 
function () {
});
Clazz.overrideMethod (c$, "set", 
function (vwr) {
this.vwr = vwr;
return this;
}, "JV.Viewer");
Clazz.overrideMethod (c$, "clear", 
function () {
this.dataValues.clear ();
});
Clazz.overrideMethod (c$, "setData", 
function (type, data, arrayCount, actualAtomCount, matchField, matchFieldColumnCount, field, fieldColumnCount) {
if (type == null) {
this.clear ();
return;
}type = type.toLowerCase ();
if (type.equals ("element_vdw")) {
var stringData = (data[1]).trim ();
if (stringData.length == 0) {
this.vwr.userVdwMars = null;
this.vwr.userVdws = null;
this.vwr.bsUserVdws = null;
return;
}if (this.vwr.bsUserVdws == null) this.vwr.setUserVdw (this.vwr.defaultVdw);
JU.Parser.parseFloatArrayFromMatchAndField (stringData, this.vwr.bsUserVdws, 1, 0, data[2], 2, 0, this.vwr.userVdws, 1);
for (var i = this.vwr.userVdws.length; --i >= 0; ) this.vwr.userVdwMars[i] = Clazz.doubleToInt (Math.floor (this.vwr.userVdws[i] * 1000));

return;
}if (data[2] != null && arrayCount > 0) {
var createNew = (matchField != 0 || field != -2147483648 && field != 2147483647);
var oldData = this.dataValues.get (type);
var bs;
var f = (oldData == null || createNew ?  Clazz.newFloatArray (actualAtomCount, 0) : JU.AU.ensureLengthA ((oldData[1]), actualAtomCount));
var depth = (data[3]).intValue ();
var stringData = (depth == 0 ? data[1] : null);
var floatData = (depth == 1 ? data[1] : null);
var strData = null;
if (field == -2147483648 && (strData = JU.PT.getTokens (stringData)).length > 1) field = 0;
if (field == -2147483648) {
bs = data[2];
JV.DataManager.setSelectedFloats (JU.PT.parseFloat (stringData), bs, f);
} else if (field == 0 || field == 2147483647) {
bs = data[2];
if (floatData != null) {
if (floatData.length == bs.cardinality ()) for (var i = bs.nextSetBit (0), pt = 0; i >= 0; i = bs.nextSetBit (i + 1), pt++) f[i] = floatData[pt];

 else for (var i = bs.nextSetBit (0); i >= 0; i = bs.nextSetBit (i + 1)) f[i] = floatData[i];

} else {
JU.Parser.parseFloatArrayBsData (strData == null ? JU.PT.getTokens (stringData) : strData, bs, f);
}} else if (matchField <= 0) {
bs = data[2];
JU.Parser.parseFloatArrayFromMatchAndField (stringData, bs, 0, 0, null, field, fieldColumnCount, f, 1);
} else {
var iData = data[2];
JU.Parser.parseFloatArrayFromMatchAndField (stringData, null, matchField, matchFieldColumnCount, iData, field, fieldColumnCount, f, 1);
bs =  new JU.BS ();
for (var i = iData.length; --i >= 0; ) if (iData[i] >= 0) bs.set (iData[i]);

}if (oldData != null && Clazz.instanceOf (oldData[2], JU.BS) && !createNew) bs.or ((oldData[2]));
data[3] = Integer.$valueOf (1);
data[2] = bs;
data[1] = f;
if (type.indexOf ("property_atom.") == 0) {
var tok = JS.T.getSettableTokFromString (type = type.substring (14));
if (tok == 0) {
JU.Logger.error ("Unknown atom property: " + type);
return;
}var nValues = bs.cardinality ();
var fValues =  Clazz.newFloatArray (nValues, 0);
for (var n = 0, i = bs.nextSetBit (0); n < nValues; i = bs.nextSetBit (i + 1)) fValues[n++] = f[i];

this.vwr.setAtomProperty (bs, tok, 0, 0, null, fValues, null);
return;
}}this.dataValues.put (type, data);
}, "~S,~A,~N,~N,~N,~N,~N,~N");
c$.setSelectedFloats = Clazz.defineMethod (c$, "setSelectedFloats", 
 function (f, bs, data) {
var isAll = (bs == null);
var i0 = (isAll ? 0 : bs.nextSetBit (0));
for (var i = i0; i >= 0 && i < data.length; i = (isAll ? i + 1 : bs.nextSetBit (i + 1))) data[i] = f;

}, "~N,JU.BS,~A");
Clazz.overrideMethod (c$, "getData", 
function (type) {
if (this.dataValues.size () == 0 || type == null) return null;
if (!type.equalsIgnoreCase ("types")) return this.dataValues.get (type);
var info =  new Array (2);
info[0] = "types";
info[1] = "";
var n = 0;
for (var name, $name = this.dataValues.keySet ().iterator (); $name.hasNext () && ((name = $name.next ()) || true);) info[1] += (n++ > 0 ? "\n" : "") + name;

return info;
}, "~S");
Clazz.overrideMethod (c$, "getDataFloatA", 
function (label) {
if (this.dataValues.size () == 0) return null;
var data = this.getData (label);
if (data == null || (data[3]).intValue () != 1) return null;
return data[1];
}, "~S");
Clazz.overrideMethod (c$, "getDataFloat", 
function (label, atomIndex) {
if (this.dataValues.size () > 0) {
var data = this.getData (label);
if (data != null && (data[3]).intValue () == 1) {
var f = data[1];
if (atomIndex < f.length) return f[atomIndex];
}}return NaN;
}, "~S,~N");
Clazz.overrideMethod (c$, "getDataFloat2D", 
function (label) {
if (this.dataValues.size () == 0) return null;
var data = this.getData (label);
if (data == null || (data[3]).intValue () != 2) return null;
return data[1];
}, "~S");
Clazz.overrideMethod (c$, "getDataFloat3D", 
function (label) {
if (this.dataValues.size () == 0) return null;
var data = this.getData (label);
if (data == null || (data[3]).intValue () != 3) return null;
return data[1];
}, "~S");
Clazz.overrideMethod (c$, "deleteModelAtoms", 
function (firstAtomIndex, nAtoms, bsDeleted) {
if (this.dataValues.size () == 0) return;
for (var name, $name = this.dataValues.keySet ().iterator (); $name.hasNext () && ((name = $name.next ()) || true);) {
if (name.indexOf ("property_") == 0) {
var obj = this.dataValues.get (name);
JU.BSUtil.deleteBits (obj[2], bsDeleted);
switch ((obj[3]).intValue ()) {
case 1:
obj[1] = JU.AU.deleteElements (obj[1], firstAtomIndex, nAtoms);
break;
case 2:
obj[1] = JU.AU.deleteElements (obj[1], firstAtomIndex, nAtoms);
break;
default:
break;
}
}}
}, "~N,~N,JU.BS");
Clazz.overrideMethod (c$, "getDefaultVdwNameOrData", 
function (type, bs) {
var sb =  new JU.SB ();
sb.append (type.getVdwLabel ()).append ("\n");
var isAll = (bs == null);
var i0 = (isAll ? 1 : bs.nextSetBit (0));
var i1 = (isAll ? JU.Elements.elementNumberMax : bs.length ());
for (var i = i0; i < i1 && i >= 0; i = (isAll ? i + 1 : bs.nextSetBit (i + 1))) sb.appendI (i).appendC ('\t').appendF (type === J.c.VDW.USER ? this.vwr.userVdws[i] : JU.Elements.getVanderwaalsMar (i, type) / 1000).appendC ('\t').append (JU.Elements.elementSymbolFromNumber (i)).appendC ('\n');

return (bs == null ? sb.toString () : "\n  DATA \"element_vdw\"\n" + sb.append ("  end \"element_vdw\";\n\n").toString ());
}, "J.c.VDW,JU.BS");
Clazz.overrideMethod (c$, "getDataState", 
function (sc, sb) {
if (this.dataValues.size () == 0) return false;
var haveData = false;
for (var name, $name = this.dataValues.keySet ().iterator (); $name.hasNext () && ((name = $name.next ()) || true);) {
if (name.indexOf ("property_") == 0) {
var obj = this.dataValues.get (name);
if (obj.length > 4 && obj[4] === Boolean.FALSE) continue;
haveData = true;
var data = obj[1];
if (data != null && (obj[3]).intValue () == 1) {
sc.getAtomicPropertyStateBuffer (sb, 14, obj[2], name, data);
sb.append ("\n");
} else {
sb.append ("\n").append (JU.Escape.encapsulateData (name, data, 0));
}} else if (name.indexOf ("data2d") == 0) {
var obj = this.dataValues.get (name);
var data = obj[1];
if (data != null && (obj[3]).intValue () == 2) {
haveData = true;
sb.append ("\n").append (JU.Escape.encapsulateData (name, data, 2));
}} else if (name.indexOf ("data3d") == 0) {
var obj = this.dataValues.get (name);
var data = obj[1];
if (data != null && (obj[3]).intValue () == 3) {
haveData = true;
sb.append ("\n").append (JU.Escape.encapsulateData (name, data, 3));
}}}
return haveData;
}, "JV.JmolStateCreator,JU.SB");
Clazz.defineStatics (c$,
"DATA_TYPE_STRING", 0,
"DATA_TYPE_AF", 1,
"DATA_ARRAY_FF", 2,
"DATA_ARRAY_FFF", 3,
"DATA_VALUE", 1,
"DATA_SELECTION_MAP", 2,
"DATA_TYPE", 3,
"DATA_SAVE_IN_STATE", 4);
});
