Clazz.declarePackage ("J.dssx");
Clazz.load (["J.api.JmolDSSRParser"], "J.dssx.DSSRParser", ["java.lang.Boolean", "$.Character", "$.Float", "java.util.Hashtable", "JU.BS", "$.Lst", "$.P3", "$.PT", "$.Rdr", "$.SB", "JM.HBond", "JM.BasePair", "$.NucleicPolymer", "JS.SV", "JU.Logger"], function () {
c$ = Clazz.decorateAsClass (function () {
this.reader = null;
this.line = null;
this.dssr = null;
this.htTemp = null;
this.message = null;
this.htPar = null;
this.basePairs = null;
this.next = null;
Clazz.instantialize (this, arguments);
}, J.dssx, "DSSRParser", null, J.api.JmolDSSRParser);
Clazz.prepareFields (c$, function () {
this.next =  Clazz.newIntArray (1, 0);
});
Clazz.makeConstructor (c$, 
function () {
});
Clazz.overrideMethod (c$, "process", 
function (info, reader, line0, htGroup1) {
info.put ("dssr", this.dssr =  new java.util.Hashtable ());
this.htTemp =  new java.util.Hashtable ();
this.htPar =  new java.util.Hashtable ();
this.htPar.put ("bp", ["bpShear", "bpStretch", "bpStagger", "bpPropeller", "bpBuckle", "bpOpening"]);
this.htPar.put ("bpChiLambda", ["bpChi1", "bpLambda1", "bpChi2", "bpLambda2"]);
this.htPar.put ("bpDistTor", ["bpDistC1C1", "bpDistNN", "bpDistC6C8", "bpTorCNNC"]);
this.htPar.put ("step", ["stShift", "stSlide", "stRise", "stTilt", "stRoll", "stTwist"]);
this.htPar.put ("hel", ["heXDisp", "heYDisp", "heRise", "heIncl", "heTip", "heTwist"]);
this.reader = reader;
this.message =  new JU.SB ();
this.line = (line0 == null ? "" : line0.trim ());
this.skipTo ("DSSR:", false);
this.addToMessages (null);
var haveHeader = false;
while (this.rd () != null) {
if (this.line.startsWith ("List of")) {
var n = JU.PT.parseInt (this.line.substring (8));
if (n < 0 || this.line.endsWith ("files")) continue;
this.addMessage (this.line);
this.line = JU.PT.rep (JU.PT.trim (this.line, "s"), " interaction", "");
var pt = "pair elix lice plet stem tack loop ulge tion ment otif pper turn bond file".indexOf (this.line.trim ().substring (this.line.length - 4));
switch (pt) {
case 0:
this.readPairs (n);
break;
case 5:
case 10:
this.getHelixOrStem (n, "helices", "helix", true);
break;
case 15:
this.readNTList (null, "multiplets", n);
break;
case 20:
this.getHelixOrStem (n, "stems", "stem", false);
break;
case 25:
this.readStacks (n);
break;
case 30:
this.readLoops (n);
break;
case 35:
this.readBulges (n);
break;
case 40:
this.readJunctions (n);
break;
case 45:
this.readNTList (null, "singleStranded", n);
break;
case 50:
this.readMotifs (n);
break;
case 55:
this.readNTList (null, "riboseZippers", n);
break;
case 60:
this.readTurns (n);
break;
case 65:
this.readHBonds (n);
break;
case 70:
break;
default:
this.addMessage ("DSSRParser ignored: " + this.line);
break;
}
} else if (!haveHeader && this.line.startsWith ("Date and time")) {
haveHeader = true;
this.addToMessages ("");
} else if (this.line.startsWith ("Secondary structures in dot-bracket")) {
this.readStructure ();
} else if (this.line.startsWith ("Mapping of")) {
this.mapGroups (htGroup1);
}}
this.dssr.put ("summary", this.message.toString ());
return this.message.toString ();
}, "java.util.Map,javajs.api.GenericLineReader,~S,java.util.Map");
Clazz.defineMethod (c$, "mapGroups", 
 function (map) {
this.rd ();
var n = 0;
var s = "";
if (JM.NucleicPolymer.htGroup1 == null) JM.NucleicPolymer.htGroup1 =  new java.util.Hashtable ();
while (this.rd () != null && this.line.length > 21) {
var g3 = this.line.substring (9, 12).trim ();
var g1 = this.line.substring (21, 22);
if ("ACGTU".indexOf (g1) < 0) {
JM.NucleicPolymer.htGroup1.put (g3, g1);
var key = " " + g3 + "(" + g1 + ")";
if (s.indexOf (key) < 0) {
n++;
s += key;
}if (map != null) map.put (g3, g1);
}}
if (n > 0) this.addMessage (n + " nonstandard base" + (n > 1 ? "s" : "") + ":" + s);
}, "java.util.Map");
Clazz.defineMethod (c$, "readStructure", 
 function () {
this.addMessage ("");
this.addMessage (this.line);
this.addMessage (this.rd ());
this.dssr.put ("seq", this.rd ());
this.addMessage (this.line);
this.dssr.put ("dbn", this.rd ());
this.addToMessages (this.line);
this.addMessage ("");
});
Clazz.defineMethod (c$, "readHBonds", 
 function (n) {
var list = this.newList ("hBonds");
for (var i = 0; i < n; i++) {
var data =  new java.util.Hashtable ();
var tokens = JU.PT.getTokens (this.rd ());
data.put ("atno1", Integer.$valueOf (JU.PT.parseInt (tokens[0])));
data.put ("atno2", Integer.$valueOf (JU.PT.parseInt (tokens[1])));
data.put ("id", tokens[2]);
data.put ("hbType", tokens[3]);
data.put ("distAng", Float.$valueOf (tokens[4]));
var pt = (tokens.length > 8 ? 6 : 5);
data.put ("energy", Float.$valueOf (pt == 6 ? tokens[5] : "0"));
data.put ("label", tokens[pt++]);
data.put ("atom1", this.fix (tokens[pt++], true));
data.put ("atom2", this.fix (tokens[pt++], true));
if (pt < tokens.length) data.put ("primary", Boolean.$valueOf (tokens[pt++].equals ("primary")));
list.addLast (data);
}
}, "~N");
Clazz.defineMethod (c$, "addToMessages", 
 function (s) {
if (s != null) this.addMessage (s);
while (this.line != null && this.line.length > 0 && this.line.indexOf ("****") < 0) {
this.addMessage (s == null ? this.line.trim () : this.line);
this.rd ();
}
}, "~S");
Clazz.defineMethod (c$, "addMessage", 
 function (s) {
this.message.append (s).append ("\n");
}, "~S");
Clazz.defineMethod (c$, "newList", 
 function (name) {
var list =  new JU.Lst ();
if (name != null) this.dssr.put (name, list);
return list;
}, "~S");
Clazz.defineMethod (c$, "readStacks", 
 function (n) {
var list = this.newList ("coaxialStacks");
for (var i = 0; i < n; i++) {
var data =  new java.util.Hashtable ();
var tokens = JU.PT.getTokens (this.rd ());
data.put ("helix", tokens[1]);
data.put ("stemCount", Integer.$valueOf (tokens[3]));
data.put ("stems", tokens[5]);
data.put ("basePairs", this.getLinkNTList (tokens[5], "stem", null));
list.addLast (data);
}
}, "~N");
Clazz.defineMethod (c$, "readInfo", 
 function (key, n) {
var list =  new JU.Lst ();
if (key != null) this.dssr.put (key, list);
for (var i = 0; i < n; i++) list.addLast (this.rd ());

return list;
}, "~S,~N");
Clazz.defineMethod (c$, "readLoops", 
 function (n) {
if (this.line.indexOf ("internal") >= 0) {
this.readSets ("internalLoops", n, 2, 4);
} else if (this.line.indexOf ("hairpin") >= 0) {
this.readSets ("hairpinLoops", n, 1, 3);
} else if (this.line.indexOf ("kissing") >= 0) {
this.readSets ("kissingLoops", n, -1, -1);
}}, "~N");
Clazz.defineMethod (c$, "readJunctions", 
 function (n) {
this.readSets ("junctions", n, 0, 3);
}, "~N");
Clazz.defineMethod (c$, "readBulges", 
 function (n) {
this.readSets ("bulges", n, 2, 2);
}, "~N");
Clazz.defineMethod (c$, "readSets", 
 function (key, n, nway, ptnts) {
var sets = this.newList (key);
var isJunction = (nway == 0);
var isKissingLoop = (ptnts == -1);
for (var i = 0; i < n; i++) {
var set =  new java.util.Hashtable ();
var tokens = JU.PT.getTokens (this.rd ());
set.put ("id", tokens[0]);
this.htTemp.put (key + tokens[0], set);
var lst =  new JU.Lst ();
set.put ("desc", this.line);
if (isKissingLoop) {
this.getNTs (this.getLinkNTList (tokens[2], "stem", null), lst, true, false);
this.getNTs (this.getLinkNTList (tokens[6], "hairpinLoops", null), lst, false, false);
this.getNTs (this.getLinkNTList (tokens[8], "hairpinLoops", null), lst, false, false);
set.put ("nts", lst);
lst =  new JU.Lst ();
this.getNTs (this.getLinkNTList (tokens[2], "stem", null), lst, true, true);
this.getNTs (this.getLinkNTList (tokens[6], "hairpinLoops", null), lst, false, true);
this.getNTs (this.getLinkNTList (tokens[8], "hairpinLoops", null), lst, false, true);
set.put ("resnos", lst);
} else {
set.put ("dssrType", tokens[1]);
if (isJunction) nway = JU.PT.parseInt (tokens[1].substring (0, tokens[1].indexOf ("-")));
set.put ("nway", Integer.$valueOf (nway));
set.put ("n", Integer.$valueOf (JU.PT.trim (tokens[ptnts], ";").substring (4)));
set.put ("linkedBy", this.getLinkNTList (tokens[ptnts + 4], "stem", lst));
set.put ("basePairs", this.readNTList (key + "#" + (i + 1), null, nway + 1));
}sets.addLast (set);
}
}, "~S,~N,~N,~N");
Clazz.defineMethod (c$, "getNTs", 
 function (linkNTList, lst, isStem, isResno) {
var o = linkNTList.get (0);
var n = o.size ();
var key = (!isResno ? "nt" : isStem ? "res" : "resno");
var nts = (isStem ?  new JU.Lst () : null);
for (var i = 0; i < n; i++) {
var m = o.get (i);
if (isStem) {
nts.addLast (m.get (key + "1"));
nts.addLast (m.get (key + "2"));
} else {
lst.addLast (m.get (key + "s"));
}}
if (isStem) {
lst.addLast (nts);
}}, "JU.Lst,JU.Lst,~B,~B");
Clazz.defineMethod (c$, "getLinkNTList", 
 function (linkStr, type, list) {
if (list == null) list =  new JU.Lst ();
var tokens = JU.PT.getTokens (JU.PT.replaceAllCharacters (linkStr, "[,]", " "));
for (var i = 0; i < tokens.length; i++) list.addLast (this.htTemp.get ((tokens[i].startsWith ("-") ? "" : type) + tokens[i]));

return list;
}, "~S,~S,JU.Lst");
Clazz.defineMethod (c$, "readMotifs", 
 function (n) {
var motifs = this.newList ("aMinorMotifs");
for (var i = 0; i < n; i++) {
var motif =  new java.util.Hashtable ();
var tokens = JU.PT.getTokens (this.rd ());
motif.put ("motiftype", this.after (tokens[1], "=") + " " + tokens[2]);
motif.put ("info", this.line);
motif.put ("data", this.readInfo (null, 2));
motifs.addLast (motif);
}
}, "~N");
Clazz.defineMethod (c$, "readTurns", 
 function (n) {
var turns = this.newList ("kinkTurns");
for (var i = 0; i < n; i++) {
var turn =  new java.util.Hashtable ();
var tokens = JU.PT.getTokens (this.rd ());
turn.put ("turnType", tokens[1]);
turn.put ("info", this.line);
turn.put ("details", this.rd ());
turn.put ("basePairs", this.readNTList (null, null, 2));
turns.addLast (turn);
}
}, "~N");
Clazz.defineMethod (c$, "readPairs", 
 function (n) {
var pairs;
if (this.line.indexOf ("lone ") >= 0) {
this.rd ();
this.skipHeader ();
pairs = this.newList ("lonePairs");
for (var i = 0; i < n; i++) {
var tokens = JU.PT.getTokens (this.line);
var data = this.htTemp.get (tokens[1] + tokens[2]);
this.htTemp.put ("#" + this.line.substring (0, 5).trim (), data);
data.put ("lonePair", Boolean.TRUE);
pairs.addLast (data);
this.rd ();
}
return;
}this.basePairs = this.newList ("basePairs");
this.skipHeader ();
for (var i = 0; i < n; i++) this.getBPData (0, null, true);

}, "~N");
Clazz.defineMethod (c$, "getBPData", 
 function (i0, type, readParams) {
var tokens = JU.PT.getTokens (this.line);
var nt12 = tokens[1] + tokens[2];
var data;
data = this.htTemp.get (nt12);
var i = i0;
if (data == null) {
data =  new java.util.Hashtable ();
i = JU.PT.parseInt (tokens[0]);
if (type != null) i = -((this.htTemp.get (tokens[2] + tokens[1])).get ("id")).intValue ();
data.put ("id", Integer.$valueOf (i));
var nt1 = this.fix (tokens[1], true);
var nt2 = this.fix (tokens[2], true);
data.put ("key", nt1 + " " + nt2);
data.put ("nt1", nt1);
data.put ("nt2", nt2);
data.put ("nt2", this.fix (tokens[2], true));
data.put ("res1", this.fix (tokens[1], false));
data.put ("res2", this.fix (tokens[2], false));
var bp = tokens[3];
data.put ("bp", bp);
data.put ("g1", bp.substring (0, 1));
data.put ("g2", bp.substring (2, 3));
var pt = (tokens.length == 8 ? 5 : 4);
data.put ("name", pt == 5 ? tokens[4] : "?");
var pt1 = tokens[pt].indexOf ("-");
data.put ("Saenger", Integer.$valueOf (pt1 > 0 ? tokens[pt].substring (0, pt1) : "0"));
data.put ("LW", tokens[++pt]);
data.put ("DSSR", tokens[++pt]);
this.htTemp.put (nt12, data);
this.basePairs.addLast (data);
}if (type != null) data.put (type + "Id", Integer.$valueOf (i0));
if (readParams) this.readMore (data, type == null, i < 0);
 else this.skipHeader ();
return data;
}, "~N,~S,~B");
Clazz.defineMethod (c$, "readMore", 
 function (data, isBP, isRev) {
var info = "";
while (this.isHeader (this.rd ())) {
var pt = this.line.indexOf ("[");
this.line = JU.PT.rep (this.line, "_pars", "-pars");
if (isBP) {
if (this.line.indexOf ("bp-pars:") >= 0) {
this.addArray (data, "bp", JU.PT.parseFloatArray (this.line.substring (pt + 1)));
} else if (this.line.indexOf ("lambda") >= 0) {
this.extractFloats (data, this.htPar.get ("bpChiLambda"));
} else if (this.line.indexOf ("tor(") >= 0) {
this.extractFloats (data, this.htPar.get ("bpDistTor"));
}info += this.line + "\n";
} else {
if (isRev && this.line.indexOf ("bp1-pars:") >= 0) {
this.addArray (data, "bp", JU.PT.parseFloatArray (this.line.substring (pt + 1)));
} else if (this.line.indexOf ("heli-pars:") >= 0) {
this.addArray (data, "hel", JU.PT.parseFloatArray (this.line.substring (pt + 1)));
} else if (this.line.indexOf ("step-pars:") >= 0) {
this.addArray (data, "step", JU.PT.parseFloatArray (this.line.substring (pt + 1)));
} else if ((pt = this.line.indexOf ("h-rise=")) >= 0) {
this.addFloat (data, "heRiseC1", pt + 7);
this.addFloat (data, "heTwistC1", this.line.indexOf ("h-twist=") + 8);
} else if ((pt = this.line.indexOf ("rise=")) >= 0) {
this.addFloat (data, "stRiseC1", pt + 5);
this.addFloat (data, "stTwistC1", this.line.indexOf ("twist=") + 6);
}}}
if (isBP) data.put ("info", info);
}, "java.util.Map,~B,~B");
Clazz.defineMethod (c$, "extractFloats", 
 function (data, names) {
this.line = this.line.$replace ('[', '=').$replace ('(', ' ').$replace (']', ' ');
this.next[0] = -1;
var n = names.length;
for (var i = 0, pt = 0; i < n; i++) {
if ((this.next[0] = pt = this.line.indexOf ("=", pt) + 1) == 0) break;
data.put (names[i], Float.$valueOf (JU.PT.parseFloatNext (this.line, this.next)));
}
}, "java.util.Map,~A");
Clazz.defineMethod (c$, "addArray", 
 function (data, key, f) {
var keys = this.htPar.get (key);
var n = Math.min (f.length, keys == null ? f.length : keys.length);
for (var i = 0; i < n; i++) data.put (keys == null ? key + (i + 1) : keys[i], Float.$valueOf (f[i]));

}, "java.util.Map,~S,~A");
Clazz.defineMethod (c$, "addFloat", 
 function (data, key, pt) {
data.put (key, Float.$valueOf (JU.PT.parseFloat (this.line.substring (pt, Math.min (this.line.length, pt + 10)))));
}, "java.util.Map,~S,~N");
Clazz.defineMethod (c$, "readNTList", 
 function (ntsKey, type, n) {
var isHairpin = (n == 2);
var list = this.newList (type);
if (ntsKey != null) this.htTemp.put (ntsKey, list);
if (isHairpin) this.rd ();
for (var i = (isHairpin ? 1 : 0); i < n; i++) list.addLast (this.getNTList ());

return list;
}, "~S,~S,~N");
Clazz.defineMethod (c$, "getHelixOrStem", 
 function (n, key, type, isHelix) {
var list = this.newList (key);
for (var i = 0; i < n; i++) {
this.skipTo ("  " + type + "#", true);
var bps = JU.PT.parseInt (this.after (this.line, "="));
var data =  new java.util.Hashtable ();
var header = this.getHeader ();
data.put ("info", header);
data.put ("bpCount", Integer.$valueOf (bps));
if (isHelix) {
var lines = JU.PT.split (header, "\n");
if (lines.length == 8) {
data.put ("helicalAxisData", this.after (lines[5], "s"));
data.put ("p1", this.getPoint (lines[6]));
data.put ("p2", this.getPoint (lines[7]));
}}list.addLast (data);
var pairs = this.newList (null);
data.put ("basePairs", pairs);
this.htTemp.put (type + "#" + (i + 1), pairs);
for (var j = 0; j < bps; j++) pairs.addLast (this.getBPData (i + 1, type, isHelix));

}
}, "~N,~S,~S,~B");
Clazz.defineMethod (c$, "getPoint", 
 function (data) {
var a = JU.PT.parseFloatArray (this.after (data, ":"));
return JU.P3.new3 (a[0], a[1], a[2]);
}, "~S");
Clazz.defineMethod (c$, "getNTList", 
 function () {
var data =  new java.util.Hashtable ();
var tokens = JU.PT.getTokens (this.rd ());
var pt = (tokens[0].startsWith ("nts") ? 0 : 1);
if (tokens.length > pt + 2) {
data.put ("nres", Integer.$valueOf (JU.PT.replaceAllCharacters (this.after (tokens[pt], "="), "*;", "")));
data.put ("seq", tokens[++pt]);
data.put ("nts", this.getNT (tokens[++pt], false));
data.put ("resnos", this.getNT (tokens[pt], true));
}return data;
});
Clazz.defineMethod (c$, "getNT", 
 function (s, isResno) {
var tokens = JU.PT.split (s, ",");
var list =  new JU.Lst ();
for (var i = 0; i < tokens.length; i++) list.addLast (this.fix (tokens[i], !isResno));

return list;
}, "~S,~B");
Clazz.defineMethod (c$, "getHeader", 
 function () {
var header =  new JU.SB ();
header.append (this.line).append ("\n");
while (this.isHeader (this.rd ())) header.append (this.line).append ("\n");

return header.toString ();
});
Clazz.defineMethod (c$, "skipHeader", 
 function () {
while (this.isHeader (this.rd ())) {
}
});
Clazz.defineMethod (c$, "isHeader", 
 function (line) {
return line.length < 6 || line.charAt (3) == ' ' || line.charAt (5) == ' ';
}, "~S");
Clazz.defineMethod (c$, "skipTo", 
 function (key, startsWith) {
while (!(startsWith ? this.line.startsWith (key) : this.line.contains (key))) {
this.rd ();
}
}, "~S,~B");
Clazz.defineMethod (c$, "fix", 
 function (nt, withName) {
var pt1;
if (nt.startsWith ("[")) {
if ((pt1 = nt.indexOf ("/")) >= 0) nt = nt.substring (0, pt1);
if (withName) return nt;
if ((pt1 = nt.indexOf (".")) >= 0) nt = nt.substring (0, pt1);
return (nt.substring (nt.indexOf ("]") + 1));
}pt1 = nt.indexOf (".");
var chain = nt.substring (0, pt1);
var pt = nt.length;
var ch;
while (Character.isDigit (ch = nt.charAt (--pt))) {
}
var ptn = chain.indexOf ("@");
if (ptn >= 0) chain = chain.substring (ptn + 1) + (withName ? "." + chain.substring (0, ptn) : "");
var pt2 = (ch == '/' ? pt : pt + 1);
return (withName ? "[" + nt.substring (pt1 + 1, pt2) + "]" : "") + nt.substring (pt + 1) + ":" + chain;
}, "~S,~B");
Clazz.defineMethod (c$, "after", 
 function (s, key) {
return s.substring (s.indexOf (key) + 1);
}, "~S,~S");
Clazz.defineMethod (c$, "rd", 
 function () {
this.line = this.reader.readNextLine ();
if (JU.Logger.debugging) JU.Logger.info (this.line);
return this.line;
});
Clazz.overrideMethod (c$, "getAtomBits", 
function (vwr, key, dssr, dssrCache) {
var s = key.toLowerCase ();
if (s.indexOf ("pairs") < 0 && s.indexOf ("kissingloops") < 0 && s.indexOf ("linkedby") < 0 && s.indexOf ("multiplets") < 0 && s.indexOf ("singlestrand") < 0) key += ".basePairs";
if (s.indexOf (".nt") < 0 && s.indexOf (".res") < 0 && s.indexOf ("[selecet res") < 0 && s.indexOf ("[selecet nt") < 0) key += ".res*";
var bs = dssrCache.get (key);
var htChains =  new java.util.Hashtable ();
if (bs == null) {
dssrCache.put (key, bs =  new JU.BS ());
var data = vwr.extractProperty (dssr, key, -1);
if (Clazz.instanceOf (data, JU.Lst)) this.getBsAtoms (vwr, null, data, bs, htChains);
}return bs;
}, "JV.Viewer,~S,~O,java.util.Map");
Clazz.defineMethod (c$, "getBsAtoms", 
 function (vwr, res, lst, bs, htChains) {
var tokens;
if (lst == null) {
tokens = JU.PT.getTokens (JU.PT.replaceAllCharacters (res.toString (), "=[,]", " "));
} else if (lst.size () == 0) {
return;
} else {
tokens =  new Array (lst.size ());
for (var i = lst.size (); --i >= 0; ) {
var o = lst.get (i);
if (Clazz.instanceOf (o, JS.SV)) o = (o).value;
if (Clazz.instanceOf (o, JU.Lst)) {
this.getBsAtoms (vwr, null, o, bs, htChains);
} else {
var s = (Clazz.instanceOf (o, JS.SV) ? (o).asString () : o.toString ());
tokens[i] = (s.startsWith ("[") ? s.substring (s.indexOf ("]") + 1) : s);
}}
}for (var j = tokens.length; --j >= 0; ) {
var t = tokens[j];
if (t == null) continue;
var pt = t.indexOf (":");
if (pt < 0 || pt + 1 == t.length) continue;
var chain = t.substring (pt + 1);
var bsChain = htChains.get (chain);
try {
if (bsChain == null) htChains.put (chain, bsChain = vwr.ms.getAtoms (1048609, Integer.$valueOf (vwr.getChainID (chain))));
var bsRes = vwr.ms.getAtoms (1095761939, Integer.$valueOf (t.substring (0, pt)));
bsRes.and (bsChain);
bs.or (bsRes);
} catch (e) {
if (Clazz.exceptionOf (e, Exception)) {
} else {
throw e;
}
}
}
}, "JV.Viewer,~S,JU.Lst,JU.BS,java.util.Map");
Clazz.overrideMethod (c$, "setAllDSSRParametersForModel", 
function (vwr, modelIndex) {
var dssr = vwr.ms.getInfo (modelIndex, "dssr");
var lst = (dssr == null ? null : dssr.get ("basePairs"));
var lst1 = (dssr == null ? null : dssr.get ("singleStranded"));
if (lst == null && lst1 == null) {
var m = vwr.ms.am[modelIndex];
var n = m.getBioPolymerCount ();
for (var i = n; --i >= 0; ) {
var bp = m.getBioPolymer (i);
if (bp.isNucleic ()) (bp).isDssrSet = true;
}
return;
}var htChains =  new java.util.Hashtable ();
var bs =  new JU.BS ();
if (lst != null) {
for (var i = lst.size (); --i >= 0; ) {
var bpInfo = lst.get (i);
JM.BasePair.add (bpInfo, this.setPhos (vwr, 1, bpInfo, bs, htChains), this.setPhos (vwr, 2, bpInfo, bs, htChains));
}
}if (lst1 != null) for (var i = lst1.size (); --i >= 0; ) {
var bp = lst1.get (i);
var resnos = bp.get ("resnos");
for (var j = resnos.size (); --j >= 0; ) this.setRes (vwr, resnos.get (j), bs, htChains);

}
}, "JV.Viewer,~N");
Clazz.defineMethod (c$, "setPhos", 
 function (vwr, n, bp, bs, htChains) {
return this.setRes (vwr, bp.get ("res" + n), bs, htChains);
}, "JV.Viewer,~N,java.util.Map,JU.BS,java.util.Map");
Clazz.defineMethod (c$, "setRes", 
 function (vwr, res, bs, htChains) {
bs.clearAll ();
this.getBsAtoms (vwr, res, null, bs, htChains);
var group = vwr.ms.at[bs.nextSetBit (0)].getGroup ();
(group.bioPolymer).isDssrSet = true;
return group;
}, "JV.Viewer,~S,JU.BS,java.util.Map");
Clazz.overrideMethod (c$, "getHBonds", 
function (ms, modelIndex, vHBonds, doReport) {
var info = ms.getInfo (modelIndex, "dssr");
if (info != null) info = (info).get ("hBonds");
if (info == null) return "no DSSR hydrogen-bond data";
var list = info;
var a0 = ms.am[modelIndex].firstAtomIndex - 1;
try {
for (var i = list.size (); --i >= 0; ) {
var hbond = list.get (i);
var a1 = (hbond.get ("atno1")).intValue () + a0;
var a2 = (hbond.get ("atno2")).intValue () + a0;
var energy = (hbond.containsKey ("energy") ? (hbond.get ("energy")).floatValue () : 0);
vHBonds.addLast ( new JM.HBond (ms.at[a1], ms.at[a2], 2048, 1, 0, energy));
}
} catch (e) {
if (Clazz.exceptionOf (e, Exception)) {
JU.Logger.error ("Exception " + e + " in DSSRParser.getHBonds");
} else {
throw e;
}
}
return "DSSR reports " + list.size () + " hydrogen bonds";
}, "JM.ModelSet,~N,JU.Lst,~B");
Clazz.overrideMethod (c$, "calculateStructure", 
function (vwr, bsAtoms) {
var bs = vwr.ms.getModelBS (bsAtoms == null ? vwr.bsA () : bsAtoms, true);
var s = "";
for (var i = bs.nextSetBit (0); i >= 0; i = bs.nextSetBit (i + 1)) s += this.getDSSRForModel (vwr, i);

return s;
}, "JV.Viewer,JU.BS");
Clazz.defineMethod (c$, "getDSSRForModel", 
 function (vwr, modelIndex) {
var info = null;
var out = null;
while (true) {
if (!vwr.ms.am[modelIndex].isBioModel) break;
info = vwr.ms.getModelAuxiliaryInfo (modelIndex);
if (info.containsKey ("dssr")) break;
var bs = vwr.getModelUndeletedAtomsBitSet (modelIndex);
bs.and (vwr.ms.getAtoms (3145742, null));
if (bs.nextClearBit (0) < 0) {
info = null;
break;
}try {
var name = vwr.setLoadFormat ("=dssrModel/", '=', false);
name = JU.PT.rep (name, "%20", " ");
JU.Logger.info ("fetching " + name + "[pdb data]");
var data = vwr.getPdbAtomData (bs, null);
data = vwr.getFileAsString (name + data, false);
this.process (info,  new JU.Rdr (JU.Rdr.getBR (data)), null, null);
} catch (e) {
if (Clazz.exceptionOf (e, Exception)) {
info = null;
out = "" + e;
} else {
throw e;
}
}
break;
}
return (info != null ? (info.get ("dssr")).get ("summary") : out == null ? "model has no nucleotides" : out);
}, "JV.Viewer,~N");
});
