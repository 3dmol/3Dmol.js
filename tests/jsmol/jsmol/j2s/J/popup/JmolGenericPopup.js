Clazz.declarePackage ("J.popup");
Clazz.load (["J.popup.GenericSwingPopup", "java.util.Properties", "JU.Lst"], "J.popup.JmolGenericPopup", ["java.lang.Boolean", "java.util.Hashtable", "JU.PT", "J.i18n.GT", "J.popup.MainPopupResourceBundle", "JU.Elements", "JV.JC"], function () {
c$ = Clazz.decorateAsClass (function () {
this.vwr = null;
this.updateMode = 0;
this.menuText = null;
this.frankPopup = null;
this.nFrankList = 0;
this.itemMax = 25;
this.titleWidthMax = 20;
this.nullModelSetName = null;
this.modelSetName = null;
this.modelSetFileName = null;
this.modelSetRoot = null;
this.currentFrankId = null;
this.configurationSelected = "";
this.altlocs = null;
this.frankList = null;
this.modelSetInfo = null;
this.modelInfo = null;
this.NotPDB = null;
this.PDBOnly = null;
this.FileUnitOnly = null;
this.FileMolOnly = null;
this.UnitcellOnly = null;
this.SingleModelOnly = null;
this.FramesOnly = null;
this.VibrationOnly = null;
this.Special = null;
this.SymmetryOnly = null;
this.ChargesOnly = null;
this.TemperatureOnly = null;
this.fileHasUnitCell = false;
this.haveBFactors = false;
this.haveCharges = false;
this.isLastFrame = false;
this.isMultiConfiguration = false;
this.isMultiFrame = false;
this.isPDB = false;
this.isSymmetry = false;
this.isUnitCell = false;
this.isVibration = false;
this.isZapped = false;
this.modelIndex = 0;
this.modelCount = 0;
this.ac = 0;
this.group3List = null;
this.group3Counts = null;
this.cnmrPeaks = null;
this.hnmrPeaks = null;
this.noZapped = null;
Clazz.instantialize (this, arguments);
}, J.popup, "JmolGenericPopup", J.popup.GenericSwingPopup);
Clazz.prepareFields (c$, function () {
this.menuText =  new java.util.Properties ();
this.frankList =  new Array (10);
this.NotPDB =  new JU.Lst ();
this.PDBOnly =  new JU.Lst ();
this.FileUnitOnly =  new JU.Lst ();
this.FileMolOnly =  new JU.Lst ();
this.UnitcellOnly =  new JU.Lst ();
this.SingleModelOnly =  new JU.Lst ();
this.FramesOnly =  new JU.Lst ();
this.VibrationOnly =  new JU.Lst ();
this.Special =  new JU.Lst ();
this.SymmetryOnly =  new JU.Lst ();
this.ChargesOnly =  new JU.Lst ();
this.TemperatureOnly =  new JU.Lst ();
this.noZapped = ["surfaceMenu", "measureMenu", "pickingMenu", "computationMenu", "saveMenu", "exportMenu", "SIGNEDJAVAcaptureMenuSPECIAL"];
});
Clazz.defineMethod (c$, "initialize", 
function (vwr, bundle, title) {
this.vwr = vwr;
this.initSwing (title, bundle, vwr.getApplet (), vwr.isJS, vwr.getBooleanProperty ("_signedApplet"), vwr.isWebGL);
}, "JV.Viewer,J.popup.PopupResource,~S");
Clazz.overrideMethod (c$, "jpiDispose", 
function () {
this.helper.menuClearListeners (this.popupMenu);
this.helper.menuClearListeners (this.frankPopup);
this.popupMenu = this.frankPopup = this.thisPopup = null;
});
Clazz.overrideMethod (c$, "jpiGetMenuAsObject", 
function () {
return this.popupMenu;
});
Clazz.overrideMethod (c$, "jpiShow", 
function (x, y) {
if (!this.vwr.haveDisplay) return;
this.show (x, y, false);
if (x < 0) {
this.getViewerData ();
this.setFrankMenu (this.currentMenuItemId);
this.thisx = -x - 50;
if (this.nFrankList > 1) {
this.thisy = y - this.nFrankList * 20;
this.menuShowPopup (this.frankPopup, this.thisx, this.thisy);
return;
}}this.appRestorePopupMenu ();
this.menuShowPopup (this.popupMenu, this.thisx, this.thisy);
}, "~N,~N");
Clazz.overrideMethod (c$, "jpiUpdateComputedMenus", 
function () {
if (this.updateMode == -1) return;
this.isTainted = true;
this.updateMode = 0;
this.getViewerData ();
this.updateSelectMenu ();
this.updateFileMenu ();
this.updateElementsComputedMenu (this.vwr.getElementsPresentBitSet (this.modelIndex));
this.updateHeteroComputedMenu (this.vwr.getHeteroList (this.modelIndex));
this.updateSurfMoComputedMenu (this.modelInfo.get ("moData"));
this.updateFileTypeDependentMenus ();
this.updatePDBComputedMenus ();
this.updateMode = 1;
this.updateConfigurationComputedMenu ();
this.updateSYMMETRYComputedMenus ();
this.updateFRAMESbyModelComputedMenu ();
this.updateModelSetComputedMenu ();
this.updateLanguageSubmenu ();
this.updateAboutSubmenu ();
});
Clazz.overrideMethod (c$, "appCheckItem", 
function (item, newMenu) {
if (item.indexOf ("!PDB") >= 0) {
this.NotPDB.addLast (newMenu);
} else if (item.indexOf ("PDB") >= 0) {
this.PDBOnly.addLast (newMenu);
}if (item.indexOf ("CHARGE") >= 0) {
this.ChargesOnly.addLast (newMenu);
} else if (item.indexOf ("BFACTORS") >= 0) {
this.TemperatureOnly.addLast (newMenu);
} else if (item.indexOf ("UNITCELL") >= 0) {
this.UnitcellOnly.addLast (newMenu);
} else if (item.indexOf ("FILEUNIT") >= 0) {
this.FileUnitOnly.addLast (newMenu);
} else if (item.indexOf ("FILEMOL") >= 0) {
this.FileMolOnly.addLast (newMenu);
}if (item.indexOf ("!FRAMES") >= 0) {
this.SingleModelOnly.addLast (newMenu);
} else if (item.indexOf ("FRAMES") >= 0) {
this.FramesOnly.addLast (newMenu);
}if (item.indexOf ("VIBRATION") >= 0) {
this.VibrationOnly.addLast (newMenu);
} else if (item.indexOf ("SYMMETRY") >= 0) {
this.SymmetryOnly.addLast (newMenu);
}if (item.indexOf ("SPECIAL") >= 0) this.Special.addLast (newMenu);
}, "~S,javajs.api.SC");
Clazz.overrideMethod (c$, "appFixLabel", 
function (label) {
return label;
}, "~S");
Clazz.overrideMethod (c$, "appFixScript", 
function (id, script) {
var pt;
if (script === "" || id.endsWith ("Checkbox")) return script;
if (script.indexOf ("SELECT") == 0) {
return "select thisModel and (" + script.substring (6) + ")";
}if ((pt = id.lastIndexOf ("[")) >= 0) {
id = id.substring (pt + 1);
if ((pt = id.indexOf ("]")) >= 0) id = id.substring (0, pt);
id = id.$replace ('_', ' ');
if (script.indexOf ("[]") < 0) script = "[] " + script;
script = script.$replace ('_', ' ');
return JU.PT.rep (script, "[]", id);
} else if (script.indexOf ("?FILEROOT?") >= 0) {
script = JU.PT.rep (script, "FILEROOT?", this.modelSetRoot);
} else if (script.indexOf ("?FILE?") >= 0) {
script = JU.PT.rep (script, "FILE?", this.modelSetFileName);
} else if (script.indexOf ("?PdbId?") >= 0) {
script = JU.PT.rep (script, "PdbId?", "=xxxx");
}return script;
}, "~S,~S");
Clazz.overrideMethod (c$, "appGetBooleanProperty", 
function (name) {
return this.vwr.getBooleanProperty (name);
}, "~S");
Clazz.overrideMethod (c$, "appGetMenuAsString", 
function (title) {
return ( new J.popup.MainPopupResourceBundle (this.strMenuStructure, null)).getMenuAsText (title);
}, "~S");
Clazz.overrideMethod (c$, "appIsSpecialCheckBox", 
function (item, basename, what, TF) {
if (this.appGetBooleanProperty (basename) == TF) return true;
if (!basename.endsWith ("P!")) return false;
if (basename.indexOf ("??") >= 0) {
what = this.menuSetCheckBoxOption (item, basename, what);
} else {
if (!TF) return true;
what = "set picking " + basename.substring (0, basename.length - 2);
}this.appRunScript (what);
return true;
}, "javajs.api.SC,~S,~S,~B");
Clazz.overrideMethod (c$, "appRestorePopupMenu", 
function () {
this.thisPopup = this.popupMenu;
if (this.vwr.isJS || this.nFrankList < 2) return;
for (var i = this.nFrankList; --i > 0; ) {
var f = this.frankList[i];
this.helper.menuInsertSubMenu (f[0], f[1], (f[2]).intValue ());
}
this.nFrankList = 1;
});
Clazz.overrideMethod (c$, "appRunScript", 
function (script) {
this.vwr.evalStringQuiet (script);
}, "~S");
Clazz.overrideMethod (c$, "appUpdateSpecialCheckBoxValue", 
function (item, what, TF) {
if (what.indexOf ("#CONFIG") >= 0) {
this.configurationSelected = what;
this.updateConfigurationComputedMenu ();
this.updateModelSetComputedMenu ();
}}, "javajs.api.SC,~S,~B");
Clazz.defineMethod (c$, "setFrankMenu", 
 function (id) {
if (this.currentFrankId != null && this.currentFrankId === id && this.nFrankList > 0) return;
if (this.frankPopup == null) this.frankPopup = this.helper.menuCreatePopup ("Frank", this.vwr.getApplet ());
this.thisPopup = this.frankPopup;
this.menuRemoveAll (this.frankPopup, 0);
this.menuCreateItem (this.frankPopup, this.getMenuText ("mainMenuText"), "MAIN", "");
this.currentFrankId = id;
this.nFrankList = 0;
this.frankList[this.nFrankList++] = [null, null, null];
if (id != null) for (var i = id.indexOf (".", 2) + 1; ; ) {
var iNew = id.indexOf (".", i);
if (iNew < 0) break;
var menu = this.htMenus.get (id.substring (i, iNew));
this.frankList[this.nFrankList++] = [menu.getParent (), menu, Integer.$valueOf (this.vwr.isJS ? 0 : this.menuGetListPosition (menu))];
this.menuAddSubMenu (this.frankPopup, menu);
i = iNew + 1;
}
this.thisPopup = this.popupMenu;
}, "~S");
c$.checkBoolean = Clazz.defineMethod (c$, "checkBoolean", 
 function (info, key) {
return (info != null && info.get (key) === Boolean.TRUE);
}, "java.util.Map,~S");
Clazz.defineMethod (c$, "getViewerData", 
 function () {
this.modelSetName = this.vwr.getModelSetName ();
this.modelSetFileName = this.vwr.getModelSetFileName ();
var i = this.modelSetFileName.lastIndexOf (".");
this.isZapped = ("zapped".equals (this.modelSetName));
if (this.isZapped || "string".equals (this.modelSetFileName) || "files".equals (this.modelSetFileName) || "string[]".equals (this.modelSetFileName)) this.modelSetFileName = "";
this.modelSetRoot = this.modelSetFileName.substring (0, i < 0 ? this.modelSetFileName.length : i);
if (this.modelSetRoot.length == 0) this.modelSetRoot = "Jmol";
this.modelIndex = this.vwr.am.cmi;
this.modelCount = this.vwr.getModelCount ();
this.ac = this.vwr.getAtomCountInModel (this.modelIndex);
this.modelSetInfo = this.vwr.getModelSetAuxiliaryInfo ();
this.modelInfo = this.vwr.getModelAuxiliaryInfo (this.modelIndex);
if (this.modelInfo == null) this.modelInfo =  new java.util.Hashtable ();
this.isPDB = J.popup.JmolGenericPopup.checkBoolean (this.modelSetInfo, "isPDB");
this.isMultiFrame = (this.modelCount > 1);
this.isSymmetry = J.popup.JmolGenericPopup.checkBoolean (this.modelInfo, "hasSymmetry");
this.isUnitCell = this.modelInfo.containsKey ("notionalUnitcell");
this.fileHasUnitCell = (this.isPDB && this.isUnitCell || J.popup.JmolGenericPopup.checkBoolean (this.modelInfo, "fileHasUnitCell"));
this.isLastFrame = (this.modelIndex == this.modelCount - 1);
this.altlocs = this.vwr.getAltLocListInModel (this.modelIndex);
this.isMultiConfiguration = (this.altlocs.length > 0);
this.isVibration = (this.vwr.modelHasVibrationVectors (this.modelIndex));
this.haveCharges = (this.vwr.havePartialCharges ());
this.haveBFactors = (this.vwr.getBooleanProperty ("haveBFactors"));
this.cnmrPeaks = this.modelInfo.get ("jdxAtomSelect_13CNMR");
this.hnmrPeaks = this.modelInfo.get ("jdxAtomSelect_1HNMR");
});
Clazz.overrideMethod (c$, "appCheckSpecialMenu", 
function (item, subMenu, word) {
if ("modelSetMenu".equals (item)) {
this.nullModelSetName = word;
this.menuEnable (subMenu, false);
}}, "~S,javajs.api.SC,~S");
Clazz.overrideMethod (c$, "appUpdateForShow", 
function () {
if (this.updateMode == -1) return;
this.isTainted = true;
this.getViewerData ();
this.updateMode = 2;
this.updateSelectMenu ();
this.updateSpectraMenu ();
this.updateFRAMESbyModelComputedMenu ();
this.updateSceneComputedMenu ();
this.updateModelSetComputedMenu ();
this.updateAboutSubmenu ();
for (var i = this.Special.size (); --i >= 0; ) this.updateSpecialMenuItem (this.Special.get (i));

});
Clazz.defineMethod (c$, "updateFileMenu", 
 function () {
var menu = this.htMenus.get ("fileMenu");
if (menu == null) return;
var text = this.getMenuText ("writeFileTextVARIABLE");
menu = this.htMenus.get ("writeFileTextVARIABLE");
var ignore = (this.modelSetFileName.equals ("zapped") || this.modelSetFileName.equals (""));
if (ignore) {
this.menuSetLabel (menu, "");
this.menuEnable (menu, false);
} else {
this.menuSetLabel (menu, J.i18n.GT.o (J.i18n.GT._ (text), this.modelSetFileName));
this.menuEnable (menu, true);
}});
Clazz.defineMethod (c$, "getMenuText", 
 function (key) {
var str = this.menuText.getProperty (key);
return (str == null ? key : str);
}, "~S");
Clazz.defineMethod (c$, "updateSelectMenu", 
 function () {
var menu = this.htMenus.get ("selectMenuText");
if (menu == null) return;
this.menuEnable (menu, this.ac != 0);
this.menuSetLabel (menu, this.gti ("selectMenuText", this.vwr.getSelectionCount ()));
});
Clazz.defineMethod (c$, "updateElementsComputedMenu", 
 function (elementsPresentBitSet) {
var menu = this.htMenus.get ("elementsComputedMenu");
if (menu == null) return;
this.menuRemoveAll (menu, 0);
this.menuEnable (menu, false);
if (elementsPresentBitSet == null) return;
for (var i = elementsPresentBitSet.nextSetBit (0); i >= 0; i = elementsPresentBitSet.nextSetBit (i + 1)) {
var elementName = JU.Elements.elementNameFromNumber (i);
var elementSymbol = JU.Elements.elementSymbolFromNumber (i);
var entryName = elementSymbol + " - " + elementName;
this.menuCreateItem (menu, entryName, "SELECT " + elementName, null);
}
for (var i = 4; i < JU.Elements.altElementMax; ++i) {
var n = JU.Elements.elementNumberMax + i;
if (elementsPresentBitSet.get (n)) {
n = JU.Elements.altElementNumberFromIndex (i);
var elementName = JU.Elements.elementNameFromNumber (n);
var elementSymbol = JU.Elements.elementSymbolFromNumber (n);
var entryName = elementSymbol + " - " + elementName;
this.menuCreateItem (menu, entryName, "SELECT " + elementName, null);
}}
this.menuEnable (menu, true);
}, "JU.BS");
Clazz.defineMethod (c$, "updateSpectraMenu", 
 function () {
var menuh = this.htMenus.get ("hnmrMenu");
var menuc = this.htMenus.get ("cnmrMenu");
if (menuh != null) this.menuRemoveAll (menuh, 0);
if (menuc != null) this.menuRemoveAll (menuc, 0);
var menu = this.htMenus.get ("spectraMenu");
if (menu == null) return;
this.menuRemoveAll (menu, 0);
var isOK =  new Boolean (this.setSpectraMenu (menuh, this.hnmrPeaks) | this.setSpectraMenu (menuc, this.cnmrPeaks)).valueOf ();
if (isOK) {
if (menuh != null) this.menuAddSubMenu (menu, menuh);
if (menuc != null) this.menuAddSubMenu (menu, menuc);
}this.menuEnable (menu, isOK);
});
Clazz.defineMethod (c$, "setSpectraMenu", 
 function (menu, peaks) {
if (menu == null) return false;
this.menuEnable (menu, false);
var n = (peaks == null ? 0 : peaks.size ());
if (n == 0) return false;
for (var i = 0; i < n; i++) {
var peak = peaks.get (i);
var title = JU.PT.getQuotedAttribute (peak, "title");
var atoms = JU.PT.getQuotedAttribute (peak, "atoms");
if (atoms != null) this.menuCreateItem (menu, title, "select visible & (@" + JU.PT.rep (atoms, ",", " or @") + ")", "Focus" + i);
}
this.menuEnable (menu, true);
return true;
}, "javajs.api.SC,JU.Lst");
Clazz.defineMethod (c$, "updateHeteroComputedMenu", 
 function (htHetero) {
var menu = this.htMenus.get ("PDBheteroComputedMenu");
if (menu == null) return;
this.menuRemoveAll (menu, 0);
this.menuEnable (menu, false);
if (htHetero == null) return;
var n = 0;
for (var hetero, $hetero = htHetero.entrySet ().iterator (); $hetero.hasNext () && ((hetero = $hetero.next ()) || true);) {
var heteroCode = hetero.getKey ();
var heteroName = hetero.getValue ();
if (heteroName.length > 20) heteroName = heteroName.substring (0, 20) + "...";
var entryName = heteroCode + " - " + heteroName;
this.menuCreateItem (menu, entryName, "SELECT [" + heteroCode + "]", null);
n++;
}
this.menuEnable (menu, (n > 0));
}, "java.util.Map");
Clazz.defineMethod (c$, "updateSurfMoComputedMenu", 
 function (moData) {
var menu = this.htMenus.get ("surfMoComputedMenuText");
if (menu == null) return;
this.menuRemoveAll (menu, 0);
var mos = (moData == null ? null : (moData.get ("mos")));
var nOrb = (mos == null ? 0 : mos.size ());
var text = this.getMenuText ("surfMoComputedMenuText");
if (nOrb == 0) {
this.menuSetLabel (menu, J.i18n.GT.o (J.i18n.GT._ (text), ""));
this.menuEnable (menu, false);
return;
}this.menuSetLabel (menu, J.i18n.GT.i (J.i18n.GT._ (text), nOrb));
this.menuEnable (menu, true);
var subMenu = menu;
var nmod = (nOrb % this.itemMax);
if (nmod == 0) nmod = this.itemMax;
var pt = (nOrb > this.itemMax ? 0 : -2147483648);
for (var i = nOrb; --i >= 0; ) {
if (pt >= 0 && (pt++ % nmod) == 0) {
if (pt == nmod + 1) nmod = this.itemMax;
var id = "mo" + pt + "Menu";
subMenu = this.menuNewSubMenu (Math.max (i + 2 - nmod, 1) + "..." + (i + 1), this.menuGetId (menu) + "." + id);
this.menuAddSubMenu (menu, subMenu);
this.htMenus.put (id, subMenu);
pt = 1;
}var mo = mos.get (i);
var entryName = "#" + (i + 1) + " " + (mo.containsKey ("type") ? mo.get ("type") + " " : "") + (mo.containsKey ("symmetry") ? mo.get ("symmetry") + " " : "") + (mo.containsKey ("occupancy") ? "(" + (mo.get ("occupancy")).intValue () + ") " : "") + (mo.containsKey ("energy") ? mo.get ("energy") : "");
var script = "mo " + (i + 1);
this.menuCreateItem (subMenu, entryName, script, null);
}
}, "java.util.Map");
Clazz.defineMethod (c$, "updateFileTypeDependentMenus", 
 function () {
for (var i = this.NotPDB.size (); --i >= 0; ) this.menuEnable (this.NotPDB.get (i), !this.isPDB);

for (var i = this.PDBOnly.size (); --i >= 0; ) this.menuEnable (this.PDBOnly.get (i), this.isPDB);

for (var i = this.UnitcellOnly.size (); --i >= 0; ) this.menuEnable (this.UnitcellOnly.get (i), this.isUnitCell);

for (var i = this.FileUnitOnly.size (); --i >= 0; ) this.menuEnable (this.FileUnitOnly.get (i), this.isUnitCell || this.fileHasUnitCell);

for (var i = this.FileMolOnly.size (); --i >= 0; ) this.menuEnable (this.FileMolOnly.get (i), this.isUnitCell || this.fileHasUnitCell);

for (var i = this.SingleModelOnly.size (); --i >= 0; ) this.menuEnable (this.SingleModelOnly.get (i), this.isLastFrame);

for (var i = this.FramesOnly.size (); --i >= 0; ) this.menuEnable (this.FramesOnly.get (i), this.isMultiFrame);

for (var i = this.VibrationOnly.size (); --i >= 0; ) this.menuEnable (this.VibrationOnly.get (i), this.isVibration);

for (var i = this.SymmetryOnly.size (); --i >= 0; ) this.menuEnable (this.SymmetryOnly.get (i), this.isSymmetry && this.isUnitCell);

for (var i = this.ChargesOnly.size (); --i >= 0; ) this.menuEnable (this.ChargesOnly.get (i), this.haveCharges);

for (var i = this.TemperatureOnly.size (); --i >= 0; ) this.menuEnable (this.TemperatureOnly.get (i), this.haveBFactors);

this.updateSignedAppletItems ();
});
Clazz.defineMethod (c$, "updateSceneComputedMenu", 
 function () {
var menu = this.htMenus.get ("sceneComputedMenu");
if (menu == null) return;
this.menuRemoveAll (menu, 0);
this.menuEnable (menu, false);
var scenes = this.vwr.getSceneList ();
if (scenes == null) return;
for (var i = 0; i < scenes.length; i++) this.menuCreateItem (menu, scenes[i], "restore scene " + JU.PT.esc (scenes[i]) + " 1.0", null);

this.menuEnable (menu, true);
});
Clazz.defineMethod (c$, "updatePDBComputedMenus", 
 function () {
var menu = this.htMenus.get ("PDBaaResiduesComputedMenu");
if (menu == null) return;
this.menuRemoveAll (menu, 0);
this.menuEnable (menu, false);
var menu1 = this.htMenus.get ("PDBnucleicResiduesComputedMenu");
if (menu1 == null) return;
this.menuRemoveAll (menu1, 0);
this.menuEnable (menu1, false);
var menu2 = this.htMenus.get ("PDBcarboResiduesComputedMenu");
if (menu2 == null) return;
this.menuRemoveAll (menu2, 0);
this.menuEnable (menu2, false);
if (this.modelSetInfo == null) return;
var n = (this.modelIndex < 0 ? 0 : this.modelIndex + 1);
var lists = (this.modelSetInfo.get ("group3Lists"));
this.group3List = (lists == null ? null : lists[n]);
this.group3Counts = (lists == null ? null : (this.modelSetInfo.get ("group3Counts"))[n]);
if (this.group3List == null) return;
var nItems = 0;
for (var i = 1; i < 24; ++i) nItems += this.updateGroup3List (menu, JV.JC.predefinedGroup3Names[i]);

nItems += this.augmentGroup3List (menu, "p>", true);
this.menuEnable (menu, (nItems > 0));
this.menuEnable (this.htMenus.get ("PDBproteinMenu"), (nItems > 0));
nItems = this.augmentGroup3List (menu1, "n>", false);
this.menuEnable (menu1, nItems > 0);
this.menuEnable (this.htMenus.get ("PDBnucleicMenu"), (nItems > 0));
nItems = this.augmentGroup3List (menu2, "c>", false);
this.menuEnable (menu2, nItems > 0);
this.menuEnable (this.htMenus.get ("PDBcarboMenu"), (nItems > 0));
});
Clazz.defineMethod (c$, "updateGroup3List", 
 function (menu, name) {
var nItems = 0;
var n = this.group3Counts[Clazz.doubleToInt (this.group3List.indexOf (name) / 6)];
var script = null;
if (n > 0) {
script = "SELECT " + name;
name += "  (" + n + ")";
nItems++;
}var item = this.menuCreateItem (menu, name, script, this.menuGetId (menu) + "." + name);
if (n == 0) this.menuEnable (item, false);
return nItems;
}, "javajs.api.SC,~S");
Clazz.defineMethod (c$, "augmentGroup3List", 
 function (menu, type, addSeparator) {
var pt = 138;
var nItems = 0;
while (true) {
pt = this.group3List.indexOf (type, pt);
if (pt < 0) break;
if (nItems++ == 0 && addSeparator) this.menuAddSeparator (menu);
var n = this.group3Counts[Clazz.doubleToInt (pt / 6)];
var heteroCode = this.group3List.substring (pt + 2, pt + 5);
var name = heteroCode + "  (" + n + ")";
this.menuCreateItem (menu, name, "SELECT [" + heteroCode + "]", this.menuGetId (menu) + "." + name);
pt++;
}
return nItems;
}, "javajs.api.SC,~S,~B");
Clazz.defineMethod (c$, "updateSYMMETRYComputedMenus", 
 function () {
this.updateSYMMETRYSelectComputedMenu ();
this.updateSYMMETRYShowComputedMenu ();
});
Clazz.defineMethod (c$, "updateSYMMETRYShowComputedMenu", 
 function () {
var menu = this.htMenus.get ("SYMMETRYShowComputedMenu");
if (menu == null) return;
this.menuRemoveAll (menu, 0);
this.menuEnable (menu, false);
if (!this.isSymmetry || this.modelIndex < 0) return;
var info = this.vwr.getProperty ("DATA_API", "spaceGroupInfo", null);
if (info == null) return;
var infolist = info.get ("operations");
if (infolist == null) return;
var name = info.get ("spaceGroupName");
this.menuSetLabel (menu, name == null ? J.i18n.GT._ ("Space Group") : name);
var subMenu = menu;
var nmod = this.itemMax;
var pt = (infolist.length > this.itemMax ? 0 : -2147483648);
for (var i = 0; i < infolist.length; i++) {
if (pt >= 0 && (pt++ % nmod) == 0) {
var id = "drawsymop" + pt + "Menu";
subMenu = this.menuNewSubMenu ((i + 1) + "..." + Math.min (i + this.itemMax, infolist.length), this.menuGetId (menu) + "." + id);
this.menuAddSubMenu (menu, subMenu);
this.htMenus.put (id, subMenu);
pt = 1;
}var sym = infolist[i][1];
if (sym.indexOf ("x1") < 0) sym = infolist[i][0];
var entryName = (i + 1) + " " + infolist[i][2] + " (" + sym + ")";
this.menuEnable (this.menuCreateItem (subMenu, entryName, "draw SYMOP " + (i + 1), null), true);
}
this.menuEnable (menu, true);
});
Clazz.defineMethod (c$, "updateSYMMETRYSelectComputedMenu", 
 function () {
var menu = this.htMenus.get ("SYMMETRYSelectComputedMenu");
if (menu == null) return;
this.menuRemoveAll (menu, 0);
this.menuEnable (menu, false);
if (!this.isSymmetry || this.modelIndex < 0) return;
var list = this.modelInfo.get ("symmetryOperations");
if (list == null) return;
var cellRange = this.modelInfo.get ("unitCellRange");
var haveUnitCellRange = (cellRange != null);
var subMenu = menu;
var nmod = this.itemMax;
var pt = (list.length > this.itemMax ? 0 : -2147483648);
for (var i = 0; i < list.length; i++) {
if (pt >= 0 && (pt++ % nmod) == 0) {
var id = "symop" + pt + "Menu";
subMenu = this.menuNewSubMenu ((i + 1) + "..." + Math.min (i + this.itemMax, list.length), this.menuGetId (menu) + "." + id);
this.menuAddSubMenu (menu, subMenu);
this.htMenus.put (id, subMenu);
pt = 1;
}var entryName = "symop=" + (i + 1) + " # " + list[i];
this.menuEnable (this.menuCreateItem (subMenu, entryName, "SELECT symop=" + (i + 1), null), haveUnitCellRange);
}
this.menuEnable (menu, true);
});
Clazz.defineMethod (c$, "updateFRAMESbyModelComputedMenu", 
 function () {
var menu = this.htMenus.get ("FRAMESbyModelComputedMenu");
if (menu == null) return;
this.menuEnable (menu, (this.modelCount > 0));
this.menuSetLabel (menu, (this.modelIndex < 0 ? this.gti ("allModelsText", this.modelCount) : this.gto ("modelMenuText", (this.modelIndex + 1) + "/" + this.modelCount)));
this.menuRemoveAll (menu, 0);
if (this.modelCount < 1) return;
if (this.modelCount > 1) this.menuCreateCheckboxItem (menu, J.i18n.GT._ ("All"), "frame 0 ##", null, (this.modelIndex < 0), false);
var subMenu = menu;
var nmod = this.itemMax;
var pt = (this.modelCount > this.itemMax ? 0 : -2147483648);
for (var i = 0; i < this.modelCount; i++) {
if (pt >= 0 && (pt++ % nmod) == 0) {
var id = "model" + pt + "Menu";
subMenu = this.menuNewSubMenu ((i + 1) + "..." + Math.min (i + this.itemMax, this.modelCount), this.menuGetId (menu) + "." + id);
this.menuAddSubMenu (menu, subMenu);
this.htMenus.put (id, subMenu);
pt = 1;
}var script = "" + this.vwr.getModelNumberDotted (i);
var entryName = this.vwr.getModelName (i);
var spectrumTypes = this.vwr.getModelAuxiliaryInfoValue (i, "spectrumTypes");
if (spectrumTypes != null && entryName.startsWith (spectrumTypes)) spectrumTypes = null;
if (!entryName.equals (script)) {
var ipt = entryName.indexOf (";PATH");
if (ipt >= 0) entryName = entryName.substring (0, ipt);
if (entryName.indexOf ("Model[") == 0 && (ipt = entryName.indexOf ("]:")) >= 0) entryName = entryName.substring (ipt + 2);
entryName = script + ": " + entryName;
}if (entryName.length > 60) entryName = entryName.substring (0, 55) + "...";
if (spectrumTypes != null) entryName += " (" + spectrumTypes + ")";
this.menuCreateCheckboxItem (subMenu, entryName, "model " + script + " ##", null, (this.modelIndex == i), false);
}
});
Clazz.defineMethod (c$, "updateConfigurationComputedMenu", 
 function () {
var menu = this.htMenus.get ("configurationComputedMenu");
if (menu == null) return;
this.menuEnable (menu, this.isMultiConfiguration);
if (!this.isMultiConfiguration) return;
var nAltLocs = this.altlocs.length;
this.menuSetLabel (menu, this.gti ("configurationMenuText", nAltLocs));
this.menuRemoveAll (menu, 0);
var script = "hide none ##CONFIG";
this.menuCreateCheckboxItem (menu, J.i18n.GT._ ("All"), script, null, (this.updateMode == 1 && this.configurationSelected.equals (script)), false);
for (var i = 0; i < nAltLocs; i++) {
script = "configuration " + (i + 1) + "; hide thisModel and not selected ##CONFIG";
var entryName = "" + (i + 1) + " -- \"" + this.altlocs.charAt (i) + "\"";
this.menuCreateCheckboxItem (menu, entryName, script, null, (this.updateMode == 1 && this.configurationSelected.equals (script)), false);
}
});
Clazz.defineMethod (c$, "updateModelSetComputedMenu", 
 function () {
var menu = this.htMenus.get ("modelSetMenu");
if (menu == null) return;
this.menuRemoveAll (menu, 0);
this.menuSetLabel (menu, this.nullModelSetName);
this.menuEnable (menu, false);
for (var i = this.noZapped.length; --i >= 0; ) this.menuEnable (this.htMenus.get (this.noZapped[i]), !this.isZapped);

if (this.modelSetName == null || this.isZapped) return;
if (this.isMultiFrame) {
this.modelSetName = this.gti ("modelSetCollectionText", this.modelCount);
if (this.modelSetName.length > this.titleWidthMax) this.modelSetName = this.modelSetName.substring (0, this.titleWidthMax) + "...";
} else if (this.vwr.getBooleanProperty ("hideNameInPopup")) {
this.modelSetName = this.getMenuText ("hiddenModelSetText");
} else if (this.modelSetName.length > this.titleWidthMax) {
this.modelSetName = this.modelSetName.substring (0, this.titleWidthMax) + "...";
}this.menuSetLabel (menu, this.modelSetName);
this.menuEnable (menu, true);
this.menuEnable (this.htMenus.get ("computationMenu"), this.ac <= 100);
this.addMenuItem (menu, this.gti ("atomsText", this.ac));
this.addMenuItem (menu, this.gti ("bondsText", this.vwr.getBondCountInModel (this.modelIndex)));
if (this.isPDB) {
this.menuAddSeparator (menu);
this.addMenuItem (menu, this.gti ("groupsText", this.vwr.getGroupCountInModel (this.modelIndex)));
this.addMenuItem (menu, this.gti ("chainsText", this.vwr.getChainCountInModel (this.modelIndex)));
this.addMenuItem (menu, this.gti ("polymersText", this.vwr.getPolymerCountInModel (this.modelIndex)));
var submenu = this.htMenus.get ("BiomoleculesMenu");
if (submenu == null) {
submenu = this.menuNewSubMenu (J.i18n.GT._ (this.getMenuText ("biomoleculesMenuText")), this.menuGetId (menu) + ".biomolecules");
this.menuAddSubMenu (menu, submenu);
}this.menuRemoveAll (submenu, 0);
this.menuEnable (submenu, false);
var biomolecules;
if (this.modelIndex >= 0 && (biomolecules = this.vwr.getModelAuxiliaryInfoValue (this.modelIndex, "biomolecules")) != null) {
this.menuEnable (submenu, true);
var nBiomolecules = biomolecules.size ();
for (var i = 0; i < nBiomolecules; i++) {
var script = (this.isMultiFrame ? "" : "save orientation;load \"\" FILTER \"biomolecule " + (i + 1) + "\";restore orientation;");
var nAtoms = (biomolecules.get (i).get ("atomCount")).intValue ();
var entryName = this.gto (this.isMultiFrame ? "biomoleculeText" : "loadBiomoleculeText", [Integer.$valueOf (i + 1), Integer.$valueOf (nAtoms)]);
this.menuCreateItem (submenu, entryName, script, null);
}
}}if (this.isApplet && !this.vwr.getBooleanProperty ("hideNameInPopup")) {
this.menuAddSeparator (menu);
this.menuCreateItem (menu, this.gto ("viewMenuText", this.modelSetFileName), "show url", null);
}});
Clazz.defineMethod (c$, "gti", 
 function (s, n) {
return J.i18n.GT.i (J.i18n.GT._ (this.getMenuText (s)), n);
}, "~S,~N");
Clazz.defineMethod (c$, "gto", 
 function (s, o) {
return J.i18n.GT.o (J.i18n.GT._ (this.getMenuText (s)), o);
}, "~S,~O");
Clazz.defineMethod (c$, "updateAboutSubmenu", 
 function () {
if (this.isApplet) this.setText ("APPLETid", this.vwr.appletName);
{
}});
Clazz.defineMethod (c$, "updateLanguageSubmenu", 
 function () {
var menu = this.htMenus.get ("languageComputedMenu");
if (menu == null) return;
this.menuRemoveAll (menu, 0);
var language = J.i18n.GT.getLanguage ();
var id = this.menuGetId (menu);
var languages = J.i18n.GT.getLanguageList (null);
for (var i = 0, p = 0; i < languages.length; i++) {
if (language.equals (languages[i].code)) languages[i].display = true;
if (languages[i].display) {
var code = languages[i].code;
var name = languages[i].language;
var nativeName = languages[i].nativeLanguage;
var menuLabel = code + " - " + J.i18n.GT._ (name);
if ((nativeName != null) && (!nativeName.equals (J.i18n.GT._ (name)))) {
menuLabel += " - " + nativeName;
}if (p++ > 0 && (p % 4 == 1)) this.menuAddSeparator (menu);
this.menuCreateCheckboxItem (menu, menuLabel, "language = \"" + code + "\" ##" + name, id + "." + code, language.equals (code), false);
}}
});
Clazz.defineMethod (c$, "updateSpecialMenuItem", 
 function (m) {
m.setText (this.getSpecialLabel (m.getName (), m.getText ()));
}, "javajs.api.SC");
Clazz.defineMethod (c$, "getSpecialLabel", 
function (name, text) {
var pt = text.indexOf (" (");
if (pt < 0) pt = text.length;
var info = null;
if (name.indexOf ("captureLooping") >= 0) info = (this.vwr.am.animationReplayMode.name ().equals ("ONCE") ? "ONCE" : "LOOP");
 else if (name.indexOf ("captureFps") >= 0) info = "" + this.vwr.getInt (553648132);
 else if (name.indexOf ("captureMenu") >= 0) info = (this.vwr.captureParams == null ? J.i18n.GT._ ("not capturing") : this.vwr.getFilePath (this.vwr.captureParams.get ("captureFileName"), true) + " " + this.vwr.captureParams.get ("captureCount"));
return (info == null ? text : text.substring (0, pt) + " (" + info + ")");
}, "~S,~S");
Clazz.defineStatics (c$,
"UPDATE_NEVER", -1,
"UPDATE_ALL", 0,
"UPDATE_CONFIG", 1,
"UPDATE_SHOW", 2,
"MENUITEM_HEIGHT", 20);
});
