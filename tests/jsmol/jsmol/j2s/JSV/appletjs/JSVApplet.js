Clazz.declarePackage ("JSV.appletjs");
Clazz.load (["JSV.api.AppletFrame", "$.JSVAppletInterface"], "JSV.appletjs.JSVApplet", ["java.lang.Boolean", "java.net.URL", "java.util.Hashtable", "JU.PT", "JSV.app.JSVApp", "JSV.common.JSVersion", "JSV.js2d.JsMainPanel", "$.JsPanel", "JU.Logger"], function () {
c$ = Clazz.decorateAsClass (function () {
this.app = null;
this.viewer = null;
this.isStandalone = false;
this.viewerOptions = null;
this.htParams = null;
Clazz.instantialize (this, arguments);
}, JSV.appletjs, "JSVApplet", null, [JSV.api.JSVAppletInterface, JSV.api.AppletFrame]);
Clazz.makeConstructor (c$, 
function (viewerOptions) {
if (viewerOptions == null) viewerOptions =  new java.util.Hashtable ();
this.viewerOptions = viewerOptions;
this.htParams =  new java.util.Hashtable ();
for (var entry, $entry = viewerOptions.entrySet ().iterator (); $entry.hasNext () && ((entry = $entry.next ()) || true);) this.htParams.put (entry.getKey ().toLowerCase (), entry.getValue ());

this.init ();
}, "java.util.Map");
Clazz.defineMethod (c$, "init", 
function () {
this.app =  new JSV.app.JSVApp (this, true);
this.initViewer ();
if (this.app.appletReadyCallbackFunctionName != null && this.viewer.fullName != null) this.callToJavaScript (this.app.appletReadyCallbackFunctionName, [this.viewer.appletID, this.viewer.fullName, Boolean.TRUE, this]);
});
Clazz.defineMethod (c$, "initViewer", 
function () {
this.viewer = this.app.vwr;
this.setLogging ();
this.viewerOptions.remove ("debug");
var o = this.viewerOptions.get ("display");
{
o = document.getElementById(o);
}this.viewer.setDisplay (o);
JU.Logger.info (this.getAppletInfo ());
});
Clazz.defineMethod (c$, "setLogging", 
 function () {
var iLevel = (this.getValue ("logLevel", (this.getBooleanValue ("debug", false) ? "5" : "4"))).charCodeAt (0) - 48;
if (iLevel != 4) System.out.println ("setting logLevel=" + iLevel + " -- To change, use script \"set logLevel [0-5]\"");
JU.Logger.setLogLevel (iLevel);
});
Clazz.defineMethod (c$, "getParameter", 
function (paramName) {
var o = this.htParams.get (paramName.toLowerCase ());
return (o == null ? null :  String.instantialize (o.toString ()));
}, "~S");
Clazz.defineMethod (c$, "getBooleanValue", 
 function (propertyName, defaultValue) {
var value = this.getValue (propertyName, defaultValue ? "true" : "");
return (value.equalsIgnoreCase ("true") || value.equalsIgnoreCase ("on") || value.equalsIgnoreCase ("yes"));
}, "~S,~B");
Clazz.defineMethod (c$, "getValue", 
 function (propertyName, defaultValue) {
var stringValue = this.getParameter (propertyName);
System.out.println ("getValue " + propertyName + " = " + stringValue);
if (stringValue != null) return stringValue;
return defaultValue;
}, "~S,~S");
Clazz.overrideMethod (c$, "isPro", 
function () {
return this.app.isPro ();
});
Clazz.overrideMethod (c$, "isSigned", 
function () {
return this.app.isSigned ();
});
Clazz.overrideMethod (c$, "finalize", 
function () {
System.out.println ("JSpecView " + this + " finalized");
});
Clazz.defineMethod (c$, "destroy", 
function () {
this.app.dispose ();
this.app = null;
});
Clazz.defineMethod (c$, "getParameter", 
function (key, def) {
return this.isStandalone ? System.getProperty (key, def) : (this.getParameter (key) != null ? this.getParameter (key) : def);
}, "~S,~S");
Clazz.overrideMethod (c$, "getAppletInfo", 
function () {
return "JSpecView Applet " + JSV.common.JSVersion.VERSION;
});
Clazz.overrideMethod (c$, "getSolnColour", 
function () {
return this.app.getSolnColour ();
});
Clazz.overrideMethod (c$, "getCoordinate", 
function () {
return this.app.getCoordinate ();
});
Clazz.overrideMethod (c$, "loadInline", 
function (data) {
this.app.loadInline (data);
}, "~S");
Clazz.defineMethod (c$, "$export", 
function (type, n) {
return this.app.exportSpectrum (type, n);
}, "~S,~N");
Clazz.overrideMethod (c$, "exportSpectrum", 
function (type, n) {
return this.app.exportSpectrum (type, n);
}, "~S,~N");
Clazz.overrideMethod (c$, "setFilePath", 
function (tmpFilePath) {
this.app.setFilePath (tmpFilePath);
}, "~S");
Clazz.overrideMethod (c$, "setSpectrumNumber", 
function (i) {
this.app.setSpectrumNumber (i);
}, "~N");
Clazz.overrideMethod (c$, "toggleGrid", 
function () {
this.app.toggleGrid ();
});
Clazz.overrideMethod (c$, "toggleCoordinate", 
function () {
this.app.toggleCoordinate ();
});
Clazz.overrideMethod (c$, "toggleIntegration", 
function () {
this.app.toggleIntegration ();
});
Clazz.overrideMethod (c$, "addHighlight", 
function (x1, x2, r, g, b, a) {
this.app.addHighlight (x1, x2, r, g, b, a);
}, "~N,~N,~N,~N,~N,~N");
Clazz.overrideMethod (c$, "removeAllHighlights", 
function () {
this.app.removeAllHighlights ();
});
Clazz.overrideMethod (c$, "removeHighlight", 
function (x1, x2) {
this.app.removeHighlight (x1, x2);
}, "~N,~N");
Clazz.overrideMethod (c$, "reversePlot", 
function () {
this.app.reversePlot ();
});
Clazz.defineMethod (c$, "script", 
function (script) {
this.app.initParams (script);
}, "~S");
Clazz.overrideMethod (c$, "runScript", 
function (script) {
this.app.runScript (script);
}, "~S");
Clazz.overrideMethod (c$, "syncScript", 
function (peakScript) {
this.app.syncScript (peakScript);
}, "~S");
Clazz.overrideMethod (c$, "writeStatus", 
function (msg) {
this.app.writeStatus (msg);
}, "~S");
Clazz.overrideMethod (c$, "getPropertyAsJavaObject", 
function (key) {
return this.app.getPropertyAsJavaObject (key);
}, "~S");
Clazz.overrideMethod (c$, "getPropertyAsJSON", 
function (key) {
return this.app.getPropertyAsJSON (key);
}, "~S");
Clazz.overrideMethod (c$, "runScriptNow", 
function (script) {
return this.app.runScriptNow (script);
}, "~S");
Clazz.overrideMethod (c$, "print", 
function (fileName) {
return this.app.print (fileName);
}, "~S");
Clazz.overrideMethod (c$, "setDropTargetListener", 
function (isSigned, viewer) {
}, "~B,JSV.common.JSViewer");
Clazz.overrideMethod (c$, "validateContent", 
function (mode) {
}, "~N");
Clazz.overrideMethod (c$, "createMainPanel", 
function (viewer) {
viewer.mainPanel =  new JSV.js2d.JsMainPanel ();
}, "JSV.common.JSViewer");
Clazz.overrideMethod (c$, "newWindow", 
function (isSelected) {
}, "~B");
Clazz.overrideMethod (c$, "callToJavaScript", 
function (callback, data) {
var tokens = JU.PT.split (callback, ".");
{
try{
var o = window[tokens[0]]
for (var i = 1; i < tokens.length; i++){
o = o[tokens[i]]
}
return o(data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9]);
} catch (e) {
System.out.println(callback + " failed " + e);
}
}}, "~S,~A");
Clazz.overrideMethod (c$, "setPanelVisible", 
function (b) {
}, "~B");
Clazz.overrideMethod (c$, "getJSVPanel", 
function (viewer, specs) {
return (specs == null ? JSV.js2d.JsPanel.getEmptyPanel (viewer) : JSV.js2d.JsPanel.getPanelMany (viewer, specs));
}, "JSV.common.JSViewer,JU.Lst");
Clazz.overrideMethod (c$, "setVisible", 
function (b) {
}, "~B");
Clazz.overrideMethod (c$, "getDocumentBase", 
function () {
try {
return  new java.net.URL (Clazz.castNullAs ("java.net.URL"), this.viewerOptions.get ("documentBase"), null);
} catch (e) {
if (Clazz.exceptionOf (e, java.net.MalformedURLException)) {
return null;
} else {
throw e;
}
}
});
Clazz.overrideMethod (c$, "repaint", 
function () {
});
Clazz.overrideMethod (c$, "validate", 
function () {
});
Clazz.overrideMethod (c$, "doExitJmol", 
function () {
});
Clazz.overrideMethod (c$, "getApp", 
function () {
return this.app;
});
});
