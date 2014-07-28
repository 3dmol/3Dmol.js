Clazz.declarePackage ("J.awtjs2d");
Clazz.load (["javajs.api.GenericPlatform"], "J.awtjs2d.Platform", ["java.net.URL", "JU.AjaxURLStreamHandlerFactory", "J.api.Interface", "J.awtjs2d.Display", "$.Image", "$.JSFile", "$.JSFont", "$.Mouse"], function () {
c$ = Clazz.decorateAsClass (function () {
this.canvas = null;
this.vwr = null;
this.context = null;
Clazz.instantialize (this, arguments);
}, J.awtjs2d, "Platform", null, javajs.api.GenericPlatform);
Clazz.overrideMethod (c$, "setViewer", 
function (vwr, canvas) {
{
this.vwr = vwr;
this.canvas = canvas;
if (canvas != null) {
this.context = canvas.getContext("2d");
canvas.imgdata = this.context.getImageData(0, 0, canvas.width, canvas.height);
canvas.buf8 = canvas.imgdata.data;
}
}try {
java.net.URL.setURLStreamHandlerFactory ( new JU.AjaxURLStreamHandlerFactory ());
} catch (e) {
}
}, "javajs.api.PlatformViewer,~O");
Clazz.overrideMethod (c$, "isSingleThreaded", 
function () {
return true;
});
Clazz.overrideMethod (c$, "getJsObjectInfo", 
function (jsObject, method, args) {
{
if (method == "localName")return jsObject[0]["nodeName"];
return (args == null ? jsObject[0][method] : jsObject[0][method](args[0]));
}}, "~A,~S,~A");
Clazz.overrideMethod (c$, "isHeadless", 
function () {
return false;
});
Clazz.overrideMethod (c$, "getMouseManager", 
function (privateKey, display) {
return  new J.awtjs2d.Mouse (privateKey, this.vwr, display);
}, "~N,~O");
Clazz.overrideMethod (c$, "convertPointFromScreen", 
function (canvas, ptTemp) {
J.awtjs2d.Display.convertPointFromScreen (canvas, ptTemp);
}, "~O,JU.P3");
Clazz.overrideMethod (c$, "getFullScreenDimensions", 
function (canvas, widthHeight) {
J.awtjs2d.Display.getFullScreenDimensions (canvas, widthHeight);
}, "~O,~A");
Clazz.overrideMethod (c$, "getMenuPopup", 
function (menuStructure, type) {
var c = (type == 'j' ? "awtjs2d.JmolJSPopup" : "awtjs2d.JSModelKitPopup");
var jmolpopup = J.api.Interface.getOption (c);
try {
if (jmolpopup != null) jmolpopup.jpiInitialize (this.vwr, menuStructure);
} catch (e) {
if (Clazz.exceptionOf (e, Exception)) {
c = "Exception creating " + c + ":" + e;
System.out.println (c);
return null;
} else {
throw e;
}
}
return jmolpopup;
}, "~S,~S");
Clazz.overrideMethod (c$, "hasFocus", 
function (canvas) {
return J.awtjs2d.Display.hasFocus (canvas);
}, "~O");
Clazz.overrideMethod (c$, "prompt", 
function (label, data, list, asButtons) {
return J.awtjs2d.Display.prompt (label, data, list, asButtons);
}, "~S,~S,~A,~B");
Clazz.overrideMethod (c$, "renderScreenImage", 
function (context, size) {
J.awtjs2d.Display.renderScreenImage (this.vwr, context, size);
}, "~O,~O");
Clazz.overrideMethod (c$, "drawImage", 
function (context, canvas, x, y, width, height) {
J.awtjs2d.Image.drawImage (context, canvas, x, y, width, height);
}, "~O,~O,~N,~N,~N,~N");
Clazz.overrideMethod (c$, "requestFocusInWindow", 
function (canvas) {
J.awtjs2d.Display.requestFocusInWindow (canvas);
}, "~O");
Clazz.overrideMethod (c$, "repaint", 
function (canvas) {
J.awtjs2d.Display.repaint (canvas);
}, "~O");
Clazz.overrideMethod (c$, "setTransparentCursor", 
function (canvas) {
J.awtjs2d.Display.setTransparentCursor (canvas);
}, "~O");
Clazz.overrideMethod (c$, "setCursor", 
function (c, canvas) {
J.awtjs2d.Display.setCursor (c, canvas);
}, "~N,~O");
Clazz.overrideMethod (c$, "allocateRgbImage", 
function (windowWidth, windowHeight, pBuffer, windowSize, backgroundTransparent, isImageWrite) {
return J.awtjs2d.Image.allocateRgbImage (windowWidth, windowHeight, pBuffer, windowSize, backgroundTransparent, (isImageWrite ? null : this.canvas));
}, "~N,~N,~A,~N,~B,~B");
Clazz.overrideMethod (c$, "notifyEndOfRendering", 
function () {
});
Clazz.overrideMethod (c$, "createImage", 
function (data) {
return null;
}, "~O");
Clazz.overrideMethod (c$, "disposeGraphics", 
function (gOffscreen) {
}, "~O");
Clazz.overrideMethod (c$, "grabPixels", 
function (canvas, width, height, pixels, startRow, nRows) {
{
if (canvas.image && (width != canvas.width || height != canvas.height))
Jmol._setCanvasImage(canvas, width, height);
if (canvas.buf32) return canvas.buf32;
}var buf = J.awtjs2d.Image.grabPixels (J.awtjs2d.Image.getGraphics (canvas), width, height);
{
canvas.buf32 = buf;
}return buf;
}, "~O,~N,~N,~A,~N,~N");
Clazz.overrideMethod (c$, "drawImageToBuffer", 
function (gOffscreen, imageOffscreen, canvas, width, height, bgcolor) {
return this.grabPixels (canvas, width, height, null, 0, 0);
}, "~O,~O,~O,~N,~N,~N");
Clazz.overrideMethod (c$, "getTextPixels", 
function (text, font3d, context, image, width, height, ascent) {
return J.awtjs2d.Image.getTextPixels (text, font3d, context, width, height, ascent);
}, "~S,javajs.awt.Font,~O,~O,~N,~N,~N");
Clazz.overrideMethod (c$, "flushImage", 
function (imagePixelBuffer) {
}, "~O");
Clazz.overrideMethod (c$, "getGraphics", 
function (image) {
return (image == null ? this.context : J.awtjs2d.Image.getGraphics (image));
}, "~O");
Clazz.overrideMethod (c$, "getImageHeight", 
function (canvas) {
return (canvas == null ? -1 : J.awtjs2d.Image.getHeight (canvas));
}, "~O");
Clazz.overrideMethod (c$, "getImageWidth", 
function (canvas) {
return (canvas == null ? -1 : J.awtjs2d.Image.getWidth (canvas));
}, "~O");
Clazz.overrideMethod (c$, "getStaticGraphics", 
function (image, backgroundTransparent) {
return J.awtjs2d.Image.getStaticGraphics (image, backgroundTransparent);
}, "~O,~B");
Clazz.overrideMethod (c$, "newBufferedImage", 
function (image, w, h) {
{
if (typeof Jmol != "undefined" && Jmol._getHiddenCanvas)
return Jmol._getHiddenCanvas(this.vwr.applet, "stereoImage", w, h);
}return null;
}, "~O,~N,~N");
Clazz.overrideMethod (c$, "newOffScreenImage", 
function (w, h) {
{
if (typeof Jmol != "undefined" && Jmol._getHiddenCanvas)
return Jmol._getHiddenCanvas(this.vwr.applet, "textImage", w, h);
}return null;
}, "~N,~N");
Clazz.overrideMethod (c$, "waitForDisplay", 
function (echoNameAndPath, zipBytes) {
{
if (typeof Jmol == "undefined" || !Jmol._getHiddenCanvas) return false;
var vwr = this.vwr;
var sc = vwr.getEvalContextAndHoldQueue(vwr.eval);
var echoName = echoNameAndPath[0];
return Jmol._loadImage(this, echoNameAndPath, zipBytes,
function(canvas, pathOrError) { vwr.loadImageData(canvas, pathOrError, echoName, sc) }
);
}}, "~O,~O");
Clazz.overrideMethod (c$, "fontStringWidth", 
function (font, text) {
return J.awtjs2d.JSFont.stringWidth (font, this.context, text);
}, "javajs.awt.Font,~S");
Clazz.overrideMethod (c$, "getFontAscent", 
function (context) {
return J.awtjs2d.JSFont.getAscent (context);
}, "~O");
Clazz.overrideMethod (c$, "getFontDescent", 
function (context) {
return J.awtjs2d.JSFont.getDescent (context);
}, "~O");
Clazz.overrideMethod (c$, "getFontMetrics", 
function (font, context) {
return J.awtjs2d.JSFont.getFontMetrics (font, context);
}, "javajs.awt.Font,~O");
Clazz.overrideMethod (c$, "newFont", 
function (fontFace, isBold, isItalic, fontSize) {
return J.awtjs2d.JSFont.newFont (fontFace, isBold, isItalic, fontSize, "px");
}, "~S,~B,~B,~N");
Clazz.overrideMethod (c$, "getDateFormat", 
function (isoiec8824) {
{
if (isoiec8824) {
var d = new Date();
var x = d.toString().split(" ");
var MM = "0" + d.getMonth(); MM = MM.substring(MM.length - 2);
var dd = "0" + d.getDate(); dd = dd.substring(dd.length - 2);
return x[3] + MM + dd + x[4].replace(/\:/g,"") + x[5].substring(3,6) + "'" + x[5].substring(6,8) + "'"
}
return ("" + (new Date())).split(" (")[0];
}}, "~B");
Clazz.overrideMethod (c$, "newFile", 
function (name) {
return  new J.awtjs2d.JSFile (name);
}, "~S");
Clazz.overrideMethod (c$, "getBufferedFileInputStream", 
function (name) {
return null;
}, "~S");
Clazz.overrideMethod (c$, "getBufferedURLInputStream", 
function (url, outputBytes, post) {
return J.awtjs2d.JSFile.getBufferedURLInputStream (url, outputBytes, post);
}, "java.net.URL,~A,~S");
Clazz.overrideMethod (c$, "getLocalUrl", 
function (fileName) {
return null;
}, "~S");
});
