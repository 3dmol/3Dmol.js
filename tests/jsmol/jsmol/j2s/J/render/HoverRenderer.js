Clazz.declarePackage ("J.render");
Clazz.load (["J.render.ShapeRenderer"], "J.render.HoverRenderer", ["JU.P3", "J.render.TextRenderer"], function () {
c$ = Clazz.decorateAsClass (function () {
this.tempXY = null;
this.ptTemp = null;
Clazz.instantialize (this, arguments);
}, J.render, "HoverRenderer", J.render.ShapeRenderer);
Clazz.prepareFields (c$, function () {
this.tempXY =  Clazz.newFloatArray (3, 0);
});
Clazz.overrideMethod (c$, "render", 
function () {
if (this.tm.isNavigating ()) return false;
if (this.ptTemp == null) this.ptTemp =  new JU.P3 ();
var hover = this.shape;
var antialias = this.g3d.isAntialiased ();
var text = hover.hoverText;
if (hover.atomIndex >= 0) {
var atom = this.ms.at[hover.atomIndex];
var label = (hover.specialLabel != null ? hover.specialLabel : hover.atomFormats != null && hover.atomFormats[hover.atomIndex] != null ? this.vwr.ms.getLabeler ().formatLabel (this.vwr, atom, hover.atomFormats[hover.atomIndex], this.ptTemp) : hover.labelFormat != null ? this.vwr.ms.getLabeler ().formatLabel (this.vwr, atom, this.fixLabel (atom, hover.labelFormat), this.ptTemp) : null);
if (label == null) return false;
text.setText (label);
text.setXYZs (atom.sX, atom.sY, 1, -2147483648);
} else if (hover.text != null) {
text.setText (hover.text);
text.setXYZs (hover.xy.x, hover.xy.y, 1, -2147483648);
} else {
return true;
}J.render.TextRenderer.render (text, this.vwr, this.g3d, 0, antialias ? 2 : 1, false, null, this.tempXY);
return true;
});
Clazz.defineMethod (c$, "fixLabel", 
function (atom, label) {
if (label == null) return null;
return (this.vwr.ms.isJmolDataFrameForModel (atom.getModelIndex ()) && label.equals ("%U") ? "%W" : label);
}, "JM.Atom,~S");
});
