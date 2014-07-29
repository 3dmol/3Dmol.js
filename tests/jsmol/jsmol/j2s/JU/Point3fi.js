Clazz.declarePackage ("JU");
Clazz.load (["JU.P3"], "JU.Point3fi", null, function () {
c$ = Clazz.decorateAsClass (function () {
this.mi = -1;
this.i = 0;
this.sX = 0;
this.sY = 0;
this.sZ = 0;
this.sD = -1;
Clazz.instantialize (this, arguments);
}, JU, "Point3fi", JU.P3, Cloneable);
c$.set2 = Clazz.defineMethod (c$, "set2", 
function (p3f, p3i) {
p3f.x = p3i.x;
p3f.y = p3i.y;
p3f.z = p3i.z;
}, "JU.P3,JU.P3i");
Clazz.defineMethod (c$, "copy", 
function () {
try {
return this.clone ();
} catch (e) {
if (Clazz.exceptionOf (e, CloneNotSupportedException)) {
return null;
} else {
throw e;
}
}
});
});
