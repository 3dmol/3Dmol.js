Clazz.declarePackage ("J.shapesurface");
Clazz.load (["J.shapesurface.Isosurface"], "J.shapesurface.Pmesh", null, function () {
c$ = Clazz.declareType (J.shapesurface, "Pmesh", J.shapesurface.Isosurface);
Clazz.defineMethod (c$, "initShape", 
function () {
Clazz.superCall (this, J.shapesurface.Pmesh, "initShape", []);
this.myType = "pmesh";
});
});
