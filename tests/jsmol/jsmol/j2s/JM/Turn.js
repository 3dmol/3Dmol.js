Clazz.declarePackage ("JM");
Clazz.load (["JM.ProteinStructure"], "JM.Turn", ["J.c.STR"], function () {
c$ = Clazz.declareType (JM, "Turn", JM.ProteinStructure);
Clazz.makeConstructor (c$, 
function (apolymer, monomerIndex, monomerCount) {
Clazz.superConstructor (this, JM.Turn, []);
this.setupPS (apolymer, J.c.STR.TURN, monomerIndex, monomerCount);
this.subtype = J.c.STR.TURN;
}, "JM.AlphaPolymer,~N,~N");
});
