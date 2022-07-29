 

viewer.addCylinder({
    radius:.5,
    start:{x:0, y:0, z:0},
    end:{x:-10, y:10, z:10},
    fromCap:$3Dmol.CAP.FLAT,
    toCap: $3Dmol.CAP.ROUND,
    color:"green",
});
viewer.zoomTo();
viewer.render();
