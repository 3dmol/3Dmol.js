 

viewer.addCylinder({
    radius:.5,
    start:{x:0, y:0, z:0},
    end:{x:-10, y:10, z:10},
    fromCap:1,
    toCap:2,
    color:"green",
});
viewer.zoomTo();
viewer.render();
