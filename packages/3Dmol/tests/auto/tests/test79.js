
viewer.addCylinder({
    dashed:true,
    radius:.5,
    dashLength:2,
    gapLength:2,
    start:{x:0, y:0, z:0},
    end:{x:-10, y:10, z:10},
    fromCap:2,
    toCap:2,
    color:"green",
});

viewer.addCylinder({
    dashed:true,
    radius:.5,
    dashLength:2,
    gapLength:2,
    start:{x:5, y:5, z:5},
    end:{x:-10, y:10, z:10},
    fromCap:2,
    toCap:2,
    color:"blue",
});

viewer.addCylinder({
    dashed:true,
    radius:1,
    dashLength:2,
    gapLength:4,
    start:{x:-10, y:10, z:10},
    end:{x:-20, y:20, z:20},
    fromCap:2,
    toCap:2,
    color:"red",
});

viewer.render();
