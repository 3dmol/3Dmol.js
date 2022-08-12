


viewer.addLine({
    dashed:true,
    linewidth:1,
    dashLength:0.25,
    gapLength:0.25,
    start:{x:0, y:0, z:0},
    end:{x:100, y:100, z:100}
});


viewer.addLine({
    dashed:true,
    linewidth:2,
    dashLength:0.25,
    gapLength:0.25,
    start:{x:0, y:0, z:0},
    end:{x:-100, y:100, z:100}
});

viewer.addLine({
    dashed:true,
    linewidth:4,
    dashLength:0.25,
    gapLength:0.25,
    start:{x:0, y:0, z:0},
    end:{x:100, y:-100, z:100}
});


viewer.setBackgroundColor(0xffffff);

viewer.render()
          
