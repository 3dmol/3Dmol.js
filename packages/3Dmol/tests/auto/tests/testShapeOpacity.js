 

viewer.addCylinder({
    radius:.5,
    start:{x:0, y:0, z:0},
    end:{x:-10, y:10, z:10},
    fromCap:$3Dmol.CAP.FLAT,
    toCap: $3Dmol.CAP.ROUND,
    color:"green",
    alpha: 0.5
});

viewer.addCurve({points: [{x:0.0,y:0.0,z:0.0}, {x:5.0,y:3.0,z:0.0}, {x:5.0,y:7.0,z:0.0}, {x:0.0,y:10.0,z:0.0}],
          radius:0.5,
          smooth: 10,
          fromArrow:false,
          toArrow: true,
          color:'orange',                                  
          alpha: 0.5
          });
          
viewer.addBox({center:{x:0,y:0,z:0},dimensions: {w:3,h:4,d:2},color:'magenta'});
viewer.addSphere({center:{x:0,y:0,z:2},radius:1.0,color:'red',alpha:0.4});

viewer.addArrow({
                      start: {x:-10.0, y:0.0, z:0.0},
                      end: {x:0.0, y:-10.0, z:0.0},
                      radius: 1.0,
                      radiusRadio:1.0,
                      mid:1.0,
                      alpha: 0.7
                  });
                                  
viewer.zoomTo();
viewer.render();
