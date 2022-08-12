 

viewer.addCylinder({
    radius:.5,
    start:{x:0, y:0, z:0},
    end:{x:-10, y:10, z:10},
    fromCap:$3Dmol.CAP.FLAT,
    toCap: $3Dmol.CAP.ROUND,
    color:"green",
    wireframe: true
});

viewer.addCurve({points: [{x:0.0,y:0.0,z:0.0}, {x:5.0,y:3.0,z:0.0}, {x:5.0,y:7.0,z:0.0}, {x:0.0,y:10.0,z:0.0}],
          radius:0.5,
          smooth: 10,
          fromArrow:false,
          toArrow: true,
          color:'orange',                                  
          wireframe: true
          });
          
viewer.addBox({center:{x:0,y:0,z:0},dimensions: {w:3,h:4,d:2},color:'magenta',wireframe:true});
viewer.addSphere({center:{x:0,y:0,z:0},radius:10.0,color:'red',wireframe:true});

viewer.addArrow({
                      start: {x:-10.0, y:0.0, z:0.0},
                      end: {x:0.0, y:-10.0, z:0.0},
                      radius: 1.0,
                      radiusRadio:1.0,
                      mid:1.0,
                      wireframe: true
                  });
                                  
viewer.zoomTo();
viewer.render();
