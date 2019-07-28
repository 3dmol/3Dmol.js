
        var shape = viewer.addShape({color:'red'});
        shape.addBox({corner: {x:1,y:2,z:0}, dimensions: {w: 4, h: 2, d: 6}});
        shape.addBox({corner: {x:-5,y:-3,z:0},
                      dimensions: { w: {x:1,y:1,z:0},
                                    h: {x:-1,y:1,z:0},
                                    d: {x:0,y:0,z:1} }});
        shape.addCylinder({start:{x:0.0,y:2.0,z:0.0},
                                        end:{x:0.0,y:10.0,z:0.0},
                                        radius:0.5,
                                        fromCap:false,
                                        toCap:true,
                                        color:'teal'});
        var state = shape.getInternalState();
        viewer.removeAllShapes();
        var shape = viewer.addShape();
        shape.setInternalState(state);
        viewer.zoomTo();
        viewer.render();
    

      
   
  