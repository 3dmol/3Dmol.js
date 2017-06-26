
          $3Dmol.download("pdb:4UAA",viewer,{},function(){  
            viewer.setStyle({},{stick:{}});

                viewer.addSphere({center:{},radius:10.0,color:'red'});
                viewer.addCylinder({start:{},
                                  end:{x:20.0,y:0.0,z:0.0},
                                  radius:1.0,
                                  color:'black',
                                  fromCap:false,
                                  toCap:false});
                                   viewer.addLine({dashed:true,start:{x:0,y:0,z:0},end:{}});

                viewer.addArrow({
                      start: {},
                      end: {x:-10.0, y:0.0, z:0.0},
                      radius: 1.0,
                      radiusRadio:1.0,
                      mid:1.0,
                      clickable:true,
                      callback:function(){
                          this.color.setHex(0xFF0000FF);
                          viewer.render( );
                      }
                  });
                
            viewer.zoomTo();
            viewer.render();
              });   