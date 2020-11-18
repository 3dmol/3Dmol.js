var v = viewer;
              viewer.addCylinder({start:{x:0.0,y:0.0,z:0.0},
                                  end:{x:10.0,y:0.0,z:0.0},
                                  radius:1.0,
                                  fromCap:1,
                                  toCap:2,
                                  color:'red',
                                  hoverable:true,
                                  clickable:true,
                                  callback:function(){ this.color.setHex(0x00FFFF00);v.render();},
                                  hover_callback: function(){ v.render();},
                                  unhover_callback: function(){ this.color.setHex(0xFF0000);v.render();}
                                 });
              viewer.addCylinder({start:{x:0.0,y:2.0,z:0.0},
                                  end:{x:0.0,y:10.0,z:-15.0},
                                  radius:0.5,
                                  fromCap:false,
                                  toCap:true,
                                  color:'teal'});
              viewer.addCylinder({start:{x:15.0,y:0.0,z:0.0},
                                  end:{x:15.0,y:0.0,z:10.0},
                                  radius:1.0,
                                  color:'black',
                                  fromCap:false,
                                  toCap:false});
                                  
              viewer.zoomTo();
              viewer.render();
