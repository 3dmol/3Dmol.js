$3Dmol.download("pdb:4DM7",viewer,{},function(){
                  
                  viewer.addArrow({
                      start: {x:-10.0, y:0.0, z:0.0},
                      end: {x:0.0, y:-10.0, z:0.0},
                      radius: 1.0,
                      radiusRadio:1.0,
                      mid:1.0,
                      clickable:true,
                      callback:function(){
                          this.color.setHex(0xFF0000FF);
                          viewer.render( /*no callback*/);
                      }
                  });
                  viewer.setStyle({chain:'A'},{line:{hidden:true}});
                  viewer.setStyle({chain:'B'},{line:{hidden:false,
                                                     colorscheme:'greenCarbon'}});
                  viewer.setStyle({chain:'C'},{line:{hidden:false,
                                                     colorscheme:'whiteCarbon'}});
                  viewer.setStyle({chain:'D'},{line:{hidden:false,
                                                     color:'black'}});

                  viewer.center({chain:'B'});
                  viewer.render();
                });
