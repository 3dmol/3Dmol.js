$3Dmol.download("pdb:2EJ0",viewer,{},function(){
                  
                  viewer.addLabel("Aromatic", {position: {x:-6.89, y:0.75, z:0.35}, backgroundColor: 0x800080, backgroundOpacity: 0.8});
                  viewer.addLabel("Label",{font:'sans-serif',fontSize:18,fontColor:'white',fontOpacity:1,borderThickness:1.0,
                                           borderColor:'red',borderOpacity:0.5,backgroundColor:'black',backgroundOpacity:0.5,
                                           position:{x:50.0,y:0.0,z:0.0},inFront:true,showBackground:true});
                  viewer.setStyle({chain:'A'},{cross:{hidden:true}});
                  viewer.setStyle({chain:'B'},{cross:{hidden:false,
                                                      linewidth:1.0,
                                                      colorscheme:'greenCarbon'}});
                  viewer.setStyle({chain:'C'},{sphere:{hidden:false,
                                                      radius:0.5}});
                  viewer.setStyle({chain:'D'},{cross:{hidden:false}});
                  viewer.setStyle({chain:'E'},{cross:{hidden:false,
                                                      color:'black'}});
                  
                  viewer.render();

                  
                });
