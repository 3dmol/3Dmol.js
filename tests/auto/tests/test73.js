$3Dmol.download("pdb:4DM7",viewer,{},function(){
                  viewer.setStyle({or:[{chain:'C'},{chain:'D'}]},{line:{hidden:false,
                     		                         linewidth:1.0,
                                                    colorscheme:'greenCarbon'}});
                  viewer.setStyle({or:[{chain:'A'},{chain:'B'}]},{line:{hidden:false,
                                                    linewidth:1.0,
                                                    colorscheme:'magentaCarbon'}});
                  
                  viewer.render();
                });
