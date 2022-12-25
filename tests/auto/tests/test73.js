$3Dmol.download("pdb:4DM7",viewer,{},function(){
                  viewer.setStyle({or:[{chain:'C'},{chain:'D'}]},{line:{hidden:false,
                                                    colorscheme:'greenCarbon'}});
                  viewer.setStyle({or:[{chain:'A'},{chain:'B'}]},{line:{hidden:false,
                                                    colorscheme:'magentaCarbon'}});
                  
                  viewer.render();
                });
