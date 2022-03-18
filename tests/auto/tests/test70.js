$3Dmol.download("pdb:4DM7",viewer,{},()=> {
                  
                  
                  viewer.setStyle({not:{chain:'A'}},{line:{hidden:false,
                                                    colorscheme:'magentaCarbon'}});

                 viewer.render();
                });
