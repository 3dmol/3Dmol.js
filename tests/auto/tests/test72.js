$3Dmol.download("pdb:4DM7",viewer,{},function(){
                  
                  viewer.setStyle({and:[{within:{distance:1,sel:{chain:'A'}}},{within:{distance:1,sel:{chain:'B'}}}]}
                  				  ,{line:{hidden:false,
                                    linewidth:1.0,
                                    colorscheme:'greenCarbon'}});

                 viewer.render();
                });