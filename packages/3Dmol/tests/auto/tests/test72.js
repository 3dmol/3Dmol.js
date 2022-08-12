$3Dmol.download("pdb:4DM7",viewer,{},function(){
                  viewer.setStyle({stick:{radius:.1,colorscheme:'rasmol'}});
                  viewer.setStyle({byres:true, and:[{within:{distance:5,sel:{chain:'A'}}},{within:{distance:5,sel:{chain:'B'}}}]}
                  				  ,{stick:{hidden:false,
                                    colorscheme:'greenCarbon'}});

                 viewer.render();
                });
