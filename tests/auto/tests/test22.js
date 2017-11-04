$3Dmol.download("pdb:4UAA",viewer,{},function(){
                  
                  viewer.setStyle({chain:'A'},{cartoon:{opacity:0.5,colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.Sinebow($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
                  viewer.setStyle({chain:'B'},{line:{colorscheme:{prop:'partialCharge',gradient: new $3Dmol.Gradient.Sinebow($3Dmol.getPropertyRange(viewer.selectedAtoms(),'partialCharge'))}}});
                  viewer.render();
              });
