$3Dmol.download("pdb:5BP0",viewer,{},function(){
                  
                  viewer.setStyle({chain:'A'},{sphere:{colorscheme:{prop:'partialCharge',gradient: new $3Dmol.Gradient.Sinebow($3Dmol.getPropertyRange(viewer.selectedAtoms(),'partialCharge'))}}});
                  viewer.setStyle({chain:'B'},{line:{linewidth:2.0,colorscheme:{prop:'partialCharge',gradient: new $3Dmol.Gradient.RWB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'partialCharge'))}}});
                  viewer.setStyle({chain:'C'},{cross:{hidden:true,colorscheme:{prop:'partialCharge',gradient: new $3Dmol.Gradient.ROYGB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'partialCharge'))}}});
                  viewer.setStyle({chain:'D'},{cross:{linewidth:2.0,colorscheme:{prop:'partialCharge',gradient: new $3Dmol.Gradient.Sinebow($3Dmol.getPropertyRange(viewer.selectedAtoms(),'partialCharge'))}}});
                  viewer.setStyle({chain:'E'},{cross:{radius:2.0,colorscheme:{prop:'partialCharge',gradient: new $3Dmol.Gradient.RWB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'partialCharge'))}}});
                  viewer.setStyle({chain:'F'},{stick:{hidden:true,colorscheme:{prop:'partialCharge',gradient: new $3Dmol.Gradient.ROYGB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'partialCharge'))}}});
                  viewer.setStyle({chain:'G'},{stick:{radius:.2,colorscheme:{prop:'partialCharge',gradient: new $3Dmol.Gradient.Sinebow($3Dmol.getPropertyRange(viewer.selectedAtoms(),'partialCharge'))}}});
                  viewer.setStyle({chain:'H'},{stick:{singleBonds:true,colorscheme:{prop:'partialCharge',gradient: new $3Dmol.Gradient.RWB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'partialCharge'))}}});
                  viewer.setStyle({chain:'I'},{sphere:{hidden:true,colorscheme:{prop:'partialCharge',gradient: new $3Dmol.Gradient.ROYGB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'partialCharge'))}}});
                  viewer.setStyle({chain:'J'},{sphere:{radius:3.0,colorscheme:{prop:'partialCharge',gradient: new $3Dmol.Gradient.Sinebow($3Dmol.getPropertyRange(viewer.selectedAtoms(),'partialCharge'))}}});
                  viewer.render();
              });
