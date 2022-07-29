$3Dmol.download("pdb:1ejg",viewer,{altLoc:'*'},function(){
                  
                  viewer.setStyle({stick:{}});
                  viewer.setStyle({altLoc:'B'},{stick:{colorscheme:'greenCarbon'}});
                  viewer.setStyle({altLoc:'C'},{stick:{colorscheme:'cyanCarbon'}});
                  viewer.render();
              });
