
              $3Dmol.download("pdb:4UB9",viewer,{},function(){
                  
                  
                  var atoms = viewer.selectedAtoms();
                  for(var i = 0; i < atoms.length; i++) {
                    var a = atoms[i];
                    a.properties.structured = (a.ss == 'h' || a.ss == 's');                    
                  }
                  
                  viewer.setStyle({properties: {structured: true}}, {cartoon: {}});
                  viewer.render();
              });
