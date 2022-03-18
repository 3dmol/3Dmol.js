
              $3Dmol.download("pdb:4UB9",viewer,{},()=> {
                  
                  
                  const atoms = viewer.selectedAtoms();
                  for(let i = 0; i < atoms.length; i++) {
                    const a = atoms[i];
                    a.properties.structured = (a.ss == 'h' || a.ss == 's');                    
                  }
                  
                  viewer.setStyle({properties: {structured: true}}, {cartoon: {}});
                  viewer.render();
              });
