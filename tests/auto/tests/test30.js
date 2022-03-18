$3Dmol.download("pdb:4UAA",viewer,{},()=> {
                  
                  const colorAsSnake = function(atom) {
                    return atom.resi % 2 ? 'white': 'green'
                  };

                  viewer.setStyle( {}, { cartoon: {colorfunc: colorAsSnake }});

                  viewer.render();
              });
