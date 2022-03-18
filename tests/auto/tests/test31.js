$3Dmol.download("pdb:1BNA", viewer, {}, ()=> {
                  
                  viewer.setStyle({}, {cartoon:{style:"oval", ribbon:true}});
                  viewer.render();
              });
