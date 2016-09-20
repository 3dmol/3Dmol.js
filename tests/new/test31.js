$3Dmol.download("pdb:1BNA", viewer, {}, function(){
                  viewer.setBackgroundColor(0xffffffff);
                  viewer.setStyle({}, {cartoon:{style:"oval", ribbon:true}});
                  viewer.render();
              });