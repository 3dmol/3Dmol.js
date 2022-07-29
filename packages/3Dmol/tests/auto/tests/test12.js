$3Dmol.download("pdb:3VOV",viewer,{},function(){
                  
                  viewer.setStyle({chain:'A'},{cartoon:{color:'spectrum'}});
                  viewer.setStyle({chain:'B'},{cartoon:{style:'trace'}});
                  viewer.setStyle({chain:'C'},{cartoon:{color:'red'}});
                  viewer.setStyle({chain:'D'},{cartoon:{thickness:2.0}});
                  viewer.render();
              });
