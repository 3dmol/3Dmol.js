$3Dmol.download("pdb:4YGY", viewer, {}, function(){
                  
                  viewer.setStyle({chain:'A'}, {cartoon:{arrows:true, opacity:0.8, color:'spectrum'}});
                  viewer.setStyle({chain:'B'}, {cartoon:{style:"trace"}});
                  viewer.render();
              });
