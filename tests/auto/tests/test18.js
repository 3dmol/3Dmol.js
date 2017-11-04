$3Dmol.download("pdb:4UND",viewer,{},function(){
                  

                  viewer.setStyle({chain:'A',within:{distance:10.0,sel:{chain:'B'}}},{sphere:{color:'blue',radius:1.0}});
                  viewer.setStyle({chain:'B'},{stick:{color:'red'}});

                  viewer.render();
              });
