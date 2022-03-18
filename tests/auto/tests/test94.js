             
              $3Dmol.download("pdb:2V0E",viewer,{multimodel:true, frames:true})
              .then((m1)=> {
                  $3Dmol.download("mmtf:4HHB",viewer,{multimodel:true, frames:true})
                  .then((m2)=> {
                      m1.setStyle({}, {stick:{color:"red"}});
                      m2.setStyle({}, {cartoon:{color:"blue"}});
                      viewer.animate({loop:"forward",reps:1});
                      viewer.render();
                  });
              });
