             
              const m1 = $3Dmol.download("pdb:2V0E",viewer,{multimodel:true, frames:true},()=> {
                const m2 = $3Dmol.download("pdb:5ZB6",viewer,{multimodel:true, frames:true},()=> {
                  
                  m1.setStyle({}, {stick:{color:"red"}});
                  m2.setStyle({}, {cartoon:{color:"blue"}});
                  viewer.animate({loop:"forward",reps:1});
                  viewer.render();
                });
              });
