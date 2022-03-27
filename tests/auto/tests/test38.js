
              $3Dmol.download("pdb:1MO8",viewer,{multimodel:true, frames:true},function(){
                  
                  viewer.setStyle({}, {cartoon:{color:"spectrum"}});
                  viewer.animate({loop:"backward",reps:1});
                  viewer.render();  
                  
              });

