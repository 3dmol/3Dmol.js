
              $3Dmol.download("pdb:1MO8",viewer,{multimodel:true, frames:true},function(){
                  
                  viewer.setStyle({}, {cartoon:{color:"spectrum"}});
                  viewer.animate({loop:"backward"});
                  viewer.render();  
                  
              });

              function toggleAnimation() {
                var button = document.getElementById("button");
                if (button.value == "Stop Animation") {
                    button.value = "Start Animation";
                    viewer.stopAnimate();
                } else {
                    button.value = "Stop Animation";
                    viewer.animate({loop:"backward"});
                }
           }