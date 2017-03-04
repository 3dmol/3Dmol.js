$3Dmol.download("pdb:4UB9",viewer,{},function(){
                  
        viewer.getModel(0).setClickable({}, true, function(){console.log(console.log(this))});

                  viewer.setStyle({chain:'A'},{cartoon:{color:"red",tubes:true}});
                  viewer.setStyle({chain:'B'},{cartoon:{color:"green",ribbon:true}});
                  viewer.setStyle({chain:'C'},{cartoon:{color:"orange"}});
                  viewer.setStyle({chain:'D'},{cartoon:{color:"yellow"}});
                  viewer.setStyle({chain:'E'},{cartoon:{color:"pink"}});
                  viewer.setStyle({chain:'F'},{cartoon:{color:"red"}});
                  viewer.setStyle({chain:'G'},{cartoon:{color:"red"}});
                  viewer.setStyle({chain:'H'},{cartoon:{color:"blue"}});
                  viewer.setStyle({chain:'A'},{cartoon:{color:"red",tubes:true}});
                  viewer.setStyle({chain:'B'},{cartoon:{color:"green",ribbon:true}});
                  viewer.setStyle({resn:"GLY"},{cartoon:{arrows:true,color:"blue"}})
  				  viewer.zoomTo({chain:'C'})
  	                viewer.render();
              });

