$3Dmol.download("pdb:4DM7",viewer,{},function(){
                  

                  viewer.setStyle({chain:'A',resn:'GLU'},{sphere:{color:'cyan'}});
                  viewer.setStyle({chain:'B',atom:'CG'},{sphere:{color:'teal'}});
                  viewer.setStyle({chain:'C',elem:'O'},{sphere:{color:'maroon'}});
                  viewer.setStyle({chain:'D',bonds:0},{sphere:{color:'gray'}});
                  viewer.render();
              });
