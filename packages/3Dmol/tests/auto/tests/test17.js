$3Dmol.download("pdb:4UND",viewer,{},function(){
                  

                  viewer.setStyle({chain:'A',invert:true},{sphere:{color:'blue',radius:1.0}});
                  viewer.setStyle({chain:'A',resi:'669',expand:5.0},{sphere:{color:0xC0C0A902}});
                  viewer.setStyle({chain:'A',resi:860,expand:5.0,byres:true},{sphere:{color:'purple'}});
                  viewer.render();
              });
