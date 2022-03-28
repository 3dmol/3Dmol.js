$3Dmol.download("pdb:2ABJ",viewer,{},function(){
                  
                  viewer.setViewStyle({style:"outline"});
                  viewer.setStyle({chain:'A'},{sphere:{hidden:true}});
                  viewer.setStyle({chain:'D'},{sphere:{radius:3.0}});
                  viewer.setStyle({chain:'G'},{sphere:{colorscheme:'greenCarbon'}});
                  viewer.setStyle({chain:'J'},{sphere:{color:'blue'}});
                  viewer.addLine({dashed:true,start:{x:0,y:0,z:0},end:{x:100,y:100,z:100}});

         
                  viewer.render();
              });
