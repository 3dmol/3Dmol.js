$3Dmol.download("pdb:2ABJ",viewer,{},function(){
                  
                  viewer.setViewStyle({style:"ambientOcclusion",strength:1.2});
                  viewer.setStyle({chain:'A'},{sphere:{hidden:true}});
                  viewer.setStyle({chain:'D'},{sphere:{radius:3.0}});
                  viewer.setStyle({chain:'G'},{stick:{colorscheme:'greenCarbon',radius:1.5}});
                  viewer.setStyle({chain:'J'},{cartoon:{color:'blue'}});
        
                  viewer.render();
              });
