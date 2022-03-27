viewer.setDefaultCartoonQuality(20);
$3Dmol.download("pdb:2BTF",viewer,{},function(){
                  
                  viewer.setStyle({chain:'A'},{cartoon:{color:'blue',style:'parabola',arrows:true}});
                  viewer.setStyle({chain:'P'},{cartoon:{color:'green',style:'oval',arrows:true}});
         
                  viewer.render();
              });
