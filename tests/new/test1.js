var viewer = $3Dmol.createViewer($("#test1"));
              $3Dmol.download("pdb:4UND",viewer,{},function(){
                  viewer.setBackgroundColor(0xffffffff);
                  viewer.render(callback);
              });