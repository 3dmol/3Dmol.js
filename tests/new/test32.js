$3Dmol.download("pdb:1YCR", viewer, {}, function(){
                  viewer.setBackgroundColor(0xffffffff);
                  viewer.setStyle({}, {cartoon:{tubes:true, style:"rectangle"}});
                  viewer.render();
              });