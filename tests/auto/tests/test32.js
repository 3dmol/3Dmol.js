$3Dmol.download("pdb:1YCR", viewer, {}, ()=> {
                  
                  viewer.setStyle({}, {cartoon:{tubes:true, style:"rectangle"}});
                  viewer.render();
              });
