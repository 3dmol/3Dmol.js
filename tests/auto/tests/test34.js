$3Dmol.download("pdb:1DC9", viewer, {}, ()=> {
                  
                  viewer.setViewStyle({style:"outline", color:"black", width:0.1});
                  viewer.setStyle({}, {cartoon:{arrows:true, tubes:true, style:"oval", color:'white'}});
                  viewer.render();
              });
