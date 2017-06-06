$3Dmol.download("pdb:5TZ8", viewer, {}, function(){
        			
                  viewer.setStyle({chain:'A'}, {cartoon:{arrows:true, opacity:0.8, color:'blue'}});
                  viewer.setStyle({chain:'B'}, {cartoon:{style:"trace",color:"green"}});
                  viewer.setStyle({chain:'C'}, {cartoon:{arrows:false, opacity:0.8, color:'red'}});

                  viewer.setStyle({chain: 'A',resi:'20-30'},{cartoon:{color:'black'}});
                  viewer.setStyle({chain:'A',resi:[1,3,5,7,9,11]},{cartoon:{color:'yellow'}});
                  viewer.render();
              });
