             
              $3Dmol.download("pdb:2V0E",viewer)
              .then((m1)=> {
                  $3Dmol.download("mmtf:4HHB",viewer)
                  .then((m2)=> {
                      m1.setStyle('sphere:colorscheme~greenCarbon');
                      m2.setStyle({}, 'stick');
                      viewer.render();
                  });
              });
