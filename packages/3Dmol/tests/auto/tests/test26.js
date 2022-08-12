

              $.get('data/1fas.pqr', function(data){
                  viewer.addModel(data, "pqr");
                  $.get("data/1fas.cube",function(volumedata){
                      viewer.addSurface($3Dmol.SurfaceType.VDW, {
                          opacity:0.85,
                          voldata: new $3Dmol.VolumeData(volumedata, "cube"),
                          volscheme: new $3Dmol.Gradient.RWB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'partialCharge'))
                      },{});
                      
                  viewer.render();
                  });
                  viewer.zoomTo();
              });
