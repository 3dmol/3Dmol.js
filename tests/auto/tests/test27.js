

              $.get('data/1fas.pqr', function(data){
                  viewer.addModel(data, "pqr");
                  $.get("data/1fas.cube",function(volumedata){
                      var range = $3Dmol.getPropertyRange(viewer.selectedAtoms(),'partialCharge');
                      viewer.addSurface($3Dmol.SurfaceType.VDW, {
                          opacity:0.85,
                          voldata: new $3Dmol.VolumeData(volumedata, "cube"),
                          volscheme: new $3Dmol.Gradient.ROYGB(range[1],range[0]) //reverse color order
                      },{});
                      
                  viewer.render();
                  });
                  viewer.zoomTo();
              });
