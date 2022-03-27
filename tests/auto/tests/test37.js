

              $.get('../test_structs/multiple.sdf', function(data){
                  viewer.addAsOneMolecule(data, "sdf");
                  viewer.zoomTo();
                  viewer.render();
              });
