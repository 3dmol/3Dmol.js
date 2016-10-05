
              $.get('test_structs/multiple.pqr', function(data){
                  viewer.addModels(data, "pqr");
                  viewer.zoomTo();
                  viewer.render();
              });