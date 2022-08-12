
// Apparently aromatic bonds can have a bond order of 4, which we weren't handling
              $.get('../test_structs/aromaticsdf.sdf', function(data){
                  viewer.addModel(data, "sdf");
                  viewer.setStyle({stick:{}});
                  viewer.zoomTo();
                  viewer.render();
              });
