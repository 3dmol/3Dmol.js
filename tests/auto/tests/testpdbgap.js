
// there are gaps in the backbone of this pdb
              $.get('../test_structs/jk.pdb', (data)=> {
                  viewer.addModel(data, "pdb");
                  viewer.setStyle({stick:{}});
                  viewer.zoomTo();
                  viewer.render();
              });
