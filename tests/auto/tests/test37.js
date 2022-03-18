

              $.get('../test_structs/multiple.sdf', (data)=> {
                  viewer.addAsOneMolecule(data, "sdf");
                  viewer.zoomTo();
                  viewer.render();
              });
