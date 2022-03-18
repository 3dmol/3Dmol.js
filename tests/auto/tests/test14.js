 $.get('../test_structs/benzene-homo.cube', (data)=> {
                  const voldata = new $3Dmol.VolumeData(data, "cube");
                  viewer.addIsosurface(voldata, {isoval: 0.01,
                                                 color: "blue",
                                                 alpha: 0.5,
                                                 smoothness: 10});
                  viewer.addIsosurface(voldata, {isoval: -0.01,
                                                 color: "red",
                                                 smoothness: 5,
                                                 opacity:0.5,
                                                 wireframe:true,
                                                 clickable:true,
                                                 callback() {
                                                     this.opacity = 0.0;
                                                     viewer.render( );
                                                 }});
                  viewer.setStyle({}, {stick:{}});
                  viewer.zoomTo();
                  viewer.render();
                });
