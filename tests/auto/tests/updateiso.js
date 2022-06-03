 $.get('../test_structs/benzene-homo.cube', function(data){
                  var voldata = new $3Dmol.VolumeData(data, "cube");
                  var blueiso = viewer.addIsosurface(voldata, {isoval: 0.01,
                                                 color: "blue",
                                                 alpha: 0.5,
                                                 smoothness: 10});
                  var rediso = viewer.addIsosurface(voldata, {isoval: -0.01,
                                                 color: "red",
                                                 smoothness: 5,
                                                 opacity:0.5,
                                                 wireframe:true,
                                                 clickable:true,
                                                 callback:
                                                 function() {
                                                     this.opacity = 0.0;
                                                     viewer.render( );
                                                 }});
                  viewer.setStyle({}, {stick:{}});
                  viewer.zoomTo();
                  rediso.updateStyle({color:"purple"});
                  blueiso.updateStyle({color:"cyan",alpha:0});
                  viewer.render();
                });
