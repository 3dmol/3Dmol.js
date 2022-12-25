
viewer.addBox({center:{x:0,y:0,z:0},dimensions: {w:3,h:3,d:2},color:'darkblue'});
viewer.zoomTo();
viewer.zoom(0.5);
viewer.addSphere({center:{x:10,y:0,z:0},radius:3,color:'grey'})
viewer.render( /*no callback */ );
viewer.center(undefined, 1000, true);
viewer.rotate(90);
viewer.render();
