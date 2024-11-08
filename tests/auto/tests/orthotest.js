viewer.setCameraParameters({orthographic:true,fov: 60});
viewer.addModel("2\n\nO  0 0 0\nC  1.5 0 0\n","xyz");
viewer.setStyle({stick:{radius:0.2}});
viewer.zoomTo({serial:1}); 
viewer.zoom(8);   
//viewer.rotate(60,"y");
viewer.render();
