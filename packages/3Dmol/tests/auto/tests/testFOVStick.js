

viewer.setCameraParameters({orthographic:false,fov: 60});
viewer.addModel("2\n\nO  0 0 0\nC  1.5 0 0\n","xyz");
viewer.setStyle({stick:{radius:0.2}});
viewer.zoomTo(); 
viewer.zoom(4);   
viewer.rotate(10,"z");
viewer.render();

