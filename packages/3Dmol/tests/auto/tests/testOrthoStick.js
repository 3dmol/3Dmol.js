

viewer.setCameraParameters({orthographic: true});
viewer.addModel("2\n\nO  0 0 0\nC  1.5 0 0\n","xyz");
viewer.setStyle({stick:{radius:0.2}});
viewer.zoomTo(); 
viewer.zoom(8);   
viewer.rotate(90);
viewer.rotate(10,"z");
viewer.render();

