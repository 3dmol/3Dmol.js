 $3Dmol.getbin('data/Dt.cube.gz', function(density) {
   $3Dmol.getbin('data/ESP.cube.gz', function(esp) {
     viewer.addModel(density, "cube.gz");
     viewer.setStyle("stick");
     viewer.addVolumetricData(density, "cube.gz", {isoval: 0.005, smoothness: 2, opacity:.9, 
            voldata: esp, volformat: 'cube.gz', 
           volscheme: {gradient:'rwb', min:-.1, max:.1}});
     viewer.zoomTo();
     viewer.render();
   });
 });
