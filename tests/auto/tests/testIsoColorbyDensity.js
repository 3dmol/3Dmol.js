$.get('data/chem/density.cub', (density) => {
  viewer.addModel(density,'cube');
  viewer.setStyle({stick:{}});
  viewer.zoomTo();
  viewer.render( /* no callback */);

 $.get('data/chem/esp.cub', (esp) => {
     const densvol = new $3Dmol.VolumeData(density, 'cube');
     const espvol = new $3Dmol.VolumeData(esp, 'cube');
     viewer.addIsosurface(densvol, {isoval: 0.001, opacity: 0.9,
            smoothness: 4, voldata: espvol,
            volscheme: {gradient:'rwb', min:-.025, max:.025}});
     viewer.render();
 });
});

