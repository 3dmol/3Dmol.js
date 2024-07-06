$3Dmol.get("data/1wav.cif",function(data){
  viewer.addModel(data,'cif');

  viewer.setStyle({stick:{colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.ROYGB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
  viewer.render();
});
