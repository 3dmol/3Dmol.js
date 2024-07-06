$3Dmol.download("bcif:1wav",viewer, {}, function(){

  viewer.setStyle({stick:{colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.ROYGB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
  viewer.render();
});
