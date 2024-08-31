
$3Dmol.download("pdb:4DLN",viewer,{},function(){
  viewer.setStyle({'cartoon':{colorscheme: 'ssPyMol'}});
  var sb = viewer.addSurface("SAS", {colorscheme: 'ssJmol'} ,{hetflag:false},null,null, function(sb) {
    viewer.render( /* no callback */);
    viewer.setSurfaceMaterialStyle(sb, {colorscheme: 'ssPyMol',opacity:.7});
    viewer.render();
  });  
});
