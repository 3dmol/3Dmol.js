$.get('data/6g6j.pdb', function(pdb) {
      viewer.addModel(pdb,'pdb',{'hbondCutoff':2});
      viewer.setStyle({cartoon:{color:'spectrum',colorscheme:'roygb'}});
      viewer.zoomTo();
      viewer.render();
});
