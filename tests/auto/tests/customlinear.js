$3Dmol.get('../test_structs/af.pdb', function(data){
    viewer.addModel(data);
    viewer.setStyle({cartoon:{colorscheme:{prop: 'b', 
      gradient:'linear_red_orange_yellow_white', min: 70, max: 100}}});
    viewer.zoomTo();
    viewer.render();
  });
