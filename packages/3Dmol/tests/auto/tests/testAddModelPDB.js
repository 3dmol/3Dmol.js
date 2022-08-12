
    
$.get('data/4csv.pdb', function(data) {
  viewer.addModel(data,'pdb');
  viewer.setStyle({cartoon:{},stick:{}});
  viewer.zoomTo();
  viewer.render();
  
});
