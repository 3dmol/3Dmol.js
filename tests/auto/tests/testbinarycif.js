
    $3Dmol.getbin('data/4hhb.bcif.gz', function(data) {
      viewer.addModel(data,'bcif.gz');
      viewer.setStyle({cartoon:{},stick:{}});
      viewer.zoomTo();
      viewer.render();
    });
    
 
