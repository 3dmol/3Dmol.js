
    $.get('data/4csv.pdb', function(data) {
      viewer.addModel(data,'pdb');
      viewer.setStyle({cartoon:{}});
      viewer.zoomTo();
      //viewer.render( /*no callback*/);
      $3Dmol.getbin('data/4csv.ccp4.gz', function(data) {
        viewer.addSurface("SAS", {opacity:0.9, voldata: data, volformat: 'ccp4.gz', 
          volscheme: {gradient:'rwb', min:-.25, max:.25}});
        viewer.render();
      });
    });
    
 
