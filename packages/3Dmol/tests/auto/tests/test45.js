
    $.get('data/4csv.pdb', function(data) {
      viewer.addModel(data,'pdb');
      viewer.setStyle({cartoon:{},stick:{}});
      viewer.zoomTo();
      //viewer.render( /*no callback*/);
      $3Dmol.getbin('data/4csv.ccp4.gz', function(data) {
        var voldata = new $3Dmol.VolumeData(data, 'ccp4.gz');
        viewer.addSurface("SES", {opacity:0.9, voldata: voldata, volscheme: new $3Dmol.Gradient.RWB(-.25,.25)});
        viewer.render();
      });
    });
    
 
