
    
    $.get('volData/4csv.pdb', function(data) {
      viewer.addModel(data,'pdb');
      viewer.setStyle({cartoon:{},stick:{}});
      viewer.zoomTo();
      viewer.render();
    });
    
    //can't use jquery with binary data
    var req = new XMLHttpRequest();
    req.open('GET', 'volData/4csv.ccp4.gz', true);
    req.responseType = "arraybuffer";
    req.onload = function (aEvt) {      
       var voldata = new $3Dmol.VolumeData(req.response, 'ccp4.gz');
                          
      //viewer.translate(10,10,1000);         
      viewer.zoomTo({resn:'STI'});
      //viewer.zoom(10,1000);
      //viewer.rotate(90,"y",1000);
      viewer.render();
    };