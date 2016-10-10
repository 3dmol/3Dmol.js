
    
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
      viewer.addIsosurface(voldata, {isoval: 0.25,
                                       color: "blue",
                                       wireframe: true,
                                       linewidth:0.005,
                                       selectedRegion: viewer.selectedAtoms({}),
                                       selectedOffset: 3,
                                       radius: 3.0                                   
                                    });

                                                     
      viewer.render();
    };
    req.send(null);


