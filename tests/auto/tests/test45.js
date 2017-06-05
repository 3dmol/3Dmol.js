
    $.get('volData/4csv.pdb', function(data) {
      viewer.addModel(data,'pdb');
      viewer.setStyle({cartoon:{},stick:{}});
      viewer.zoomTo();
      $3Dmol.getbin('volData/4csv.ccp4.gz', function(data) {
        var voldata = new $3Dmol.VolumeData(data, 'ccp4.gz');
        viewer.addIsosurface(voldata, {isoval: 0.25,
                                       color: "blue",
                                       wireframe: true,
                                       linewidth:0.005,
                                       selectedRegion: viewer.selectedAtoms({}),
                                       selectedOffset: 3,
                                       radius: 3.0                                   
                                    });
        viewer.render();
      });
    });
    
 
