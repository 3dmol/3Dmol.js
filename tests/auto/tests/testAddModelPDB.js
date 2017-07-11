
    
    $.get('volData/4csv.pdb', function(data) {
      viewer.addModel(data,'pdb');
      viewer.setStyle({cartoon:{},stick:{}});
      viewer.zoomTo();
      viewer.render( /*no callback*/);
    
    //can't use jquery with binary data
    $3Dmol.getbin('volData/4csv.ccp4.gz', function(data) {     
       var voldata = new $3Dmol.VolumeData(data, 'ccp4.gz');
      viewer.addIsosurface(voldata, {isoval: 0.25,
                                       color: "blue",
                                       wireframe: true,
                                       selectedRegion: viewer.selectedAtoms({}),
                                       selectedOffset: 3,
                                       radius: 3.0                                   
                                    });
                                                     
      viewer.render();
    });
});
