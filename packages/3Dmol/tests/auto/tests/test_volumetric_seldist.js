
    
    $.get('data/4csv.pdb', function(data) {
      viewer.addModel(data,'pdb');
      viewer.setStyle({cartoon:{},stick:{}});
      viewer.zoomTo();
      viewer.render( /*no callback*/);
    
    //can't use jquery with binary data
    $3Dmol.getbin('data/4csv.ccp4.gz', function(data) {     
       var voldata = new $3Dmol.VolumeData(data, 'ccp4.gz');
      viewer.addVolumetricRender(voldata, {        transferfn:[
            { color: "#ffffff", opacity: 0.0, value:  0 }, 
            { color: "#0000ff", opacity: .1, value: .5 } 
            ],
           coords: viewer.selectedAtoms({}),
           seldist: 3.0                                   
        });
                                                     
      viewer.render();
    });
});
