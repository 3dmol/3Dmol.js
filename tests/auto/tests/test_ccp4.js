
    
$.get('data/1lo6.pdb', function(data) {
      viewer.addModel(data,'pdb');
      viewer.setStyle({cartoon:{},stick:{}});

    //can't use jquery with binary data
    $3Dmol.getbin('data/1lo6_2FOFC.ccp4', function(data) {     
       var voldata = new $3Dmol.VolumeData(data, 'ccp4');
      viewer.addIsosurface(voldata, {isoval: 0.1,
                                       color: "green",
                                       opacity: .5,
                                       coords: viewer.selectedAtoms({resi: 18}),
                                       seldist: 3.0                                   
                                    });
      viewer.addIsosurface(voldata, {isoval: 0.3,
                                       color: "blue",
                                       wireframe: true,
                                       coords: viewer.selectedAtoms({resi: 18}),
                                       seldist: 3.0                                   
                                    });                                                     
      viewer.zoomTo({resi: 18});
      viewer.render();
    });
});
