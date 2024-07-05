
       $3Dmol.getbin('data/5ire.bcif.gz',function(data) {
           viewer.addModel(data,'bcif',{doAssembly: true, noComputeSecondaryStructure: true});
           viewer.setStyle({cartoon:{color:'spectrum'}});
       viewer.zoomTo();
       viewer.zoom(2);
       viewer.render();
    });
