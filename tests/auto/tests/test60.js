

    $3Dmol.download('pdb:1pfl',viewer,{},function(m) {
        m.setStyle({'cartoon':{colorscheme:'ssPyMol'}});
       viewer.zoomTo();
       viewer.render();
    });
