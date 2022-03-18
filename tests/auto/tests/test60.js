

    $3Dmol.download('pdb:1pfl',viewer,{},(m) => {
        m.setStyle({'cartoon':{colorscheme:'ssPyMol'}});
       viewer.zoomTo();
       viewer.render();
    });
