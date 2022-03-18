

    $3Dmol.download('pdb:1pfl',viewer,{},(m) => {
        m.setStyle({'cartoon':{colorscheme:{prop:'ss',map:$3Dmol.ssColors.Jmol}}});
       viewer.zoomTo();
       viewer.render();
    });
