
       $3Dmol.download('pdb:2nbd',viewer,{onemol: true,multimodel: true},(m) => {
        m.setStyle({'cartoon':{colorscheme:{prop:'ss',map:$3Dmol.ssColors.Jmol}}});
       viewer.zoomTo();
       viewer.render();
    });
