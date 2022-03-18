
    $.get("https://3dmol.csb.pitt.edu/mdsrv/file/data/md.gro",  (data)=> {
        const m = viewer.addModel(data, "gro");
        m.setStyle({stick:{}});
        viewer.zoomTo();
        viewer.render( /* no callback */ );
        const url = "https://3dmol.csb.pitt.edu/mdsrv/";
        const pathToFile = "data/md.xtc";
        
        m.setCoordinatesFromURL(url, pathToFile)
        .then(() => {
            viewer.setStyle({line:{},cartoon:{}});
            viewer.zoomTo();
            viewer.animate({loop:"forward",reps:1,step:10});
        }).catch().then(() => {viewer.render()});
        
    });
