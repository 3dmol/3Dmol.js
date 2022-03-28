
    $.get("https://3dmol.csb.pitt.edu/mdsrv/file/data/md.gro",  function (data){
        var m = viewer.addModel(data, "gro");
        m.setStyle({stick:{}});
        viewer.zoomTo();
        viewer.render( /*no callback*/ );
        var url = "https://3dmol.csb.pitt.edu/mdsrv/";
        var pathToFile = "data/md.xtc";
        
        m.setCoordinatesFromURL(url, pathToFile)
        .then(function() {
            viewer.setStyle({line:{},cartoon:{}});
            viewer.zoomTo();
            viewer.animate({loop:"forward",reps:1,step:10});
        }).catch().then(function() {viewer.render()});
        
    });
