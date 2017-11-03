
    $.get("http://proteinformatics.charite.de/tool-mdsrv/file/MDsrv/3PQR/GaCT/md02/md02.gro",  function (data){
        var m = viewer.addModel(data, "gro");
        m.setStyle({stick:{}});
        viewer.zoomTo();
        viewer.render( /*no callback*/ );
        var url = "http://proteinformatics.charite.de/tool-mdsrv/";
        var pathToFile = "MDsrv/3PQR/GaCT/md02/md02.xtc";
        
        m.setCoordinatesFromURL(url, pathToFile)
        .then(function() {
            viewer.setStyle({line:{},cartoon:{}});
            viewer.zoomTo();
            viewer.animate({loop:"forward",reps:1,step:1000});
        }).catch().then(function() {viewer.render()});
        
    });
