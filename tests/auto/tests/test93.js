
    $.get("volData/model1.prmtop",  function (data){
        var m = viewer.addModel(data, "prmtop");
        var url = "127.0.0.1:800/mdsrv";
        var pathToFile = "root/3Dmol.js/tests/auto/volData/model1_md2.nc";
        $3Dmol.runGenerator(m.setCoordinatesFromURL(url,pathToFile), function(){
            viewer.setStyle({},{sphere:{}});
            viewer.zoomTo();
            viewer.animate({loop:"forward",reps:1});
        });
    });
