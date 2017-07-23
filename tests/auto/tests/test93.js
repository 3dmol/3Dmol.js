
    $.get("volData/md.xyz",  function (data){
        var m = viewer.addModel(data, "xyz");
        var url = "proteinformatics.charite.de/mdsrvdev";
        var pathToFile = "data2/md.xtc";
        m.setCoordinatesFromURL(url, pathToFile, function(){
            viewer.setStyle({},{sphere:{}});
            viewer.zoomTo();
            viewer.animate({loop:"forward",reps:1});
            viewer.render();
        });
    });
