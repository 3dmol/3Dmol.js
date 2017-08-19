

    var stereoViewer = new $3Dmol.createStereoViewer("gldiv", 10);

    $.get("volData/TC5b.prmtop",
        function(data) {

            var models = stereoViewer.addModel(data, "prmtop");
            $.get("volData/heat1.mdcrd", 
            function(data) {
                stereoViewer.setCoordinates(models, data, "mdcrd");
                stereoViewer.setStyle({},{sphere:{}});
                stereoViewer.zoomTo();
                stereoViewer.animate({loop:"forward",reps:1});
                stereoViewer.render(callback);
            });
        });
