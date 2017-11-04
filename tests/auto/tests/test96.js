
/*
    var stereoViewer = new $3Dmol.createStereoViewer("gldiv", 10);

    $.get("data/TC5b.prmtop",
        function(data) {
            var models = stereoViewer.addModel(data, "prmtop");
            $.get("data/heat1.mdcrd", 
            function(data) {
                stereoViewer.setCoordinates(models, data, "mdcrd");
                stereoViewer.setStyle({},{sphere:{}});
                stereoViewer.zoomTo();
                var dist = stereoViewer.getPerceivedDistance();
                stereoViewer.setPerceivedDistance(dist[0]-10); //useful for zooming in or out
                stereoViewer.setAutoEyeSeparation(); //optional. changes the camera x value
                stereoViewer.animate({loop:"forward",reps:1});
                stereoViewer.render(callback);
            });
        });
*/
