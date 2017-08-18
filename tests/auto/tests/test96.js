

        var setStyle = function(style) {
                glviewer1.setStyle(style); 
                glviewer2.setStyle(style); 

        };

    var glviewer1;
    var glviewer2;
    [glviewer1, glviewer2] = $3Dmol.createStereoViewer("gldiv", 10);

    $.get("volData/TC5b.prmtop",
        function(data) {

            var m1 = glviewer1.addModel(data, "prmtop");
            var m2 = glviewer2.addModel(data,"prmtop");

            $.get("volData/heat1.mdcrd", 
            function(data) {
                m1.setCoordinates(data, "mdcrd");
                m2.setCoordinates(data, "mdcrd");
                glviewer1.setStyle({},{sphere:{}});
                glviewer2.setStyle({},{sphere:{}});
                glviewer1.zoomTo();
                glviewer2.zoomTo();
                glviewer1.animate({loop:"forward",reps:1});
                glviewer2.animate({loop:"forward",reps:1});
                glviewer1.render(callback);
                glviewer2.render(callback);
            });
        });
