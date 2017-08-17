        var setStyle = function(style) {
                glviewer1.setStyle(style); 
                glviewer2.setStyle(style); 

        };

    var glviewer1;
    var glviewer2;

    $.get("volData/1fas.pqr",
        function(data) {
            
            var m = viewer.addModel(data, "1fas.pqr");

            [glviewer1, glviewer2] = $3Dmol.createStereoViewer("gldiv", m.selectedAtoms({}));
            glviewer1.addModel(data, "1fas.pqr");
            glviewer2.addModel(data,"1fas.pqr");
            setStyle({cartoon:{color:'spectrum'}})

            glviewer1.zoomTo();
            glviewer2.zoomTo();

            glviewer1.render();
            glviewer2.render();

        });
