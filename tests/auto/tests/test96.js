        var setStyle = function(style) {
                glviewer1.setStyle(style); 
                glviewer2.setStyle(style); 
                glviewer1.render();
                glviewer2.render();
        };

    var viewers = $3Dmol.createStereoViewer("gldiv");
    var glviewer1 = viewers[0];
    var glviewer2 = viewers[1];

    $.get("volData/1fas.pqr",
        function(data) {

            glviewer1.addModel(data, "pqr");
            glviewer2.addModel(data, "pqr");            
           
            setStyle({cartoon:{color: 'spectrum'}})

            glviewer1.zoomTo();
            glviewer2.zoomTo();
        });
