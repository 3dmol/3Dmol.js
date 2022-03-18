

    const viewer = new $3Dmol.createStereoViewer("gldiv");

    $.get("data/TC5b.prmtop",
        (data) => {
            const models = viewer.addModel(data, "prmtop");
            $.get("data/heat1.mdcrd", 
            (data) => {
                viewer.setCoordinates(models, data, "mdcrd");
                viewer.setStyle({},{sphere:{}});
                viewer.zoomTo();
                const dist = viewer.getPerceivedDistance();
                viewer.setPerceivedDistance(dist[0]-10); // useful for zooming in or out
                viewer.setAutoEyeSeparation(); // optional. changes the camera x value
                viewer.render();
            });
        });

