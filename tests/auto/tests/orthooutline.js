        

    viewer.setViewStyle({style:"outline"});
    viewer.setCameraParameters({orthographic: true});
    $.get('data/1fas.pqr', function(data){
        viewer.addModel(data, "pqr");
        viewer.setStyle({'cartoon':{},'stick':{},'sphere':{radius:0.5}});
        viewer.zoomTo({resi:1});
        viewer.render();
    });


