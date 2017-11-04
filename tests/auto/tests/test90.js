
    $.get("data/Al.lammpstrj",
        function (data){
                var m = viewer.addModel(data, "lammpstrj");	
                viewer.setStyle({},{sphere:{radius:0.05}});
                viewer.zoomTo();
                viewer.zoom(4);
                viewer.animate({loop:"forward",reps:1});
                viewer.render();
        });

