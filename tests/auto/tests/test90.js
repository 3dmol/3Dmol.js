
    $.get("volData/Al.lammpstrj",
    	function (data){
		m = viewer.addModel(data, "lammpstrj");	
	        viewer.setStyle({},{sphere:{radius:0.05}});
	        viewer.zoomTo();
                viewer.zoom(4);
                viewer.animate({loop:"forward"});
	        viewer.render();
	    });
