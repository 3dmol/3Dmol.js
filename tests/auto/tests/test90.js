    var rawFile = new XMLHttpRequest();
    rawFile.open("GET", "volData/Al.lammpstrj", false);
    rawFile.onreadystatechange = function (){
        if(rawFile.readyState === 4){
            if(rawFile.status === 200 || rawFile.status == 0){
                var data = rawFile.responseText;
		var m = viewer.addModel(data, "lammpstrj");	
	        viewer.setStyle({},{sphere:{radius:0.05}});
	        viewer.zoomTo();
                viewer.zoom(4)
                viewer.animate({loop:"forward"});
	        viewer.render();
	    }	
        }
    }
    rawFile.send(null);
