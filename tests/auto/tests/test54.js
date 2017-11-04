
    
    var m;
    var rawFile = new XMLHttpRequest();
    rawFile.open("GET", "data/2water.gro", false);
    rawFile.onreadystatechange = function (){
        if(rawFile.readyState === 4){
            if(rawFile.status === 200 || rawFile.status == 0){
                var data = rawFile.responseText;
		m = viewer.addModel(data, "gro");	
		viewer.setStyle({},{sphere:{}});
		viewer.zoomTo();
		viewer.render();
	    }	
        }
    }
    rawFile.send(null);    
