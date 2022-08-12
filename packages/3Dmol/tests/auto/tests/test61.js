
    var m;
    var rawFile = new XMLHttpRequest();
    rawFile.open("GET", "data/4csv.pdb", false);
    rawFile.onreadystatechange = function (){
        if(rawFile.readyState === 4){
            if(rawFile.status === 200 || rawFile.status == 0){
                var data = rawFile.responseText;
		        m = viewer.addModel(data, "pdb");	
                var req = new XMLHttpRequest();
   		req.open("GET", "data/4csv.pdb", false);
    		req.onload = function(aEvt){
		    m.setCoordinates(req.response, "pdb");
		    viewer.zoomTo();
		    viewer.render();
		};
		req.send(null);
	    }	
        }
    }
    rawFile.send(null);    
