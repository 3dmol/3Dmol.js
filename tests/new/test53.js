
    var m;
    var rawFile = new XMLHttpRequest();
    rawFile.open("GET", "volData/TC5b.prmtop", false);
    rawFile.onreadystatechange = function (){
        if(rawFile.readyState === 4){
            if(rawFile.status === 200 || rawFile.status == 0){
                var data = rawFile.responseText;
		m = viewer.addModel(data, "prmtop");	
                var req = new XMLHttpRequest();
   		req.open("GET", "volData/heat1.mdcrd.gz", true);
		req.responseType = "arraybuffer";
    		req.onload = function(aEvt){
		    m.setCoordinates(req.response, "mdcrd.gz");
		    viewer.setStyle({},{sphere:{}});
		    viewer.zoomTo();
		    viewer.animate({loop:"forward"});
		    viewer.render();
		};
		req.send(null);
	    }	
        }
    }
    rawFile.send(null);    