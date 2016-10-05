var viewer = $3Dmol.createViewer($("#gldiv"), {
		defaultcolors : $3Dmol.rasmolElementColors
		});
    
    var rawFile = new XMLHttpRequest();
    rawFile.open("GET", "volData/TRPcage.prmtop", false);
    rawFile.onreadystatechange = function (){
        if(rawFile.readyState === 4){
            if(rawFile.status === 200 || rawFile.status == 0){
                var data = rawFile.responseText;
		var m = viewer.addModel(data, "prmtop");	
                var File = new XMLHttpRequest();
   		File.open("GET", "volData/TRPcage.inpcrd", false);
    		File.onreadystatechange = function(){
		    if(File.readyState === 4){
	    		if(File.status === 200 || File.status == 0){
			    var coordinatesData = File.responseText;
			    m.setCoordinates(coordinatesData, "inpcrd");
			    viewer.setStyle({},{sphere:{}});
			    viewer.zoomTo();
			    viewer.render(callback);
			}
		    }
		}
		File.send(null);
	    }	
        }
    }
    rawFile.send(null);    