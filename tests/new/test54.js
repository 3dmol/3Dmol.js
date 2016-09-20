var viewer = $3Dmol.createViewer($("#gldiv"), {
		defaultcolors : $3Dmol.rasmolElementColors
		});
    viewer.setBackgroundColor(0xffffffff);
    var m;
    var rawFile = new XMLHttpRequest();
    rawFile.open("GET", "volData/2water.gro", false);
    rawFile.onreadystatechange = function (){
        if(rawFile.readyState === 4){
            if(rawFile.status === 200 || rawFile.status == 0){
                var data = rawFile.responseText;
		m = viewer.addModel(data, "gro");	
		viewer.setStyle({},{sphere:{}});
		viewer.zoomTo();
		viewer.render(callback);
	    }	
        }
    }
    rawFile.send(null);    