var File = new XMLHttpRequest();
File.open("GET", "volData/Al.lammpstrj", false);
File.onreadystatechange = function(){
    if(File.readyState === 4){
	if(File.status === 200 || File.status == 0){
	    var coordinatesData = File.responseText;
	    m.setCoordinates(coordinatesData, "lammpstrj");
	    viewer.setStyle({},{sphere:{}});
	    viewer.zoomTo();
	    viewer.render();
	}
    }
}
File.send(null);
