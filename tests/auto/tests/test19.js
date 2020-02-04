var setStyles = function(volumedata){
	var data = new $3Dmol.VolumeData(volumedata, "cube");
	viewer.zoomTo({chain:'A'});
        viewer.rotate(-110);
	viewer.render( /*no callback, why? */);
        viewer.addSurface("VDW", {opacity:0.85, voldata: data, volscheme: new $3Dmol.Gradient.RWB(-10,10)},{chain:'A'});
	//make sure we see the volume mapped color
        viewer.render( /*no callback */);
};
$3Dmol.download("pdb:4DLN",viewer,{},function(){
	$.get("data/1fas.cube",setStyles);
});
