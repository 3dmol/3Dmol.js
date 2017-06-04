var oldval = $3Dmol.syncSurface;

var setStyles = function(volumedata){
    $3Dmol.syncSurface = true;    
    var data = new $3Dmol.VolumeData(volumedata, "cube.gz");
    viewer.addSurface("VDW", {color:'red'} ,{chain:'A'},null,null, function(sa) {
        viewer.setSurfaceMaterialStyle(sa, {color: 0x00ff00, opacity: 0.5});
      
    });
    var sb = viewer.addSurface("SAS", {color:'red'} ,{chain:'B'},null,null, function(sb) {
        viewer.setSurfaceMaterialStyle(sb, {opacity:1.0, voldata: data, volscheme: new $3Dmol.Gradient.RWB(-10,10)});
      
    });
    $3Dmol.syncSurface = oldval; //otherwise we change this for all other tests
    viewer.render();
};
$3Dmol.download("pdb:4DLN",viewer,{},function(){
  $3Dmol.syncSurface = true;    
  $3Dmol.getbin("volData/4dln.cube.gz",setStyles);
  viewer.render( /*no callback*/);
});
