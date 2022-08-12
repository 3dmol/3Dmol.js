
var setStyles = function(volumedata){
    var data = new $3Dmol.VolumeData(volumedata, "cube.gz");
    viewer.addSurface("VDW", {color:'red'} ,{chain:'A'},null,null, function(sa) {
        viewer.render( /* no callback */);
        viewer.setSurfaceMaterialStyle(sa, {color: 0x00ff00, opacity: 0.5});
        var sb = viewer.addSurface("SAS", {color:'red'} ,{chain:'B'},null,null, function(sb) {
          viewer.render( /* no callback */);
          viewer.setSurfaceMaterialStyle(sb, {opacity:1.0, voldata: data, volscheme: new $3Dmol.Gradient.RWB(-10,10)});
          viewer.render();
        });      
    });


};
$3Dmol.download("pdb:4DLN",viewer,{},function(){
  $3Dmol.getbin("data/4dln.cube.gz",setStyles);
  viewer.render( /*no callback*/);
});
