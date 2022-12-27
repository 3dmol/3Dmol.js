/* @script
   $3Dmol.initShapes = function(viewer) { 
    $3Dmol.get('../test_structs/benzene-homo.cube', function(data){
        var voldata = new $3Dmol.VolumeData(data, "cube");
        viewer.addIsosurface(voldata, {isoval: 0.01, color: "blue", alpha: 0.95, smoothness: 10});              
        viewer.addIsosurface(voldata, {isoval: -0.01, color: "red", alpha: 0.95, smoothness: 10}); 
        viewer.zoomTo();
        viewer.zoom(.75);
        viewer.render();
    },'text');
};
*/
  /* @div
  <div  class='viewer_3Dmoljs'  style="width: 400px; height: 400px;" data-href="../test_structs/benzene.sdf" data-style="stick" data-callback='$3Dmol.initShapes' data-backgroundcolor='0xf6f6f6'></div>
*/


