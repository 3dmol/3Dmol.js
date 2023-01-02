/* @script
   $3Dmol.testclicky = function(viewer) { 
    viewer.setHeight(400);
    viewer.setWidth(400);
    viewer.zoomTo();
    viewer.setClickable({},true,function(atom) {
      viewer.removeAllShapes();
      viewer.addSphere({center:atom,radius:1.0,color:'purple',alpha:0.4});
      viewer.addLabel("Label",{position: atom, bold: true});
    viewer.render();
});
    viewer.render( );
   
    viewer._handleMouseDown({pageX: 127, pageY: 129, preventDefault: function(){}});
    viewer._handleMouseUp({pageX: 127, pageY: 129, preventDefault: function(){}});

  };
*/
  /* @div
  <div  class='viewer_3Dmoljs'  style="width: 400px; height: 400px;" data-backgroundColor="white" data-href="../test_structs/benzene.sdf" data-style="stick" data-callback='$3Dmol.testclicky'></div>
*/


