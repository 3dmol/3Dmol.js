/* @script
   $3Dmol.testmove = function(viewer) { 

    viewer.setClickable({},true,function(atom) {
      console.log('clicked');
      viewer.removeAllShapes();
      viewer.addSphere({center:atom,radius:1.0,color:'purple',alpha:0.4});
      viewer.addLabel("Label",{position: atom, bold: true});
    viewer.render( );
});
    viewer.render( );

   
    viewer._handleMouseDown({pageX: 127, pageY: 129, preventDefault: function(){}});
    viewer._handleMouseMove({pageX: 187, pageY: 129, preventDefault: function(){}});
    viewer._handleMouseUp({pageX: 187, pageY: 129, preventDefault: function(){}});

    viewer._handleMouseDown({pageX: 127, pageY: 129, preventDefault: function(){}});
    viewer._handleMouseMove({pageX: 127, pageY: 90, shiftKey: true, preventDefault: function(){}});
    viewer._handleMouseUp({pageX: 127, pageY: 139, preventDefault: function(){}});    


    viewer._handleMouseDown({pageX: 127, pageY: 129, preventDefault: function(){}});
    viewer._handleMouseMove({pageX: 247, pageY: 139, ctrlKey: true, preventDefault: function(){}});
    viewer._handleMouseUp({pageX: 247, pageY: 139, preventDefault: function(){}});        
  };
*/
  /* @div
  <div  class='viewer_3Dmoljs'  style="width: 400px; height: 400px;" data-backgroundColor="white" data-href="../test_structs/benzene.sdf" data-style="stick" data-callback='$3Dmol.testmove'></div>
*/


