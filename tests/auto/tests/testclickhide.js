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
              let cyl = viewer.addCylinder({start:{x:0.0,y:2.0,z:0.0},
                                  end:{x:0.0,y:10.0,z:-15.0},
                                  radius:3.5,
                                  fromCap:false,
                                  toCap:true,
                                  color:'teal'});
              let cyl2 = viewer.addCylinder({start:{x:0.0,y:0.0,z:0.0},
                                  end:{x:10.0,y:0.0,z:0.0},
                                  radius:1.0,
                                  fromCap:1,
                                  toCap:2,
                                  color:'red',
                                  hoverable:true,
                                  clickable:true,
                                  callback:function(){ cyl.color.setHex(0x00FFFF00);viewer.render();}
                                  }
                                 );

              viewer.addCylinder({start:{x:15.0,y:0.0,z:0.0},
                                  end:{x:15.0,y:0.0,z:10.0},
                                  radius:1.0,
                                  color:'black',
                                  fromCap:false,
                                  toCap:false});
                                  
              viewer.zoomTo();    
    viewer.render( );
    cyl2.updateStyle({hidden:true});
    viewer.render( );
    
    viewer._handleMouseDown({pageX: 220, pageY: 265, preventDefault: function(){}});
    viewer._handleMouseUp({pageX: 220, pageY: 265, preventDefault: function(){}});

  };
*/
  /* @div
  <div  class='viewer_3Dmoljs'  style="width: 400px; height: 400px;" data-backgroundColor="white"  data-callback='$3Dmol.testclicky'></div>
*/


