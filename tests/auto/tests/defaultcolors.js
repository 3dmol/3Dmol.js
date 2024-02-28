   var viewer = $3Dmol.createViewer(
     'gldiv', //id of div to create canvas in
     {
       defaultcolors: $3Dmol.elementColors.Jmol
     }
   );
       
    var xyz_data = `1
test
N 0.0 0.0 0.0
1
test
N 0.0 0.0 0.0`;
    viewer.addModelsAsFrames(xyz_data, "xyz");
    viewer.setStyle({}, {stick: {radius: 0.15},sphere: {scale: 1}});
    viewer.zoomTo();
    viewer.zoom(2.5);
    viewer.setViewStyle({style: 'outline', color: 'black', width: 0.02});
    viewer.render();
