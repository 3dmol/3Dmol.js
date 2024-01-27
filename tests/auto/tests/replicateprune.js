
//check cif reading and individual styling of duplicated assembly atoms
        
$.get('data/Al1Cd2_H.cif', function(data) {
      let m = viewer.addModel(data,'cif',{doAssembly:true,duplicateAssemblyAtoms:true,normalizeAssembly:true});
      viewer.setStyle("sphere");
      viewer.addUnitCell();
      viewer.replicateUnitCell(1);
      viewer.zoomTo();
      viewer.render(callback);
});                  
