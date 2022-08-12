
//check cif reading and individual styling of duplicated assembly atoms
        
$.get('data/254385.cif', function(data) {
      let m = viewer.addModel(data,'cif',{doAssembly:true,duplicateAssemblyAtoms:true,normalizeAssembly:true});
      viewer.setStyle({sphere:{scale:.25},stick:{}});
      viewer.addUnitCell();
      viewer.replicateUnitCell(3);
      viewer.zoomTo();
      viewer.render(callback);
});                  
