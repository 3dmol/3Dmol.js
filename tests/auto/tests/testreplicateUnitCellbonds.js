
//check cif reading and individual styling of duplicated assembly atoms
        
$.get('data/58093.cif', function(data) {
      let m = viewer.addModel(data,'cif',{doAssembly:true,duplicateAssemblyAtoms:true,normalizeAssembly:true});
      viewer.setStyle({sphere:{scale:.25},stick:{}});
      viewer.addUnitCell();
      viewer.replicateUnitCell(3,2,1,0,true);
      viewer.zoomTo();
      viewer.render(callback);
});                  
