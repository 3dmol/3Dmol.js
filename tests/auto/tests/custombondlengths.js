
//check cif reading and individual styling of duplicated assembly atoms
        
$.get('data/bao3ti.cif', function(data) {
      let l = $3Dmol.bondLength('Ba');
      $3Dmol.setBondLength('Ba',l*0.1);  
      let m = viewer.addModel(data,'cif',{doAssembly:true,duplicateAssemblyAtoms:true,normalizeAssembly:true});
      viewer.setStyle({sphere:{scale:.4,colorscheme:'Jmol'},stick:{colorscheme:'Jmol'}});
      viewer.replicateUnitCell(2,2,2,0,true);
      viewer.zoomTo();
      viewer.render(callback);
});                  
