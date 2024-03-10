
//Issue #770
//make sure atoms are wrapped even from the reference coordinates (identity)
$.get('data/MIL-101.cif', function(cif) {
    let m = viewer.addModel(cif, 'cif', 
      {doAssembly: true, duplicateAssemblyAtoms: true, 
        wrapAtoms: true});
    
    viewer.addStyle('stick');
    viewer.addUnitCell(m);

    viewer.zoomTo();
    viewer.render();
});
