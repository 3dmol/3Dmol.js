
//Issue #770
//normalize per-atom when duplicate is true
$.get('data/HKUST-1.cif', function(cif) {
    let m = viewer.addModel(cif, 'cif', 
      {doAssembly: true, duplicateAssemblyAtoms: true, 
        wrapAtoms: true});
    
    viewer.addStyle('stick');
    viewer.addUnitCell(m);

    viewer.zoomTo();
    viewer.render();
});
