
//check cif reading and individual styling of duplicated assembly atoms
$.get('data/254385.cif', function(cif) {
    viewer.addModel(cif,'cif',{doAssembly:true,normalizeAssembly:true,duplicateAssemblyAtoms:true,dontConnectDuplicatedAtoms:true});
    viewer.setStyle({sphere:{colorscheme:'Jmol',scale:0.25},stick:{colorscheme:'Jmol'}});
        viewer.setStyle({sym:0},{sphere:{colorscheme:'rasmol',scale:0.25},stick:{colorscheme:'rasmol'}});

    viewer.addUnitCell(null,{box:{},astyle:{hidden:true},bstyle:{hidden:true},cstyle:{hidden:true}});

    viewer.zoomTo();
    viewer.render();
});
