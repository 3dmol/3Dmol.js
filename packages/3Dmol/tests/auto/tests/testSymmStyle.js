
//check cif reading and individual styling of duplicated assembly atoms
$.get('data/9002806.cif', function(cif) {
    viewer.addModel(cif,'cif',{doAssembly:true,duplicateAssemblyAtoms:true, dontConnectDuplicatedAtoms: true});
    viewer.setStyle({sphere:{colorscheme:'Jmol',scale:0.5},stick:{colorscheme:'Jmol'}});
    viewer.addUnitCell();

    viewer.setStyle({sym:2},{sphere:{scale:.5,color:'blue'},stick:{color:'cyan'}});

    viewer.zoomTo();
    viewer.render();
});
