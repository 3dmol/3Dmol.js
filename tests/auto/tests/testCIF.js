
//check cif reading and individual styling of duplicated assembly atoms
$.get('data/254385.cif', function(cif) {
    viewer.addModel(cif,'cif',{doAssembly:true,duplicateAssemblyAtoms:true});
    viewer.setStyle({sphere:{colorscheme:'Jmol',scale:0.5},stick:{colorscheme:'Jmol'}});
    viewer.addUnitCell();

    viewer.zoomTo();
    viewer.render();
});
