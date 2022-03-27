
//check cif reading and individual styling of duplicated assembly atoms
$.get('data/254385.cif', function(cif) {
    viewer.addModel(cif,'cif',{doAssembly:true,normalizeAssembly:true,duplicateAssemblyAtoms:true});
    viewer.addUnitCell();
    viewer.setStyle({stick:{}});
    viewer.setStyle({sym:0},{stick:{colorscheme:'greenCarbon'}});
    viewer.setStyle({sym:1},{stick:{colorscheme:'magentaCarbon'}});
    viewer.setStyle({sym:2},{stick:{colorscheme:'cyanCarbon'}});
    viewer.addSurface($3Dmol.VDW,{opacity:0.5});
    viewer.zoomTo();
    viewer.render();
});
