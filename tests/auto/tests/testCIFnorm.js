
//check cif reading and individual styling of duplicated assembly atoms
$.get('data/254385.cif', function(cif) {
    viewer.addModel(cif,'cif',{doAssembly:true,normalizeAssembly:true});
    viewer.setStyle({sphere:{colorscheme:'Jmol',scale:0.25},stick:{colorscheme:'Jmol'}});
    viewer.addUnitCell(null,{alabel:'x',blabel:'y',clabel:'z',box:{hidden:true}});

    viewer.zoomTo();
    viewer.render();
});
