
//check cif reading and individual styling of duplicated assembly atoms
$.get('data/icsd_200866.cif', function(cif) {
    viewer.addModel(cif,'cif',{doAssembly:true,normalizeAssembly:true});
    viewer.setStyle({sphere:{scale:.25,colorscheme:{'Cs':'grey','I':'purple','Mg':'brown'}},stick:{}});
    viewer.addUnitCell();

    viewer.zoomTo();
    viewer.render();
});
