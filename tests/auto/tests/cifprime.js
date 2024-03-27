
//Issue #772
//dump biopython cif generation
$.get('data/5fv3.cif', function(cif) {
    viewer.addModel(cif);
    viewer.addStyle({sphere:{colorscheme:"greenCarbon"}});
    viewer.zoomTo();
    viewer.render();
});
