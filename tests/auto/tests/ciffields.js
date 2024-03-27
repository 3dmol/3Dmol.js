
//Issue #772
$.get('data/6ej1_A_trans.cif', function(cif) {
    viewer.addModel(cif);
    viewer.addStyle({hetflag: true},{stick:{colorscheme:"greenCarbon"}});
    viewer.addStyle({hetflag: false}, {cartoon: {style: 'oval', color: 'white', arrows: true,}});
    viewer.zoomTo();
    viewer.render();
});
