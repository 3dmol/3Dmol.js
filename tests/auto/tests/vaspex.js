

$.get('data/mg.vasp', function(v) {
    let m = viewer.addModel(v, 'vasp');
    
    viewer.addStyle('sphere');

    viewer.zoomTo();
    viewer.render();
});
