$.get('data/v3000.mol', function(data) {
    viewer.addModel(data,'sdf');
    viewer.setStyle({stick:{}});
    viewer.zoomTo();
    viewer.render();
});
