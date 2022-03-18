$.get('data/v3000.mol', (data) => {
    viewer.addModel(data,'sdf');
    viewer.setStyle({stick:{}});
    viewer.zoomTo();
    viewer.render();
});
