
//check unit cell parsing of xyz (who knew?)
$.get('data/empty-V2CO2.xyz', function(xyz) {
    viewer.addModel(xyz,'xyz');
    viewer.setStyle({sphere:{radius:0.5}});
    viewer.addUnitCell();
    viewer.zoomTo();
    viewer.render();
});
