

$.get("data/2water.gro", function (data) {
    var m = viewer.addModel(data, "gro");
    viewer.setStyle({}, { sphere: {} });
    viewer.zoomTo();
    viewer.render();
});
