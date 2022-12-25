

var m;
$.get("data/2water.gro", function (data) {
    m = viewer.addModel(data, "gro");
    viewer.setStyle({}, { sphere: {} });
    viewer.zoomTo();
    viewer.render();
});
