$.get("data/4csv.pdb", function (data) {

    var m = viewer.addModel(data, "pdb");
    $.get("data/4csv.pdb", function (data) {
        m.setCoordinates(data, "pdb");
        viewer.zoomTo();
        viewer.render();
    });
});