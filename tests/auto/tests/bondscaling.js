$3Dmol.download("cid:10321", viewer, {}, function () {
        viewer.setStyle({stick:{doubleBondScaling: 1, tripleBondScaling: 0.5}});

        viewer.zoomTo();
        viewer.render();
});
