$3Dmol.download("cid:5975", viewer, {}, function () {
        viewer.setStyle('line');

        viewer.zoomTo();
        viewer.zoom(2);
        viewer.rotate(90,'x');
        viewer.render();
});