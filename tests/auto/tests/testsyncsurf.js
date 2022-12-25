
$3Dmol.setSyncSurface(true);
$3Dmol.download("pdb:1YCR", viewer, {}, function () {


    viewer.setStyle({ 'cartoon': {} });
    viewer.setStyle({ chain: 'B' }, { 'cartoon': { color: 'blue' } });
    viewer.addSurface("SAS", { opacity: 0.85, color: 'white' }, { chain: 'A' });
    viewer.rotate(90);
    viewer.zoomTo();
    viewer.render();

});
