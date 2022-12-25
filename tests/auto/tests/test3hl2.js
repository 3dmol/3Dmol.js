
//there are long lines in this pdb

$3Dmol.download("pdb:3hl2", viewer, {}, function() {
    viewer.setStyle({});
    viewer.setStyle({chain:'B'},{cartoon:{color:'spectrum'}});
    viewer.zoomTo({chain:'B'});
    viewer.render();
});
