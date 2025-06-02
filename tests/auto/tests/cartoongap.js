$3Dmol.download("pdb:7ksa", viewer, {}, function () {
    viewer.setStyle({cartoon: {color:"spectrum",gapcutoff:11}});
    viewer.zoomTo({resi:[208,218,258,269,278,287]});
    viewer.render();
});
