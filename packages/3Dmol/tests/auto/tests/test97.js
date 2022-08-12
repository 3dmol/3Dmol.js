


$3Dmol.download("pdb:4uhx", viewer, {}, function() {
    viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity: 0.8, color: "white"});
    viewer.render()
});
