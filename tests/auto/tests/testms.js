


$3Dmol.download("pdb:1ycr", viewer, {}, function() {
    viewer.addSurface($3Dmol.SurfaceType.MS, {opacity: 0.8, colorscheme: "greenCarbon"},{chain:'B'});
    viewer.render()
});
