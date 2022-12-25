$3Dmol.download("pdb:3VOV", viewer, {}, function () {
    viewer.setStyle({ chain: 'A' }, { cartoon: { style: 'edged', color: 'blue' } });
    viewer.setStyle({ chain: 'B' }, { cartoon: { thickness: 0, color: 'purple' } });
    viewer.setStyle({ chain: 'C' }, { cartoon: { color: 'white' } });
    viewer.setStyle({ chain: 'D' }, { cartoon: { thickness: 2.0, tubes: 'true' } });
    viewer.zoom(2);
    viewer.rotate(90);
    viewer.render();
});
