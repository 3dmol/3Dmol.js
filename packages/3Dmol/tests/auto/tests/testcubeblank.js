//make sure cube file with empty comment lines works


$.get("data/LUMO.cube", function(voldata) {
viewer.addModel(voldata, 'cube');
viewer.setStyle({'sphere': {'colorscheme': 'Jmol','scale': 0.3}, 'stick': {'colorscheme': 'Jmol', 'radius': 0.2}});
viewer.addVolumetricData(voldata, "cube", {'isoval': 0.1, 'color': "red", 'opacity': 0.9});
viewer.addVolumetricData(voldata, "cube", {'isoval': -0.1, 'color': "blue", 'opacity': 0.9});
viewer.zoomTo();
viewer.rotate(90);
viewer.render();
});
