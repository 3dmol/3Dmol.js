let atoms = [{elem: 'C', x: 0, y: 0, z: 0}];
let m = viewer.addModel();
m.addAtoms(atoms);
m.setStyle({},{sphere:{radius:2.0}});
viewer.zoomTo();
viewer.addLabel("X",{screenOffset: {x: 100, y: 100}},{});
viewer.addLabel("Z",{},{});
viewer.render();
