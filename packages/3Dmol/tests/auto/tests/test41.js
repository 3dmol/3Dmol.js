 var atoms = [{elem: 'C', x: 0, y: 0, z: 0, bonds: [1,2], bondOrder: [1,2]}, {elem: 'O', x: -1.5, y: 0, z: 0, bonds: [0]},{elem: 'O', x: 1.5, y: 0, z: 0, bonds: [0], bondOrder: [2]}];
           
            
            var m = viewer.addModel();
            m.addAtoms(atoms);
            m.setStyle({},{stick:{}});
            viewer.zoomTo();
            viewer.render();
