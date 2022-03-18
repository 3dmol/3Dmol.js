
// do not calculate bonds
$.get('../test_structs/jk.pdb', (data)=> {
    viewer.addModel(data, "pdb", {assignBonds: false});
    viewer.setStyle({stick:{},sphere:{radius:0.5}});
    viewer.zoomTo();
    viewer.render();
});
