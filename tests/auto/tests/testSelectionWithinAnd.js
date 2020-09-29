

$3Dmol.download("cid:3672",viewer,{},function(){
        viewer.setStyle({}, {stick:{hidden:true},sphere:{hidden:true}});
        viewer.setStyle({and: [{or: [{within:{sel:{elem:"O"}, distance: 2.5}}, {elem: "H"}]}, {within: {sel:{elem: "O"}, distance: 5}}]}, {stick:{radius:0.2},sphere:{radius:0.5}});
        viewer.zoomTo();
        viewer.zoom(0.8);
        viewer.render();
});