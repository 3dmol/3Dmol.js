
$3Dmol.download("cid:3672",viewer,{},function(){
    let sel = {and: [{or: [{within:{sel:{elem:"O"}, distance: 2.5}}, {elem: "H"}]}, {within: {sel:{elem: "O"}, distance: 5}}]};
    viewer.setStyle({}, {stick:{hidden:true},sphere:{hidden:true}});
    viewer.setStyle(sel, {stick:{radius:0.2},sphere:{radius:0.5}});
    if ('__cached_results' in sel.and) {
        // fail the test if the user could access cached results
        viewer.setStyle({}, {stick:{hidden:true},sphere:{hidden:true}});
    }
    viewer.zoomTo();
    viewer.zoom(0.8);
    viewer.render();
});