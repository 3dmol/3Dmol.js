

$3Dmol.download("cid:3672",viewer,{},function(){
        viewer.setStyle({},{stick:{radius:0.2},sphere:{radius:0.5}});
        viewer.addStyle({elem:"H",within:{sel:{elem:"C"},distance: 1.5}},{stick:{hidden:true},sphere:{hidden:true}});
        viewer.zoomTo();    
        viewer.render();
});
