

$3Dmol.download("cid:2519",viewer,{},function(){
        viewer.setStyle({},{stick:{radius:0.2},sphere:{radius:0.5}});
        viewer.zoomTo();    
        viewer.zoom(4);
        //turn off imposters
        viewer.render(callback,{supportsAIA:true,supportsImposters:false});
});
