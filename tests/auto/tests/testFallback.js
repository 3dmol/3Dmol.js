

$3Dmol.download("cid:2519",viewer,{},function(){
        viewer.setStyle({},{stick:{radius:0.2},sphere:{radius:0.5}});
        viewer.zoomTo();    
        //turn off imposters and instancing
        viewer.zoom(4);
        viewer.render(callback,{supportsAIA:false,supportsImposters:false} );
});
