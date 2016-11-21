

viewer.setViewStyle({style:"outline"});
$.get('volData/example.pse', function(data){
    viewer.addModel(data, "pse");
    viewer.render();
    viewer.zoomTo();
 });