
$.get("volData/TRPcage.prmtop",  function (data){
    var m = viewer.addModel(data, "prmtop");	
    $.get("volData/TRPcage.inpcrd", function(coordinatesData){
        m.setCoordinates(coordinatesData, "inpcrd");
        viewer.setStyle({},{sphere:{}});
        viewer.zoomTo();
        viewer.render();
      });
    });
