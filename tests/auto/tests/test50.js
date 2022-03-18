
$.get("data/TRPcage.prmtop",  (data)=> {
    const m = viewer.addModel(data, "prmtop");	
    $.get("data/TRPcage.inpcrd", (coordinatesData)=> {
        m.setCoordinates(coordinatesData, "inpcrd");
        viewer.setStyle({},{sphere:{}});
        viewer.zoomTo();
        viewer.render();
      });
    });
