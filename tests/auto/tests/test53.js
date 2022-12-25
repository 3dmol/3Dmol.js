
$.get("data/TC5b.prmtop",
  function(data) {
    var m = viewer.addModel(data, "prmtop");
    $.get("data/heat1.mdcrd",
    function(data) {
        m.setCoordinates(data, "mdcrd");
        viewer.setStyle({},{sphere:{}});
        viewer.zoomTo();
        viewer.animate({loop:"forward",reps:1});
        viewer.render(callback);
        }); //getbin  
    });
  
