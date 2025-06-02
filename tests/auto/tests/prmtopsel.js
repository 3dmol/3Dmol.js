
$.get("data/TC5b.prmtop",
  function(data) {
    var m = viewer.addModel(data, "prmtop");
    $.get("data/heat1.mdcrd",
    function(data) {
        m.setCoordinates(data, "mdcrd");
        viewer.addStyle({resn:'PRO'},{stick:{}});
        viewer.addStyle({resi:'1-3'},{sphere:{}});
        viewer.addStyle({resn:'TRP'},{stick:{}});
        viewer.zoomTo();
        viewer.render(callback);
        });   
    });
  
