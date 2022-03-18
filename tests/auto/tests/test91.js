
$.get("data/model1.prmtop",  (data)=> {
    const m = viewer.addModel(data, "prmtop");	
    $3Dmol.getbin("data/model1_md2.nc", (ret)=> {
        m.setCoordinates(ret, "netcdf");
        viewer.setStyle({},{sphere:{}});
        viewer.zoomTo();
        viewer.animate({loop:"forward",reps:1});
        viewer.render();
    });
});
