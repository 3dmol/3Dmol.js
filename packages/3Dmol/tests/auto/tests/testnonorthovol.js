
$.get('data/h-bn.cube', function(data){
    viewer.addModel(data, "cube");
    viewer.setStyle({}, {stick:{}});
    viewer.zoomTo();
    //$.get('test_structs/benzene-homo.cube', function(data){
      var voldata = new $3Dmol.VolumeData(data, "cube");
      if(viewer.hasVolumetricRender()) {
        viewer.addVolumetricRender(voldata, {transferfn: [{color: "blue", opacity: .075, value: 0.1},
                                                                    {color: "blue", opacity: .001, value: 0.01},
                                                                    {color: "white", opacity: 0, value: 0},
                                                                    {color: "red", opacity: .001, value: -0.01},
                                                                    {color: "red", opacity: .075, value: -0.1}]});
      } else {
        viewer.addIsosurface(voldata, {isoval: 0.01, color: "blue", alpha: 0.95, smoothness: 10});
        viewer.addIsosurface(voldata, {isoval: -0.01, color: "red", alpha: 0.95, smoothness: 10});
      }

      viewer.setStyle({}, {stick:{}});
      viewer.zoomTo();
      viewer.addUnitCell();    

      viewer.render();
    }, 'text');
