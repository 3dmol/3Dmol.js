
$.get('data/h-bn.cube', function(data){
      viewer.addModel(data, "cube");
      viewer.setStyle({}, {stick:{}});
      viewer.addUnitCell();    
      var voldata = new $3Dmol.VolumeData(data, "cube");
      viewer.addIsosurface(voldata, {isoval: 0.066, color: "blue", alpha: 0.95, smoothness: 10});
      viewer.rotate(90);
      viewer.zoomTo();
      viewer.render();
    }, 'text');
