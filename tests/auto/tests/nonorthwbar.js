
$.get('data/h-bn.cube', function(data){
    viewer.addModel(data, "cube");
    viewer.setStyle({}, {stick:{}});
    viewer.zoomTo();
      var voldata = new $3Dmol.VolumeData(data, "cube");
      if(viewer.hasVolumetricRender()) {
        viewer.addVolumetricRender(voldata, {transferfn: [{color: "blue", opacity: .1, value: 1},
                                                                    {color: "white", opacity: 0.0, value: 0}]});
      } 
      let grad =  new $3Dmol.Gradient.RWB(0,1,0);
      viewer.addLabel("0                                   1",{position:{x:0,y:0,z:0},useScreen: true, fontSize: 12, fontColor: 'black', borderThickness: 0.5, borderColor: 'black', backgroundGradient: grad});
      viewer.setStyle({}, {stick:{}});
      viewer.zoomTo();
      viewer.zoom(5);
      viewer.render();
    }, 'text');
