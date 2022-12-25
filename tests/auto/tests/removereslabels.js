$3Dmol.download("mmtf:2ll5",viewer,{},function(){
    viewer.setStyle({stick:{radius:0.15},cartoon:{}});
    var labels = viewer.addResLabels({hetflag:false}, {font: 'Arial', fontColor:'black',backgroundColor:'yellow',showBackground:true});
    for(let i = 0; i < labels.length; i++) {
      if(i%2) viewer.removeLabel(labels[i]);
    }
    viewer.zoomTo();
    viewer.render();                  
  });
