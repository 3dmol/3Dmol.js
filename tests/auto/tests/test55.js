
    
$.get('data/4csv.pdb', (data) => {
      viewer.addModel(data,'pdb');
      viewer.setStyle({cartoon:{},stick:{}});
      viewer.zoomTo({resn:'STI'},1000);
      viewer.render();
    
 });
 

