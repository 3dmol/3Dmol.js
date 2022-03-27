
    
$.get('data/temp_1_2_28.pdb', function(data) {
    viewer.addModelsAsFrames(data,'pdb');
    viewer.setStyle({cartoon:{},stick:{}});
    viewer.addResLabels({'resn':['ADE','THY','CYT','GUA']},{},true);
    viewer.zoomTo();
    viewer.animate({loop:"forward",reps:1});
    viewer.render();
    });

