

$.get('data/1fas.pqr', function(data){
    viewer.addModel(data, "pqr");    
    //reversed color scheme
    viewer.setStyle({'stick':{'colorscheme':{'prop':'partialCharge','gradient':'rwb','min':1,'max':-1,'mid':0}}});
    viewer.zoomTo();
    viewer.render();
});
