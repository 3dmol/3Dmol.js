$.get('../test_structs/multi.pdb', function(data){
    //Issue #804 + Issue #834
    let models = viewer.addModels(data,'pdb');
    viewer.setStyle({model:0},{line:{color:'lightblue'}});
    viewer.setStyle({model:1},'cartoon');
    viewer.addStyle({or:[{resn: 'GLU',model:[1]},{resn: 'GLY',model:models}]}, 'sphere');
    
    viewer.zoomTo();
    viewer.render();
});
