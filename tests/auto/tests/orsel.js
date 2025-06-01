$.get('../test_structs/multi.pdb', function(data){
    //Issue #804
    viewer.addModels(data,'pdb');
    viewer.setStyle({model:0},{line:{color:'lightblue'}});
    viewer.setStyle({model:1},'cartoon');
    viewer.addStyle({or:[{resn: 'GLU',model:1}]}, 'sphere');

    viewer.zoomTo();
    viewer.render();
});
