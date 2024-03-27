$.get('../test_structs/multi.pdb', function(data){
    viewer.addModelsAsFrames(data,'pdb');
    viewer.setStyle({'stick':{'colorscheme':'greenCarbon'}});
    viewer.zoomTo();
    viewer.animate({'loop': "forward", reps: 1, interval: 1000, 'step': 1});
    viewer.render();
});
