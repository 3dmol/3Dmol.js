$.get('../test_structs/multi.cif', function(data){
    viewer.addModelsAsFrames(data,'cif',{'doAssembly':true,'duplicateAssemblyAtoms':true,normalizeAssembly:true});
    viewer.addUnitCell();
    viewer.setStyle({'stick':{'colorscheme':'greenCarbon'}});
    viewer.setFrame(-1);
    viewer.zoomTo();
    viewer.setFrame(0);
    viewer.animate({'loop': "forward", reps: 1, 'step': 1});
    viewer.render();
});
