 $3Dmol.get("../test_structs/1vy4.cif",function(data){
    viewer.addModel(data,'cif');
    viewer.setStyle({lchain:'B'},{cartoon:{color:'spectrum',colorscheme:'roygb'}});
    viewer.setStyle({chain:'AC'},{cartoon:{color:'spectrum',colorscheme:'sinebow'}});
    viewer.zoomTo({or:[{lchain:'B'},{chain:'AC'}]});
    viewer.render();
});
