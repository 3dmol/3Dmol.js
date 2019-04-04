$3Dmol.download("pdb:2ABJ",viewer,{},function(){

    console.log("tooot")

    viewer.setStyle({stick:{colorscheme:'yellowCarbon', opacity: 0.7}});
    viewer.setStyle({chain:'A'},{stick:{colorscheme:'yellowCarbon', opacity:1.0}});
    viewer.setStyle({chain:'D'}, {sphere:{color:'blue', opacity:1.0}});    
    viewer.setStyle({chain:'G'},{cartoon:{color:'green', opacity:1.0}});  
    
    viewer.zoomTo({chain:'G'});
    
    viewer.render();

});
