$3Dmol.download("pdb:3erk",viewer,{},function(){
    viewer.setStyle("cartoon");
    viewer.setStyle({resn:'SB4'},{'stick':{'colorscheme':'purpleCarbon'}});
    viewer.zoomTo({resn:'SB4'});
    viewer.fitSlab({resn:'SB4'});
    viewer.render();
});