 $3Dmol.download("pdb:4c7j",viewer,{},function(){
      viewer.setStyle({cartoon:{color:'spectrum',colorscheme:'roygb'}});
      viewer.zoomTo();
      viewer.rotate(90,{x:1,y:1,z:1});
      viewer.rotate(90,"vx");
      viewer.render();
});
