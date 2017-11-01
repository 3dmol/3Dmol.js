 $3Dmol.download("pdb:2btf",viewer,{},function(){
        viewer.setStyle({cartoon:{color:'spectrum',colorscheme:'roygb'}});
        viewer.setStyle({hetflag:true},{sphere:{opacity:0.8}});
        viewer.render();
              });
