 $3Dmol.download("pdb:4c7j",viewer,{},function(){
        viewer.setStyle({cartoon:{color:'spectrum',colorscheme:'roygb'}});
        //overly complex selections
        viewer.setStyle({and:[{hetflag:true},{not:{bonds:0}}]}, {stick:{opacity:0.8}});
        viewer.setStyle({and:[{hetflag:true},{bonds:0}]},{sphere:{opacity:0.4}});
       viewer.zoomTo({chain:'A',resi:1286});
       viewer.rotate(90);
        viewer.render();
              });
