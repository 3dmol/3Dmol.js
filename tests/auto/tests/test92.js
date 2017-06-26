
          $3Dmol.download("pdb:4UAA",viewer,{},function(){  
            viewer.setStyle({},{stick:{}});
            viewer.addSphere({center:{},radius:10.0,color:'red'});
         
            viewer.zoomTo();
            viewer.render();
              });   