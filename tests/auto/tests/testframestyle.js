 $3Dmol.download("pdb:2n4f",viewer,{multimodel:true, frames:true},function(){
      viewer.setStyle({cartoon:{color:'black'}});
      viewer.zoomTo();
      viewer.setStyle({frame:-1},{cartoon:{color:'blue'}});

      viewer.animate({loop: "forward",reps:1});

      viewer.render();
});


