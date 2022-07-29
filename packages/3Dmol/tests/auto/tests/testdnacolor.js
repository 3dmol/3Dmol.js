 $3Dmol.download("pdb:1DNS",viewer,{},function(){
  viewer.setStyle('cartoon');
  viewer.setStyle({resn:'DA', atom:'N1'}, {cartoon:{color:'red'}});
  viewer.setStyle({resn:'DG', atom:'N1'}, {cartoon:{color:'green'}});
  viewer.setStyle({resn:'DC', atom:'N3'}, {cartoon:{color:'blue'}});
  viewer.setStyle({resn:'DT', atom:'N3'}, {cartoon:{color:'yellow'}});
  viewer.zoomTo();
  viewer.rotate(90);
  viewer.render();
});
