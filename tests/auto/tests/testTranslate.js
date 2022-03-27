
// zoomto should not distrupt translation
    $3Dmol.download('pdb:5oi8',viewer,{},function(m) {
       m.setStyle({'cartoon':{color:'spectrum'}});
       viewer.translate(50,-50);
       viewer.zoomTo();       
       viewer.render();
    });
