

    $3Dmol.download('mmtf:1mo8',viewer,{},function(m) {
       m.setStyle({'cartoon':{color:'spectrum'}});
       viewer.zoomTo();
       viewer.render();
    });
