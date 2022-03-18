

    $3Dmol.download('mmtf:1mo8',viewer,{},(m) => {
       m.setStyle({'cartoon':{color:'spectrum'}});
       viewer.zoomTo();
       viewer.render();
    });
