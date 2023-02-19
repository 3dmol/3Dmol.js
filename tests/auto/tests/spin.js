
$.get('data/1lo6.pdb', function(data) {
      viewer.addModel(data,'pdb');
      viewer.setStyle({cartoon:{color:'spectrum'}});
      viewer.zoomTo();
      viewer.spin('z');
      viewer.render( );
      var cnt = 0;
      viewer.setViewChangeCallback(() => {
         cnt += 1;
         if(cnt == 47) {
            viewer.spin(false);
            viewer.render();
         }
      });      
});
