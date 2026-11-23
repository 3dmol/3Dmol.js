
$.get('data/1lo6.pdb', function(data) {
      viewer.addModel(data,'pdb');
      viewer.setStyle({cartoon:{color:'spectrum'}});
      viewer.zoomTo();
      viewer.spin('z');
      viewer.render( );
      var cnt = 0;
      var rotation = viewer.getView().slice(4, 8);
      viewer.setViewChangeCallback((view) => {
         var nextRotation = view.slice(4, 8);
         // Resize observers can redraw without advancing the spin.
         if(nextRotation.every((value, index) => value === rotation[index])) return;
         rotation = nextRotation;
         cnt += 1;
         if(cnt == 47) {
            viewer.spin(false);
            viewer.setViewChangeCallback(null);
            viewer.render();
         }
      });      
});
