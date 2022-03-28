$.get('data/C111tiny.xyz', function(data) {
 viewer.addModel(data,'xyz');
 viewer.setStyle({stick:{}});
 viewer.zoomTo();
 let pos = {x:-0.57470, y:-7.29915, z:1.13071};
 viewer.rotate(90,'z');
 let screen = viewer.modelToScreen(pos);
 let dist = viewer.screenToModelDistance({x:screen.x, y:screen.y},pos);
 if(dist*dist < 0.01) {
     viewer.setStyle({stick:{color:'green'}});
 }
 viewer.render();
});
