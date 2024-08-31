

let l = viewer.addLabel("Hello World",{position:{x:0,y:0,z:1},useScreen: true, fontSize: 64});
let l2 = viewer.addLabel("Hello World",{position:{x:10,y:0,z:0},backgroundColor: "blue",useScreen: true, fontSize: 64});
let l3 = viewer.addLabel("Hello World",{position:{x:10,y:100,z:0},backgroundColor: "yellow",useScreen: true, fontSize: 64});

l.hide();
l3.hide();
viewer.render( );
l3.show();
viewer.render();
