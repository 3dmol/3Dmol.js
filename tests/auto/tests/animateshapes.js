
for(var i = 0; i < 10; i++) {
    viewer.addBox({center:{x:i,y:0,z:i},dimensions: {w:3,h:3,d:2},color:'darkblue',frame: i});
}

viewer.zoomTo();
viewer.animate({loop:"forward",reps:1});
viewer.render();
