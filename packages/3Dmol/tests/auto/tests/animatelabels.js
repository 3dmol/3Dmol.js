
for(var i = 0; i < 10; i++) {
    viewer.addLabel("Hello "+i,{position:{x:i,y:0,z:i}, backgroundColor: 'orange',frame:i});
}

viewer.animate({loop:"forward",reps:1});
viewer.render();
