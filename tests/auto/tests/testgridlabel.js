var benz=`
     RDKit          3D

  6  6  0  0  0  0  0  0  0  0999 V2000
   -0.9517    0.7811   -0.6622 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2847    1.3329   -0.3121 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2365    0.5518    0.3512 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9517   -0.7811    0.6644 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2847   -1.3329    0.3144 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2365   -0.5518   -0.3489 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
M  END
$$$$
`

var viewers = $3Dmol.createViewerGrid(
 'gldiv', //id of div to create canvas in
 {
   rows: 2,
   cols: 2,
   linked: true,
   control_all: true  //mouse controls all viewers
 }
);
   
var view0 = viewers[0][0];
var view1 = viewers[1][1];

view0.addModel(benz,'sdf')
view0.setStyle({'stick':{}})
view1.addModel(benz,'sdf')
view1.setStyle({'stick':{"colorscheme": "cyanCarbon"}})
view0.setHoverable(
    {},
    true,
    function(atom,viewer,event,container) {
        if(!atom.label) {
        atom.label = viewer.addLabel(atom.atom,{position: atom, backgroundColor: 'mintcream', fontColor:'black'});
        viewer.render( );
        }},
    function(atom,viewer) {
        if(atom.label) {
        viewer.removeLabel(atom.label);
        delete atom.label;
        }
    }
)

view1.setHoverable(
    {},
    true,
    function(atom,viewer,event,container) {
        if(!atom.label) {
        atom.label = viewer.addLabel(atom.atom,{position: atom, backgroundColor: 'cyan', fontColor:'black'});
        viewer.render( );
        }},
    function(atom,viewer) {
        if(atom.label) {
        viewer.removeLabel(atom.label);
        delete atom.label;
        }
    }
)

view0.setClickable({},true,function(atom,viewer) {
    viewer.removeAllShapes();
    viewer.addSphere({center:atom,radius:1.0,color:'red',alpha:0.4});
    viewer.render( );
});

view1.setClickable({},true,function(atom,viewer) {
    viewer.removeAllShapes();
    viewer.addSphere({center:atom,radius:1.0,color:'purple',alpha:0.4});
    viewer.render();
});

view0.zoomTo();
view1.zoomTo();
view0.render( );
view1.render( );
viewers[1][0].render( );
viewers[0][1].render( );

let xy = view0.modelToScreen(view0.models[0].atoms[1]);

view0._handleMouseDown({pageX: xy.x, pageY: xy.y, preventDefault: function(){}});
view0._handleMouseUp({pageX: xy.x, pageY: xy.y, preventDefault: function(){}});

xy = view1.modelToScreen(view1.models[0].atoms[0]);

view1._handleMouseDown({pageX: xy.x, pageY: xy.y, preventDefault: function(){}});
view1._handleMouseUp({pageX: xy.x, pageY: xy.y, preventDefault: function(){}});

if(isNaN(xy.x) || xy.x == 0) {
    let viewer = view1;
    viewer.render( callback); //presumably in coverage
}
