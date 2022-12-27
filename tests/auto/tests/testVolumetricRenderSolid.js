/* @data
<textarea id="benzene" style="display: none;">

 OpenBabel05131411233D

 12 12  0  0  0  0  0  0  0  0999 V2000
    0.3122    1.3530   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3278    0.4060   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0156   -0.9469   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3123   -1.3529   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3278   -0.4061    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0155    0.9469    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5536    2.3986   -0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.3541    0.7199    0.0001 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.8004   -1.6788    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5536   -2.3986    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3541   -0.7198   -0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8005    1.6787    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  7  1  0  0  0  0
  1  6  2  0  0  0  0
  2  3  2  0  0  0  0
  2  1  1  0  0  0  0
  2  8  1  0  0  0  0
  3  4  1  0  0  0  0
  3  9  1  0  0  0  0
  4  5  2  0  0  0  0
  4 10  1  0  0  0  0
  5  6  1  0  0  0  0
  6 12  1  0  0  0  0
 11  5  1  0  0  0  0
M  END
$$$$
        </textarea>
*/

$3Dmol.get('data/solid.cube', function(data){  // --iso 0.000001 or 0.1 / data: -0.2->0.2
    var voldata = new $3Dmol.VolumeData(data, "cube");
    
    viewer.addVolumetricRender(voldata, {
        transferfn:[
            { color: "#0000ff", opacity: 0.5, value: -0.2 }, //-0.03
            { color: "#ff00ff", opacity: 0.5, value: 0.2 }, // 1
        ]
        // coords: [{x: 0, y: 0, z: 0}], 
        // seldist: 1.7
    }); 
    
    var rec1 = viewer.addModel($('#benzene').val(), "sdf");
    rec1.setStyle({stick:{color:'lightgray', opacity:'1.0'}, sphere:{color:'gray', radius: 0.4}});
    viewer.addBox({color: "grey", wireframe: true, corner: {x: -4.774, y: -4.7993, z: -3.4464}, 
      dimensions: {w: 9.702, h: 9.702, d: 7.055}});
    viewer.zoomTo();
    viewer.rotate(180);
    viewer.render();
});
