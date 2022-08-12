/*
@data
<textarea style="display:none;" id="waterpdb">HETATM    1 H1   HOH A   0       2.924   3.219   1.550  1.00  0.00       B   H  
HETATM    2 H2   HOH A   0       4.098   2.445   0.948  1.00  0.00       B   H  
HETATM    3 O    HOH A   0       3.850   2.984   1.714  1.00  0.00       B   O  
HETATM    4 H1   HOH A   1       1.432   4.761   0.988  1.00  0.00       B   H  
HETATM    5 H2   HOH A   1       0.066   4.168   0.637  1.00  0.00       B   H  
HETATM    6 O    HOH A   1       0.884   3.972   1.119  1.00  0.00       B   O  
HETATM    7 H1   HOH A   2       2.280   0.928   1.171  1.00  0.00       B   H  
HETATM    8 H2   HOH A   2       0.907   1.567   1.386  1.00  0.00       B   H  
HETATM    9 O    HOH A   2       1.335   0.779   1.019  1.00  0.00       B   O  
HETATM   10 H1   HOH A   3       4.694   4.067   3.659  1.00  0.00       B   H  
HETATM   11 H2   HOH A   3       3.377   4.728   4.071  1.00  0.00       B   H  
HETATM   12 O    HOH A   3       3.945   4.569   3.303  1.00  0.00       B   O 
</textarea>
*/

let data = $('#waterpdb').val();
let m = viewer.addModel(data, "pdb", {keepH:true});
viewer.setStyle({'model': -1}, {'sphere':{}});
viewer.zoomTo();
viewer.render();
