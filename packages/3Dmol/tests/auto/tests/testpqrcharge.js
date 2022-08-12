/* @data <textarea style="display: none;" id="pqrmol">REMARK
HETATM    1  O   MOL A   1      -2.030   0.000  -0.000  -0.33957515   1.5200 O
HETATM    2  C   MOL A   1      -0.845  -0.000   0.000   0.36693862   1.7000 C
HETATM    3  O   MOL A   1      -0.068   1.100   0.105  -0.19657335   1.5200 O
HETATM    4  C   MOL A   1       1.285   0.755  -0.114  -0.05553065   1.7000 C
HETATM    5  C   MOL A   1       1.285  -0.755   0.114  -0.05533179   1.7000 C
HETATM    6  O   MOL A   1      -0.068  -1.100  -0.105  -0.19670603   1.5200 O
HETATM    7  H   MOL A   1       1.923   1.308   0.588   0.12701495   1.2000 H
HETATM    8  H   MOL A   1       1.569   1.025  -1.145   0.11137415   1.2000 H
HETATM    9  H   MOL A   1       1.569  -1.025   1.146   0.11136246   1.2000 H
HETATM   10  H   MOL A   1       1.923  -1.308  -0.588   0.12702678   1.2000 H
END
</textarea>
*/
let data = $('#pqrmol').val();
let mol = viewer.addModel(data,'pqr');
let grad =  new $3Dmol.Gradient.RWB(-.3,0.3);
mol.setStyle({},{stick:{radius:0.075,colorscheme:{prop:'partialCharge',gradient:grad}},sphere:{radius:0.35,colorscheme:{prop:'partialCharge',gradient:grad}}});
mol.setStyle({elem:'H'},{sphere:{radius:0.25,colorscheme:{prop:'partialCharge',gradient:grad}}},true);
viewer.zoomTo();
viewer.render();
