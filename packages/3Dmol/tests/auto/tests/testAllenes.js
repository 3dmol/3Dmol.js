var mol1=`11535
  -OEChem-08011811583D

 10  9  0     0  0  0  0  0  0999 V2000
   -1.8156    0.3958    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6365   -0.5182    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6061   -0.1313    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8460    0.2537    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4296    0.2139   -0.8873 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5165    1.4489    0.0009 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4304    0.2128    0.8865 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8434   -1.5844   -0.0004 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0963    1.3081   -0.0012 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.6497   -0.4732    0.0012 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  5  1  0  0  0  0
  1  6  1  0  0  0  0
  1  7  1  0  0  0  0
  2  3  2  0  0  0  0
  2  8  1  0  0  0  0
  3  4  2  0  0  0  0
  4  9  1  0  0  0  0
  4 10  1  0  0  0  0
M  END
> <PUBCHEM_COMPOUND_CID>
11535

> <PUBCHEM_CONFORMER_RMSD>
0.4

> <PUBCHEM_CONFORMER_DIVERSEORDER>
1

> <PUBCHEM_MMFF94_PARTIAL_CHARGES>
7
1 0.14
10 0.15
2 -0.22
3 -0.13
4 -0.24
8 0.15
9 0.15

> <PUBCHEM_EFFECTIVE_ROTOR_COUNT>
0

> <PUBCHEM_PHARMACOPHORE_FEATURES>
2
1 1 hydrophobe
1 4 hydrophobe

> <PUBCHEM_HEAVY_ATOM_COUNT>
4

> <PUBCHEM_ATOM_DEF_STEREO_COUNT>
0

> <PUBCHEM_ATOM_UDEF_STEREO_COUNT>
0

> <PUBCHEM_BOND_DEF_STEREO_COUNT>
0

> <PUBCHEM_BOND_UDEF_STEREO_COUNT>
0

> <PUBCHEM_ISOTOPIC_ATOM_COUNT>
0

> <PUBCHEM_COMPONENT_COUNT>
1

> <PUBCHEM_CACTVS_TAUTO_COUNT>
1

> <PUBCHEM_CONFORMER_ID>
00002D0F00000001

> <PUBCHEM_MMFF94_ENERGY>
0.0497

> <PUBCHEM_FEATURE_SELFOVERLAP>
10.149

> <PUBCHEM_SHAPE_FINGERPRINT>
139733 1 18410856546855120967
20096714 4 9583512135831906210
21015797 1 18047193226197203286
5460574 1 18410856551107883748

> <PUBCHEM_SHAPE_MULTIPOLES>
82.32
2.49
0.74
0.62
0.07
0.02
0
-0.48
0
-0.08
0
0
0
0

> <PUBCHEM_SHAPE_SELFOVERLAP>
140.087

> <PUBCHEM_SHAPE_VOLUME>
55.5

> <PUBCHEM_COORDINATE_TYPE>
2
5
10

$$$$
`
var mol2 = `M0001
  Spartan 08011809273D

  7  6  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    1.3069 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6557    0.6557    1.8766 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6557   -0.6557    1.8766 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000   -1.3069 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6557   -0.6557   -1.8766 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6557    0.6557   -1.8766 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  2  0  0  0  0
  5  6  1  0  0  0  0
  5  7  1  0  0  0  0
  4  5  2  0  0  0  0
M  END
$$$$
`   
   
   var viewers = $3Dmol.createViewerGrid(
     'gldiv', //id of div to create canvas in
     {
       rows: 2,
       cols: 1,
       control_all: true  //mouse controls all viewers
     }
   );
   
     var viewer = viewers[0][0];
     viewer.addModel(mol1,'sdf');
     viewer.setStyle({stick:{}});
     viewer.zoomTo();
     viewer.render( );


     viewer = viewers[1][0];
     viewer.addModel(mol2);
     viewer.setStyle({stick:{}});
     viewer.zoomTo();
     viewer.rotate(90);
     viewer.render( );
  


