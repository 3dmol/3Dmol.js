
//this is from Jacabo Durrant
//I don't yet have away to validate the vrml output, but it should mess up setting of hidden style

var data = "HETATM 2364  C1 A1BO A 501       9.103  27.061  27.339  0.50 40.13           C  \n";
  data += "HETATM 2365  C2 A1BO A 501       8.290  25.757  27.249  0.50 40.25           C  \n";
  data += "HETATM 2366  C3 A1BO A 501       8.877  24.861  26.140  0.50 40.02           C  \n";
  data += "HETATM 2367  C4 A1BO A 501       7.814  24.631  25.053  0.50 40.18           C  \n";
  data += "HETATM 2368  OH A1BO A 501       7.249  23.322  25.194  0.50 40.25           O  ";

  var model = viewer.addModel(data, "pdb");

  model.setStyle(
      {},
      {
          "stick": {
              "colorscheme": "default",
              "hidden": false,
              "radius": 0.3,
          },
          "line": {
              "hidden": true
          }
      },
      true
  );

  viewer.zoomTo();
  viewer.render( );
  
   model.setStyle({},
      {"stick": {
          "hidden": true,
      }},
      true
  );

  viewer.zoomTo();
  viewer.render( );
                    
                    
model.setStyle(
                {},
                {"stick": {
                    "hidden": false,
                }},
                true
            );

  viewer.zoomTo();
  viewer.render( );

  var txt = viewer.exportVRML();
  if(txt.length < 100) {
      throw "Where's the VRML?";
  }
  console.log(txt);                    

model.setStyle(
      {},
      {"stick": {
          "hidden": true,
      }},
      true
  );

  viewer.zoomTo();
  viewer.render();

