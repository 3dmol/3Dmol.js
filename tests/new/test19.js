var setStyles = function(volumedata){
                    var data = new $3Dmol.VolumeData(volumedata, "cube");
                    viewer.addSurface("VDW", {opacity:0.85, voldata: data, volscheme: new $3Dmol.Gradient.RWB(-10,10)},{chain:'A'},callback);
                    viewer.mapAtomProperties($3Dmol.applyPartialCharges);
                    viewer.addSurface($3Dmol.SurfaceType.SAS, {map:{prop:'partialCharge',scheme:new $3Dmol.Gradient.RWB(-.05,.05)}, opacity:1.0},{chain:'B'},callback);
                    viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity:0.85,voldata: data, color:'red'},{chain:'C'},callback);
                    viewer.addSurface($3Dmol.SurfaceType.SAS, {opacity:0.85,voldata: data, colorscheme:'greenCarbon'},{chain:'D'},callback);
                    
              viewer.render();
              };
              $3Dmol.download("pdb:4DLN",viewer,{},function(){
                  
                  $.get("volData/1fas.cube",setStyles);
                  viewer.render();
                });