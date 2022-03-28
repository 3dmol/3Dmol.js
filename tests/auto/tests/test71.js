$3Dmol.download("pdb:4DM7",viewer,{},function(){
                 //viewer.setStyle({elem:'C'},{line:{hidden:true}});
                 viewer.setStyle({and:[{or:[{chain:'B'},{chain:'A'}]},{elem:'C'}]},{line:{hidden:true}});//redundant because "and" is implicit unless it is inside of an "or" object
                  

                 viewer.render();
                });
