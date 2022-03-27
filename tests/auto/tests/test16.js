 $3Dmol.download("pdb:2YBB",viewer,{},function(){
                  
                  viewer.setStyle({chain:'A'},{stick:{hidden:true}});
                  viewer.setStyle({resi:91},{sphere:{color:'red',radius:5.0}});
                  viewer.setStyle({bonds:0},{sphere:{color:'fuchsia',radius:5.0}});
                  var atom={};
                  viewer.setStyle({predicate: function(atom){if(atom.chain==="K") return true;}},{sphere:{color: 'olive',radius:5.0}});
                  viewer.setStyle({chain:'B'},{sphere:{color:'blue',radius:5.0}});
                  viewer.render();
              });
