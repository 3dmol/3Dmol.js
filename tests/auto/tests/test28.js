$3Dmol.download("cid:5",viewer,{},function(){
                  
                  viewer.setViewStyle({style:"outline",color:"blue",width:0.2});
                  viewer.setStyle({},{stick:{radius:0.5,singleBonds:true,colorscheme:'greenCarbon',outline:true}});
                  viewer.render();
              });
