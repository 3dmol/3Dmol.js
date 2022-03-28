$3Dmol.download("cid:5",viewer,{},function(){
                  
                  viewer.setStyle({stick:{radius:0.5,singleBonds:true,colorscheme:'greenCarbon',outline:true}});
                  viewer.addLabel("1",{},{index:1});
                  viewer.addLabel("2",{alignment: "center"},{index:2});
                  viewer.addLabel("3",{alignment: "bottomRight"},{index:3});
                  viewer.addLabel("4",{alignment: "centerRight"},{index:4});
                  viewer.addLabel("5",{alignment: "topRight"},{index:5});
                  
                  viewer.render();
              });
