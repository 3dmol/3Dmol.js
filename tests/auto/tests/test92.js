          $3Dmol.download("pdb:4UAA",viewer,{},()=> {  
            viewer.setStyle({},{stick:{}});
            const resis=viewer.getUniqueValues('resi',{})
            const sel = {resi:resis[0]}
            const sel1 = {resi:resis[resis.length-1]}
          //  console.log(sel)
          //  console.log(sel1)
                viewer.addSphere({center:{resi:147},radius:10.0,color:'red'});
                viewer.addCylinder({start:sel1,
                                  end:{x:20.0,y:0.0,z:0.0},
                                  radius:1.0,
                                  color:'black',
                                  fromCap:false,
                                  toCap:false});
                viewer.addLine({dashed:true,start:{x:0,y:0,z:0},end:sel});

                viewer.addArrow({
                      start: sel1,
                      end: {x:-10.0, y:0.0, z:0.0},
                      radius: 1.0,
                      radiusRadio:1.0,
                      mid:1.0,
                      clickable:true,
                      callback(){
                          this.color.setHex(0xFF0000FF);
                          viewer.render( );
                      }
                  });
               
            viewer.zoomTo();
            viewer.render();
              });   
