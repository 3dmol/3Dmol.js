
       $3Dmol.download('pdb:5IRE',viewer,{doAssembly: false},function(m) {
        m.setStyle({chain:'A'},{'cartoon':{color:'spectrum'}});
        m.setStyle({chain:'C'},{'cartoon':{style:'trace',color:'blue'}});
        m.setStyle({chain:'E'},{'cartoon':{tubes:true,arrows:true,color:'green',opacity:0.75}});
        m.setStyle({chain:'B'},{'cartoon':{color:'red',opacity:0.5}});
        m.setStyle({chain:'D'},{'cartoon':{style:'trace',color:'grey',opacity:0.75}});
        m.setStyle({chain:'F'},{'cartoon':{arrows:true,color:'white'}});
       // viewer.addStyle({chain:'B'},{line:{}});
       viewer.zoomTo();
       viewer.render();
    });
