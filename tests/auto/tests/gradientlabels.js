$3Dmol.download("pdb:4YCV",viewer,{},function(){
                  
    let grad =  new $3Dmol.Gradient.RWB(40,90,70);
    viewer.setStyle({},{cartoon:{thickness:2.0,colorscheme:{prop:'b',gradient: grad}}});
    viewer.addLabel("50                                   90",{position:{x:0,y:0,z:0},useScreen: true, fontSize: 12, borderThickness: 0.5, borderColor: 'black', backgroundGradient: grad});
    viewer.zoomTo();
    viewer.render();
});
