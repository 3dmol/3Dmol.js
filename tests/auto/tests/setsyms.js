$3Dmol.download("pdb:1YCR",viewer,{},function(){
    viewer.setStyle('cartoon');
    let m2 = viewer.createModelFrom({chain:'B'},true);
    m2.setStyle({'cartoon':{'color':'red'}});
    viewer.zoomTo();
    viewer.render();
    viewer.addSurface("VDW",{opacity:.5},{model:m2}).then((surfid) => {

        let syms = m2.getSymmetries()
        let mat = syms[0];
        mat.setRotationFromEuler(new $3Dmol.Vector3(Math.PI/2,0,0));
        mat.makeTranslation(-10,-10,10);
        m2.setSymmetries(syms);
        //surfaces have separate symmetries and are potentially split
        //into pieces if the mesh is large
        let surf = viewer.getSurface(surfid); 
        surf.setSymmetries(syms);
        viewer.zoomTo();
        viewer.render( ); // calling setSymmetries will trigger redraw now
    });
});
