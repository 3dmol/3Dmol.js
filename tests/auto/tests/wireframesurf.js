
$3Dmol.download("pdb:4DLN",viewer,{},function(){
	viewer.zoomTo({chain:'A'});
        viewer.rotate(-110);

        viewer.addSurface("VDW", {wireframe:true},{chain:'A'});
        viewer.render();
});
