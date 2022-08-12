/* @script
    $3Dmol.init1YCR = function(viewer) {
        viewer.addStyle({}, {"line":{}});
        viewer.addSurface($3Dmol.SurfaceType.VDW, {'opacity':0.8,colorscheme:'whiteCarbon'}, {"chain":"A"},{"chain":"A"});
        viewer.addStyle({"chain":"B"}, {"stick":{colorscheme:'Jmol'}});
        viewer.zoomTo({"chain":"B"});
        viewer.render();
    };
*/

/*

@div
        <div style="width: 400px; height: 400px;" id='viewer3' class='viewer_3Dmoljs' data-href='data/1YCR.pdb' data-backgroundcolor='0xffffff' data-style='{"cartoon":{}}' data-callback='$3Dmol.init1YCR'></div>  


*/
