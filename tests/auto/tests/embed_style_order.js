// Suffixes are sorted as strings: style1, style10, then style2.
// Both carbon atoms are green spheres; the oxygen remains a red sphere.
/* @data
<textarea id="inline-style-order-data" style="display:none">3
Inline style ordering
C -2 0 0
C 0 0 0
O 2 0 0
</textarea>
*/
/* @div
<div id="inline-style-order"
     class="viewer_3Dmoljs"
     style="width:400px;height:400px;position:relative"
     data-element="inline-style-order-data" data-type="xyz"
     data-backgroundcolor="black"
     data-config='{"backgroundColor":"white"}'
     data-style="sphere:radius=0.6,color=red"
     data-select1="elem:C" data-style1="stick:radius=0.3,color=blue"
     data-select10="elem:C" data-style10="sphere:radius=0.6,color=orange"
     data-select2="elem:C" data-style2="sphere:radius=0.6,color=green"
     data-callback="window.checkInlineStyleOrder"></div>
*/
/* @script
window.checkInlineStyleOrder = function(viewer) {
    delete window.checkInlineStyleOrder;
    var carbons = viewer.selectedAtoms({elem: "C"});
    var oxygen = viewer.selectedAtoms({elem: "O"});
    if (carbons.length !== 2 || oxygen.length !== 1 ||
        carbons.some(function(atom) { return !atom.style.sphere || atom.style.sphere.color !== "green" || atom.style.stick; }) ||
        oxygen[0].style.sphere.color !== "red" || viewer.getConfig().backgroundColor !== "white") {
        throw new Error("Inline loading, lexical style order, or config precedence changed");
    }
    viewer.addLabel("style1 < style10 < style2", {useScreen: true,
        position: {x: 12, y: 12, z: 0}, fontSize: 16, backgroundColor: "black"});
    viewer.render();
};
*/
