function resetIntersectionShapes(atoms) {
    for (var i = 0; i < atoms.length; i++) {
        var atom = atoms[i];
        if (atom && (atom.clickable || atom.hoverable)) {
            atom.intersectionShape = { sphere: [], cylinder: [], line: [], triangle: [] };
        }
    }
}

module.exports = { resetIntersectionShapes: resetIntersectionShapes };
