// return the value of an atom property prop, or null if non existent
// looks first in properties, then in the atom itself
export default function getAtomProperty(atom, prop) {
    let val = null;
    if (atom.properties &&
            typeof (atom.properties[prop]) != "undefined") {
        val = atom.properties[prop];
    } else if(typeof(atom[prop]) != 'undefined') {
        val = atom[prop];
    }
    return val;
};