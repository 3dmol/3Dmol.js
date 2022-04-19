import getAtomProperty from "./getAtomProperty";

/* get the min and max values of the specified property in the provided
* @function $3Dmol.getPropertyRange
* @param {AtomSpec[]} atomlist - list of atoms to evaluate
* @param {string} prop - name of property 
* @return {Array} - [min, max] values
*/
export default function getPropertyRange(atomlist, prop) {
    let min = Number.POSITIVE_INFINITY;
    let max = Number.NEGATIVE_INFINITY;

    for (let i = 0, n = atomlist.length; i < n; i++) {
        const atom = atomlist[i];
        const val = getAtomProperty(atom, prop);
        
        if(val != null) {
            if (val < min)
                min = val;
            if (val > max)
                max = val;                
        }
    }

    if (!Number.isFinite(min) && !Number.isFinite(max))
        min = max = 0;
    else if (!Number.isFinite(min))
        min = max;
    else if (!Number.isFinite(max))
        max = min;

    return [ min, max ];
};
