// this is optimized for proteins where it is assumed connected
// atoms are on the same or next residue

import areConnected from './areConnected';
import assignBonds from './assignBonds';

const standardResidues = new Set([
  'ABU',
  'ACD',
  'ALA',
  'ALB',
  'ALI',
  'ARG',
  'AR0',
  'ASN',
  'ASP',
  'ASX',
  'BAS',
  'CYS',
  'CYH',
  'CYX',
  'CSS',
  'CSH',
  'GLN',
  'GLU',
  'GLX',
  'GLY',
  'HIS',
  'HIE',
  'HID',
  'HIP',
  'HYP',
  'ILE',
  'ILU',
  'LEU',
  'LYS',
  'MET',
  'PCA',
  'PGA',
  'PHE',
  'PR0',
  'PRO',
  'PRZ',
  'SER',
  'THR',
  'TYR',
  'VAL',
  'A',
  '1MA',
  'C',
  '5MC',
  'OMC',
  'G',
  '1MG',
  '2MG',
  'M2G',
  '7MG',
  'OMG',
  'YG',
  'I',
  'T',
  'U',
  '+U',
  'H2U',
  '5MU',
  'PSU',
  'ACE',
  'F0R',
  'H2O',
  'HOH',
  'WAT',
]);

/**
 * @param {import("../../specs").AtomSpec[]}
 *            atomsarray
 * @returns {import('../../specs').AtomSpec[]}
 */
export default function assignPDBBonds(atomsarray) {
  // assign bonds - yuck, can't count on connect records
  const protatoms = [];
  const hetatoms = [];
  let i;
  let n;
  for (i = 0, n = atomsarray.length; i < n; i++) {
    const atom = atomsarray[i];
    atom.index = i;
    if (atom.hetflag || !standardResidues.has(atom.resn || '')) hetatoms.push(atom);
    else protatoms.push(atom);
  }

  assignBonds(hetatoms);

  // sort by resid
  protatoms.sort((a, b) => {
    if (a.chain && b.chain && a.chain !== b.chain) return a.chain < b.chain ? -1 : 1;
    return (a.resi || 0) - (b.resi || 0);
  });

  // for identifying connected residues
  let currentResi = -1;
  let reschain = -1;
  let lastResConnected;

  for (i = 0, n = protatoms.length; i < n; i++) {
    const ai = protatoms[i];

    if (currentResi && ai.resi !== currentResi) {
      if (!ai.resi) throw new Error(`Residue number missing for atom  ${i}`);
      currentResi = ai.resi;
      if (!lastResConnected) reschain++;

      lastResConnected = false;
    }

    ai.reschain = reschain;

    for (let j = i + 1; j < protatoms.length; j++) {
      const aj = protatoms[j];
      if (aj.chain !== ai.chain) break;
      if (aj.resi && ai.resi && (aj.resi - ai.resi > 1)) break;
      if (areConnected(ai, aj)) {
        if (ai.bonds.indexOf(aj.index) === -1) {
          // only add if not already there
          ai.bonds.push(aj.index);
          ai.bondOrder.push(1);
          aj.bonds.push(ai.index);
          aj.bondOrder.push(1);
        }

        if (ai.resi !== aj.resi) lastResConnected = true;
      }
    }
  }

  return atomsarray
}
