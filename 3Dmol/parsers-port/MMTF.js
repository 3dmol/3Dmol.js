import base64ToArray from '../util/base64ToArray';
import {Matrix4} from '../WebGL/math';
import computeSecondaryStructure from './util/computeSecondaryStructure';
import processSymmetries from './util/processSymmeties';

const fromCharCode = charCodeArray =>
  String.fromCharCode.apply(null, charCodeArray).replace(/\0/g, '');

const convertSS = function (val) {
  // convert mmtf code to 3dmol code
  /*    
        0:  pi helix
        1:  bend
        2:  alpha helix
        3:  sheet extended
        4:  3-10 helix
        5:  bridge
        6:  turn
        7:  coil
       */
  if (val === 0 || val === 2 || val === 4) return 'h';
  if (val === 3) return 's';
  return 'c';
};

const mmtfHETATMtypes = new Set([
  'D-SACCHARIDE',
  'D-SACCHARIDE 1,4 AND 1,4 LINKING',
  'D-SACCHARIDE 1,4 AND 1,6 LINKING',
  'L-SACCHARIDE',
  'L-SACCHARIDE 1,4 AND 1,4 LINKING',
  'L-SACCHARIDE 1,4 AND 1,6 LINKING',
  'NON-POLYMER',
  'OTHER',
  'PEPTIDE-LIKE',
  'SACCHARIDE',
]);

/**
 * mmtf shoul be passed as a binary UInt8Array buffer or a base64 encoded string
 * @param {string|ArrayBuffer} bindata 
 * @param {import('../specs').ParserOptionsSpec} options 
 * @returns 
 */
export default function parseMMTF(bindata, options) {
  const noH = !options.keepH; // suppress hydrogens by default
  const selAltLoc = options.altLoc ? options.altLoc : 'A'; // default alternate location to select if present
  const ignoreStruct = !!options.noSecondaryStructure;
  const computeStruct = !options.noComputeSecondaryStructure;
  // extract symmetries - only take first assembly, apply to all models (ignoring changes for now)
  const noAssembly = !options.doAssembly; // don't assemble by default
  const assemblyIndex = options.assemblyIndex ? options.assemblyIndex : 0;

  if (typeof bindata === 'string') {
    // assume base64 encoded
    bindata = base64ToArray(bindata);
  }

  const mmtfData = MMTF.decode(bindata);

  /** @type {import("../specs").ParserResult} */
  const atoms = [[]];
  const modelData = (atoms.modelData = []);

  // setup index counters
  let modelIndex = 0;
  let chainIndex = 0;
  let groupIndex = 0;
  let atomIndex = 0;

  // setup optional fields
  const {secStructList} = mmtfData;
  const {insCodeList} = mmtfData;
  const {sequenceIndexList} = mmtfData;
  const {bFactorList} = mmtfData;
  const {altLocList} = mmtfData;
  const {occupancyList} = mmtfData;
  const {bondAtomList} = mmtfData;
  const {bondOrderList} = mmtfData;

  let {numModels} = mmtfData;
  if (numModels === 0) return atoms;
  if (!options.multimodel) numModels = 1; // first only
  // hoisted loop variables
  let i;
  let j;
  let k;
  let kl;
  let m;
  let n;

  const symmetries = [];
  if (!noAssembly && mmtfData.bioAssemblyList && mmtfData.bioAssemblyList.length > 0) {
    const transforms = mmtfData.bioAssemblyList[assemblyIndex].transformList;
    for (i = 0, n = transforms.length; i < n; i++) {
      const matrix = new Matrix4(transforms[i].matrix);
      matrix.transpose();
      symmetries.push(matrix);
    }
  }
  let unitCell = null;
  // unit cell info
  if (mmtfData.unitCell) {
    const u = mmtfData.unitCell;
    unitCell = {a: u[0], b: u[1], c: u[2], alpha: u[3], beta: u[4], gamma: u[5]};
  }

  const chainIsPolymer = [];
  mmtfData.entityList.forEach(entity => {
    entity.chainIndexList.forEach(ch => {
      chainIsPolymer[ch] = entity.type === 'polymer';
    });
  });
  let bondAtomListStart = 0; // for current model
  // loop over models,
  for (m = 0; m < numModels; m++) {
    const modelChainCount = mmtfData.chainsPerModel[m];
    const matoms = atoms[atoms.length - 1];
    const serialToIndex = []; // map to matoms index, needed for noh

    modelData.push({symmetries, cryst: unitCell});
    for (i = 0; i < modelChainCount; ++i) {
      const chainGroupCount = mmtfData.groupsPerChain[chainIndex];
      let chainId = fromCharCode(mmtfData.chainIdList.subarray(chainIndex * 4, chainIndex * 4 + 4));
      if (mmtfData.chainNameList) {
        chainId = fromCharCode(mmtfData.chainNameList.subarray(chainIndex * 4, chainIndex * 4 + 4));
      }

      const startGroup = groupIndex;
      let prevSS = '';
      for (j = 0; j < chainGroupCount; ++j) {
        // over residues (groups)

        const groupData = mmtfData.groupList[mmtfData.groupTypeList[groupIndex]];
        const groupAtomCount = groupData.atomNameList.length;
        let secStruct = 0;
        let secStructBegin = false;
        let secStructEnd = false;

        if (secStructList) {
          secStruct = secStructList[groupIndex];
          const sscode = convertSS(secStruct);
          if (groupIndex === 0 || sscode !== prevSS) {
            secStructBegin = true;
          }
          prevSS = sscode;
          const nextgroup = groupIndex + 1;
          if (nextgroup >= secStructList.length || convertSS(secStructList[nextgroup] !== sscode)) {
            secStructEnd = true;
          }
        }
        let insCode = null;
        if (mmtfData.insCodeList) {
          insCode = String.fromCharCode(insCodeList[groupIndex]);
        }
        let sequenceIndex = null;
        if (sequenceIndexList) {
          sequenceIndex = sequenceIndexList[groupIndex];
        }

        const groupId = mmtfData.groupIdList[groupIndex];
        const {groupName} = groupData;
        const groupType = groupData.chemCompType;
        const startAtom = atomIndex;
        // note the following is not identical to respecting HETATM records
        // this information isn't available in MMTF.
        const isHETATM = mmtfHETATMtypes.has(groupType) || !chainIsPolymer[chainIndex];

        for (k = 0; k < groupAtomCount; ++k) {
          const element = groupData.elementList[k];
          if (noH && element === 'H') {
            atomIndex += 1;
            continue;
          }

          let bFactor = '';
          if (bFactorList) {
            bFactor = bFactorList[atomIndex];
          }
          let altLoc = '';
          if (altLocList && altLocList[atomIndex]) {
            // not zero
            altLoc = String.fromCharCode(altLocList[atomIndex]);
          }
          let occupancy = '';
          if (occupancyList) {
            occupancy = occupancyList[atomIndex];
          }

          if (altLoc !== '' && altLoc !== selAltLoc && selAltLoc !== '*') {
            atomIndex += 1;
            continue;
          }

          const atomId = mmtfData.atomIdList[atomIndex];
          const atomName = groupData.atomNameList[k];
          let atomCharge = 0;
          if (groupData.atomChargeList) atomCharge = groupData.atomChargeList[k];
          const xCoord = mmtfData.xCoordList[atomIndex];
          const yCoord = mmtfData.yCoordList[atomIndex];
          const zCoord = mmtfData.zCoordList[atomIndex];

          serialToIndex[atomIndex] = matoms.length;
          matoms.push({
            resn: groupName,
            x: xCoord,
            y: yCoord,
            z: zCoord,
            elem: element,
            hetflag: isHETATM,
            chain: chainId,
            resi: groupId,
            icode: altLoc,
            rescode: groupId + (altLoc !== ' ' ? `^${altLoc}` : ''), // combo
            // resi
            // and
            // icode
            serial: atomId,
            altLoc,
            index: atomIndex,
            atom: atomName,
            bonds: [],
            ss: convertSS(secStruct),
            ssbegin: secStructBegin,
            ssend: secStructEnd,
            bondOrder: [],
            properties: {charge: atomCharge, occupancy},
            b: bFactor,
          });

          atomIndex += 1;
        }

        // intra group bonds
        const groupBondAtomList = groupData.bondAtomList;
        for (k = 0, kl = groupData.bondOrderList.length; k < kl; ++k) {
          const atomIndex1 = startAtom + groupBondAtomList[k * 2];
          const atomIndex2 = startAtom + groupBondAtomList[k * 2 + 1];
          const bondOrder = groupData.bondOrderList[k];

          // I assume bonds are only recorded once
          const i1 = serialToIndex[atomIndex1];
          const i2 = serialToIndex[atomIndex2];
          const a1 = matoms[i1];
          const a2 = matoms[i2];
          if (a1 && a2) {
            if (!a1.bonds) a1.bonds = [];
            if (!a2.bonds) a2.bonds = [];
            if (!a1.bondOrder) a1.bondOrder = [];
            if (!a2.bondOrder) a2.bondOrder = [];
            a1.bonds.push(i2);
            a1.bondOrder.push(bondOrder);
            a2.bonds.push(i1);
            a2.bondOrder.push(bondOrder);
          }
        }

        groupIndex += 1;
      }

      // reset for bonds
      groupIndex = startGroup;
      for (j = 0; j < chainGroupCount; ++j) {
        // over residues (groups)

        groupIndex += 1;
      }

      chainIndex += 1;
    }

    // inter group bonds
    if (bondAtomList) {
      for (let k = bondAtomListStart, kl = bondAtomList.length; k < kl; k += 2) {
        const atomIndex1 = bondAtomList[k];
        const atomIndex2 = bondAtomList[k + 1];
        const bondOrder = bondOrderList ? bondOrderList[k / 2] : 1;

        if (atomIndex1 >= atomIndex) {
          bondAtomListStart = k;
          break; // on next model
        }
        // I assume bonds are only recorded once
        const i1 = serialToIndex[atomIndex1];
        const i2 = serialToIndex[atomIndex2];
        const a1 = matoms[i1];
        const a2 = matoms[i2];
        if (a1 && a2) {
          if (!a1.bonds) a1.bonds = [];
          if (!a2.bonds) a2.bonds = [];
          if (!a1.bondOrder) a1.bondOrder = [];
          if (!a2.bondOrder) a2.bondOrder = [];
          a1.bonds.push(i2);
          a1.bondOrder.push(bondOrder);
          a2.bonds.push(i1);
          a2.bondOrder.push(bondOrder);
        }
      }
    }

    if (options.multimodel) {
      if (!options.onemol) atoms.push([]);
    }

    modelIndex += 1;
  }

  if (!noAssembly) {
    for (let n = 0; n < atoms.length; n++) {
      processSymmetries(modelData[n].symmetries, atoms[n], options, modelData[n].cryst);
    }
  }

  if (computeStruct && !ignoreStruct) {
    computeSecondaryStructure(atoms);
  }

  return atoms;
}
