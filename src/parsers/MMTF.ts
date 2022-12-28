import { base64ToArray } from "../utilities";
import { Matrix4 } from "../WebGL";
import { computeSecondaryStructure } from "./utils/computeSecondaryStructure";
import { processSymmetries } from "./utils/processSymmetries";

interface MMTFobj {
    decode(data: Uint8Array | ArrayBuffer): any;
}
declare var MMTF: MMTFobj;

var fromCharCode = function (charCodeArray) {
    return String.fromCharCode.apply(null, charCodeArray).replace(/\0/g, '');
};

var convertSS = function (val) {
    //convert mmtf code to 3dmol code
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
    if (val == 0 || val == 2 || val == 4) return 'h';
    if (val == 3) return 's';
    return 'c';
};

let mmtfHETATMtypes = new Set([
    "D-SACCHARIDE",
    "D-SACCHARIDE 1,4 AND 1,4 LINKING",
    "D-SACCHARIDE 1,4 AND 1,6 LINKING",
    "L-SACCHARIDE",
    "L-SACCHARIDE 1,4 AND 1,4 LINKING",
    "L-SACCHARIDE 1,4 AND 1,6 LINKING",
    "NON-POLYMER",
    "OTHER",
    "PEPTIDE-LIKE",
    "SACCHARIDE"]);

/** @param bindata - binary UInt8Array buffer or a base64 encoded string
 *  @param ParserOptionsSpec
 *  @category Parsers
*/
export function MMTFparser(bindata, options) {

    var noH = !options.keepH; // suppress hydrogens by default
    var selAltLoc = options.altLoc ? options.altLoc : 'A'; //default alternate location to select if present
    var ignoreStruct = !!options.noSecondaryStructure;
    var computeStruct = !options.noComputeSecondaryStructure;
    //extract symmetries - only take first assembly, apply to all models (ignoring changes for now)
    var noAssembly = !options.doAssembly; // don't assemble by default
    var assemblyIndex = options.assemblyIndex ? options.assemblyIndex : 0;

    if (typeof (bindata) == "string") {
        //assume base64 encoded
        bindata = base64ToArray(bindata);
    } else {
        bindata = new Uint8Array(bindata);
    }

    var mmtfData = MMTF.decode(bindata);

    var atoms: any[][] & Record<string, any> = [[]];
    var modelData: any[] = atoms.modelData = [];

    // setup index counters
    var modelIndex = 0;
    var chainIndex = 0;
    var groupIndex = 0;
    var atomIndex = 0;

    // setup optional fields
    var secStructList = mmtfData.secStructList;
    var insCodeList = mmtfData.insCodeList;
    var sequenceIndexList = mmtfData.sequenceIndexList;
    var bFactorList = mmtfData.bFactorList;
    var altLocList = mmtfData.altLocList;
    var occupancyList = mmtfData.occupancyList;
    var bondAtomList = mmtfData.bondAtomList;
    var bondOrderList = mmtfData.bondOrderList;

    var numModels = mmtfData.numModels;
    if (numModels == 0) return atoms;
    if (!options.multimodel) numModels = 1; //first only
    // hoisted loop variables
    var i, j, k, kl, m, n;



    var symmetries: Matrix4[] = [];
    if (!noAssembly && mmtfData.bioAssemblyList && mmtfData.bioAssemblyList.length > 0) {
        var transforms = mmtfData.bioAssemblyList[assemblyIndex].transformList;
        for (i = 0, n = transforms.length; i < n; i++) {
            var matrix = new Matrix4(transforms[i].matrix);
            matrix.transpose();
            symmetries.push(matrix);
        }
    }
    var unitCell = null as Record<string, number> | null;
    //unit cell info
    if (mmtfData.unitCell) {
        var u = mmtfData.unitCell;
        unitCell = { 'a': u[0], 'b': u[1], 'c': u[2], 'alpha': u[3], 'beta': u[4], 'gamma': u[5] };
    }

    let chainIsPolymer: boolean[] = [];
    mmtfData.entityList.forEach(entity => {
        entity.chainIndexList.forEach(ch => {
            chainIsPolymer[ch] = entity.type == "polymer";
        });
    });
    var bondAtomListStart = 0; //for current model
    //loop over models, 
    for (m = 0; m < numModels; m++) {
        var modelChainCount = mmtfData.chainsPerModel[m];
        var matoms = atoms[atoms.length - 1];
        var serialToIndex: number[] = []; // map to matoms index, needed for noh

        modelData.push({ symmetries: symmetries, cryst: unitCell });
        for (i = 0; i < modelChainCount; ++i) {

            var chainGroupCount = mmtfData.groupsPerChain[chainIndex];
            var chainId = fromCharCode(
                mmtfData.chainIdList.subarray(chainIndex * 4, chainIndex * 4 + 4)
            );
            if (mmtfData.chainNameList) {
                chainId = fromCharCode(
                    mmtfData.chainNameList.subarray(chainIndex * 4, chainIndex * 4 + 4)
                );
            }

            var startGroup = groupIndex;
            var prevSS = '';
            for (j = 0; j < chainGroupCount; ++j) { //over residues (groups)

                var groupData = mmtfData.groupList[mmtfData.groupTypeList[groupIndex]];
                var groupAtomCount = groupData.atomNameList.length;
                var secStruct = 0;
                var secStructBegin = false;
                var secStructEnd = false;

                if (secStructList) {
                    secStruct = secStructList[groupIndex];
                    var sscode = convertSS(secStruct);
                    if (groupIndex == 0 || sscode != prevSS) {
                        secStructBegin = true;
                    }
                    prevSS = sscode;
                    var nextgroup = groupIndex + 1;
                    if (nextgroup >= secStructList.length || convertSS(secStructList[nextgroup] != sscode)) {
                        secStructEnd = true;
                    }
                }
                var insCode = null as string | null;
                if (mmtfData.insCodeList) {
                    insCode = String.fromCharCode(insCodeList[groupIndex]);
                }
                var sequenceIndex = null;
                if (sequenceIndexList) {
                    sequenceIndex = sequenceIndexList[groupIndex];
                }

                var groupId = mmtfData.groupIdList[groupIndex];
                var groupName = groupData.groupName;
                let groupType = groupData.chemCompType;
                var startAtom = atomIndex;
                //note the following is not identical to respecting HETATM records
                //this information isn't available in MMTF.  
                let isHETATM = mmtfHETATMtypes.has(groupType) || !chainIsPolymer[chainIndex];

                for (k = 0; k < groupAtomCount; ++k) {

                    var element = groupData.elementList[k];
                    if (noH && element == 'H') {
                        atomIndex += 1;
                        continue;
                    }

                    var bFactor = '';
                    if (bFactorList) {
                        bFactor = bFactorList[atomIndex];
                    }
                    var altLoc = '';
                    if (altLocList && altLocList[atomIndex]) { //not zero
                        altLoc = String.fromCharCode(altLocList[atomIndex]);
                    }
                    var occupancy = '';
                    if (occupancyList) {
                        occupancy = occupancyList[atomIndex];
                    }

                    if (altLoc != '' && altLoc != selAltLoc && selAltLoc != '*') {
                        atomIndex += 1;
                        continue;
                    }

                    var atomId = mmtfData.atomIdList[atomIndex];
                    var atomName = groupData.atomNameList[k];
                    var atomCharge = 0;
                    if (groupData.atomChargeList) atomCharge = groupData.atomChargeList[k];
                    var xCoord = mmtfData.xCoordList[atomIndex];
                    var yCoord = mmtfData.yCoordList[atomIndex];
                    var zCoord = mmtfData.zCoordList[atomIndex];

                    serialToIndex[atomIndex] = matoms.length;
                    matoms.push({
                        'resn': groupName,
                        'x': xCoord,
                        'y': yCoord,
                        'z': zCoord,
                        'elem': element,
                        'hetflag': isHETATM,
                        'chain': chainId,
                        'resi': groupId,
                        'icode': altLoc,
                        'rescode': groupId + (altLoc != ' ' ? "^" + altLoc : ""), // combo
                        // resi
                        // and
                        // icode
                        'serial': atomId,
                        'altLoc': altLoc,
                        'index': atomIndex,
                        'atom': atomName,
                        'bonds': [],
                        'ss': convertSS(secStruct),
                        'ssbegin': secStructBegin,
                        'ssend': secStructEnd,
                        'bondOrder': [],
                        'properties': { charge: atomCharge, occupancy: occupancy },
                        'b': bFactor,
                    });

                    atomIndex += 1;
                }

                // intra group bonds
                var groupBondAtomList = groupData.bondAtomList;
                for (k = 0, kl = groupData.bondOrderList.length; k < kl; ++k) {
                    var atomIndex1 = startAtom + groupBondAtomList[k * 2];
                    var atomIndex2 = startAtom + groupBondAtomList[k * 2 + 1];
                    var bondOrder = groupData.bondOrderList[k];

                    //I assume bonds are only recorded once
                    var i1 = serialToIndex[atomIndex1];
                    var i2 = serialToIndex[atomIndex2];
                    var a1 = matoms[i1];
                    var a2 = matoms[i2];
                    if (a1 && a2) {
                        a1.bonds.push(i2);
                        a1.bondOrder.push(bondOrder);
                        a2.bonds.push(i1);
                        a2.bondOrder.push(bondOrder);
                    }
                }

                groupIndex += 1;
            }

            //reset for bonds
            groupIndex = startGroup;
            for (j = 0; j < chainGroupCount; ++j) { //over residues (groups)

                groupIndex += 1;

            }

            chainIndex += 1;
        }


        // inter group bonds
        if (bondAtomList) {
            for (let k = bondAtomListStart, kl = bondAtomList.length; k < kl; k += 2) {
                let atomIndex1 = bondAtomList[k];
                let atomIndex2 = bondAtomList[k + 1];
                let bondOrder = bondOrderList ? bondOrderList[k / 2] : 1;

                if (atomIndex1 >= atomIndex) {
                    bondAtomListStart = k;
                    break; //on next model
                }
                //I assume bonds are only recorded once
                let i1 = serialToIndex[atomIndex1];
                let i2 = serialToIndex[atomIndex2];
                let a1 = matoms[i1];
                let a2 = matoms[i2];
                if (a1 && a2) {
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
        computeSecondaryStructure(atoms, options.hbondCutoff);
    }

    return atoms;
};