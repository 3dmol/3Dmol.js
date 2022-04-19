/**
 * $3Dmol surface types
 * @enum {number}
 */
const SurfaceType = {
  VDW: 1,
  MS: 2,
  SAS: 3,
  SES: 4,
};

export const surfaceTypeMap = {
  VDW: SurfaceType.VDW,
  MS: SurfaceType.MS,
  SAS: SurfaceType.SAS,
  SES: SurfaceType.SES,
};

/**
 *
 * @param {number|string} surfaceType
 * @returns {number}
 */
export function normalizeSurfaceType(surfaceType) {
  if (typeof surfaceType == 'string') {
    if (surfaceTypeMap[surfaceType] !== undefined) {
      return surfaceTypeMap[surfaceType];
    }
    console.log(`Surface type : ${surfaceType} is not recognized`);
    return 0;
  }
  
  if (typeof surfaceType === 'number') {
    return surfaceType;
  } 
  
  return SurfaceType.VDW; // default
}

export default SurfaceType;
