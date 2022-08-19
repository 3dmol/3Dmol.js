/**
 * VolumetricRenderer style specification
 * @TODO implement pruning of data, should start with box
 * @TODO prop {list} coords - coordinates around which to include data; use viewer.selectedAtoms() to convert an AtomSelectionSpec to coordinates
 * @TODO prop {number} seldist - distance around coords to include data [default = 2.0]
*/
export type VolumetricRendererSpec = Partial<{
  /** list of objects containing @color, @opacity and @value properties to specify color per voxel data value */
  transferfn: { color: unknown; opacity: unknown; value: unknown }[];
  /** number of times to sample each voxel approximately (default 5) */
  subsamples: number;
}>;
