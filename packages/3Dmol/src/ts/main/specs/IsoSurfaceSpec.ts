import { ColorSpec } from "colors";
import { Gradient } from "Gradient";

/**
 * Isosurface style specification
 */
export type IsoSurfaceSpec = Partial<{
    /** specifies the isovalue to draw surface at */
    isoval: number;
    /** if true uses voxel style rendering */
    voxel: boolean;
    /** solid color */
    color: ColorSpec;
    /** transparency, between 0 and 1 */
    opacity: number;
    /** draw as wireframe, not surface */
    wireframe: boolean;
    /** width of line for wireframe rendering **No longer supported by most browsers** */
    linewidth: number;
    /** amount to smooth surface (default 1) */
    smoothness: number;
    /** coordinates around which to include data; use viewer.selectedAtoms() to convert an AtomSelectionSpec to coordinates */
    coords: number[][];
    /** distance around coords to include data [default = 2.0] */
    seldist: number;
    /** volumetric data for vertex coloring, can be VolumeData object or raw data if volformat is specified */
    voldata: VolumeData;
    /** coloring scheme for mapping volumetric data to vertex color, if not a Gradient object, show describe a builtin gradient one by providing an object with gradient, min, max, and (optionally) mid fields. */
    volscheme: Gradient;
    /** format of voldata if not a $3Dmol.VolumeData object */
    volformat: string;
    /** if true, user can click on object to trigger callback */
    clickable: boolean;
    /** function to call on click */
    callback: (viewer: GLViewer, event: MouseEvent) => void;
}>;