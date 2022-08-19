import { ColorSpec } from "colors";
import { Gradient } from "Gradient";

/**
 * @example
 * var setStyles = function(volumedata){
 *  var data = new $3Dmol.VolumeData(volumedata, "cube");
 *  viewer.addSurface("VDW", {opacity:0.85, voldata: data, volscheme: new $3Dmol.Gradient.RWB(-10,10)},{chain:'A'});
 *  viewer.mapAtomProperties($3Dmol.applyPartialCharges);
 *  viewer.addSurface($3Dmol.SurfaceType.SAS, {map:{prop:'partialCharge',scheme:new $3Dmol.Gradient.RWB(-.05,.05)}, opacity:1.0},{chain:'B'});
 *  viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity:0.85,voldata: data, color:'red'},{chain:'C'});
 *  viewer.addSurface($3Dmol.SurfaceType.SAS, {opacity:0.85,voldata: data, colorscheme:'greenCarbon'},{chain:'D'});
 *  viewer.render();
 * };
 * $3Dmol.download("pdb:4DLN",viewer,{},function(){
 *   $.get("data/1fas.cube",setStyles);
 * });
 */
export type SurfaceStyleSpec = Partial<{
    /** sets the transparency: 0 to hide, 1 for fully opaque */
    opacity: number;
    /** element based coloring is a ColorSchemeSpec */
    colorscheme:any;
    /** fixed coloring, overrides colorscheme */
    color: ColorSpec;
    /** volumetric data for vertex coloring, can be VolumeData object or raw data if volformat is specified */
    voldata: any;//VolumeData not ported yet;
    /** coloring scheme for mapping volumetric data to vertex color, if not a Gradient object, show describe a builtin gradient one by providing an object with gradient, min, max, and (optionally) mid fields. */
    volscheme: Gradient;
    /** format of voldata if not a $3Dmol.VolumeData object */
    volformat: string;
    /** specifies a numeric atom property (prop) and color mapping (scheme) such as {@link $3Dmol.Gradient.RWB}.  Deprecated, use colorscheme instead. */
    map: Record<string, unknown>
}>;
