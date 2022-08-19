import { ColorSpec } from "colors";
import { Vector3 } from "WebGL";

/**
 * GLShape style specification
 */
export type ShapeSpec = Partial<{
  /** solid color */
  color: ColorSpec;
  /** transparency */
  alpha: number;
  /** draw as wireframe, not surface */
  wireframe: boolean;
  /** if true, do not display object */
  hidden: boolean;
  /** width of line for wireframe rendering **No longer supported by most browsers** */
  linewidth: number;
  /** if true, user can click on object to trigger callback */
  clickable: boolean;
  /** function to call on click */
  callback: (...args: any[]) => void;
  /** if set, only display in this frame of an animation */
  frame: number;
}>;

/**
 * Specification for adding custom shape. Extends {@link ShapeSpec}.
 */
export type CustomShapeSpec = ShapeSpec &
  Partial<{
    /** List of vertex positions */
    vertexArr: Vector3[];
    /** List of normal vectors for each vertex */
    normalArr: Vector3[];
    /** List of triangles to build the custom shape. Each triangle is defined by the indices of 3 vertices in vertexArr, so the array length should be 3 times the number of faces. */
    faceArr: number[];
    /** Either a single color for the whole object or an array specifying the color at each vertex. */
    color: ColorSpec | ColorSpec[];
  }>;

/**
 * Sphere shape specification. Extends {@link ShapeSpec}.
 */
export type SphereSpec = ShapeSpec &
  Partial<{
    /** center of sphere */
    center: Vector3;
    /** radius of sphere */
    radius: number;
  }>;

/**
 * Box shape specification. Extends {@link ShapeSpec}
 */
export type BoxSpec = ShapeSpec &
  Partial<{
    /** bottom corner of box */
    corner: Vector3;
    /** center of box */
    center: Vector3;
    /** width, height, depth of box */
    dimensions: {
      w: number | Vector3;
      h: number | Vector3;
      d: number | Vector3;
    };
  }>;

/**
 * Arrow shape specification.  Extends {@link ShapeSpec}
 */
export type ArrowSpec = ShapeSpec & Partial<{ 
  start: Vector3;
  end: Vector3;
  radius: number;
  color: ColorSpec;
  hidden: boolean;
  /** ratio of arrow base to cylinder (1.618034 default) */
  radiusRatio: number;
  /** relative position of arrow base (0.618034 default) */
  mid: number;
  /** position of arrow base in length units, if negative positioned from end instead of start.  Overrides mid. */
  midpos: number;
}>;

/**
 * Cylinder shape specification.  Extends {@link ShapeSpec}
 */
export type CylinderSpec = ShapeSpec & Partial<{
  start: Vector3;
  end: Vector3;
  radius: number;
  /** 0 for none, 1 for flat, 2 for round */
  fromCap: number;
  /** 0 for none, 1 for flat, 2 for round */
  toCap: number;
  dashed: boolean;
}>;


/**
 * Curve shape specification.  Extends {@link ShapeSpec}
 */
export type CurveSpec = ShapeSpec & Partial<{
  points: Vector3[];
  /** amount of interpolation */
  smooth: number;
  /** radius of curve */
  radius: number;
  /** if an arrow should be drawn at the start */
  fromArrow: boolean;
  /** if an arrow should be drawn at the end */
  toArrow: boolean;
}>;


/**
 * Line shape specification.  Extends {@link ShapeSpec}  (but defaults to wireframe)
 */
export type LineSpec = ShapeSpec & Partial<{
  start: Vector3;
  end: Vector3;
}>;

export type AnyCase<T extends string> =
    string extends T ? string :
    T extends `${infer F1}${infer F2}${infer R}` ? (
        `${Uppercase<F1> | Lowercase<F1>}${Uppercase<F2> | Lowercase<F2>}${AnyCase<R>}`
    ) :
    T extends `${infer F}${infer R}` ? `${Uppercase<F> | Lowercase<F>}${AnyCase<R>}` :
    ""

/**
* File formats supported by 3Dmol.js
* - cdjson,json -> Chemical JSON format
* - cube -> Gaussian cube format
* - gro -> Gromacs topology format, need to add coordinates to resulting model.
* - mcif,cif -> Crystallographic Information File, the successor to PDB that makes you miss the PDB file format
* - mmtf -> Macromolecular Transmission Format, the successor to PDB that is totally awesome
* - mol2 -> Sybyl Mol2 format
* - pdb -> The venerable Protein Data Bank format
* - pqr -> Like PDB but with partial charges which are read into the partialcharge atom property
* - prmtop ->  Amber topology file, must add coordinates
* - sdf -> MDL MOL format, supports multiple models and meta data
* - vasp -> VASP format (CONTCAR, POSCAR)
* - xyz -> XYZ cartesian coordinates format
*/
export type FileFormats = AnyCase<"cdjson" | "json" | "cube" | "gro" | "mcif" | "cif" | "mmtf" | "mol2" | "pdb" | "pqr" | "prmtop" | "sdf" | "vasp" | "xyz">;

