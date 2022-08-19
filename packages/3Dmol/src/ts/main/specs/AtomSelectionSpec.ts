import { AtomSpec } from "./AtomSpec";

/**
 * Atom selection object. Used to specify what atoms should be selected.  Can include
 * any field from {@link AtomSpec} in which case atoms must equal the specified value.
 * All fields must match for the selection to hold. If values
 * are provided as a list, then only one value of the list must match.
 * 
 * @example
 * $3Dmol.download("pdb:2EJ0",viewer,{},function(){
 *  viewer.setStyle({chain:'B'},{cartoon:{color:'spectrum'}});
 *  viewer.setStyle({chain:'B',invert:true},{cartoon:{}});
 *  viewer.setStyle({bonds: 0},{sphere:{radius:0.5}}); //water molecules
 *  viewer.setStyle({resn:'PMP',byres:true,expand:5},{stick:{colorscheme:"greenCarbon"}});
 *  viewer.setStyle({resi:["91-95","42-50"]},{cartoon:{color:"green",thickness:1.0}});
 *  viewer.render();
 * });
 */
export type AtomSelectionSpec = AtomSpec &
  Partial<{
    /** a single model or list of models from which atoms should be selected.  Can also specify by numerical creation order.  Reverse indexing is allowed (-1 specifies last added model). */
    model: /**GLModel*/ unknown | number | /** GLmodel[]*/ unknown[] | number[];
    /** overloaded to select number of bonds, e.g. {bonds: 0} will select all nonbonded atoms */
    bonds: number;
    /** user supplied function that gets passed an {@link AtomSpec} and should return true if the atom should be selected */
    predicate: (atom: AtomSpec) => boolean;
    /** if set, inverts the meaning of the selection */
    invert: boolean;
    /** if set, expands the selection to include all atoms of any residue that has any atom selected */
    byres: boolean;
    /** expands the selection to include all atoms within a given distance from the selection */
    expand: number;
    /** intersects the selection with the set of atoms within a given distance from another selection */
    within: WithinSelectionSpec;
    /** take the intersection of the provided lists of {@link AtomSelectionSpec}s */
    and: AtomSelectionSpec[];
    /** take the union of the provided lists of {@link AtomSelectionSpec}s */
    or: AtomSelectionSpec[];
    /** take the inverse of the provided {@link AtomSelectionSpec} */
    not: AtomSelectionSpec;
  }>;

/**
 * Within selection object. Used to find the subset of an atom selection that is within
 * some distance from another atom selection. When added as a field of an {@link AtomSelectionSpec},
 * intersects the set of atoms in that selection with the set of atoms within a given
 * distance from the given {@link AtomSelectionSpec}.
 * @example
 * $3Dmol.download("pdb:2EJ0",viewer,{},function(){
 *  viewer.setStyle({chain: 'A', within:{distance: 10, sel:{chain: 'B'}}}, {sphere:{}});
 *  viewer.render();
 * });// stylizes atoms in chain A that are within 10 angstroms of an atom in chain B
 *
 */
export type WithinSelectionSpec = Partial<{
  /** the distance in angstroms away from the atom selection to include atoms in the parent selection */
  distance: number;
  /** if set, selects atoms not within distance range for intersection */
  invert: boolean;
  /** the selection of atoms against which to measure the distance from the parent atom selection */
  sel: AtomSelectionSpec;
}>;
