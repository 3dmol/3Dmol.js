export type AtomStyleSpec = Partial<{
  /** draw bonds as lines */
  line: LineStyleSpec;
  /** draw atoms as crossed lines (aka stars) */
  cross: CrossStyleSpec;
  /** draw bonds as capped cylinders */
  stick: StickStyleSpec;
  /** draw atoms as spheres */
  sphere: SphereStyleSpec;
  /** draw cartoon representation of secondary structure */
  cartoon: CartoonStyleSpec;
  /** invisible style for click handling only */
  clicksphere: ClickSphereStyleSpec;
}>;


