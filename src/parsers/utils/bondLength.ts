
/* Covalent radii lookup table used to identify bonds in assignBonds 
 */
export let bondTable = {
  H :0.37,                                                                                                                                He:0.32,
  Li:1.34,Be:0.90,                                                                                B :0.82,C :0.77,N :0.75,O :0.73,F :0.71,Ne:0.69,
  Na:1.54,Mg:1.30,                                                                                Al:1.18,Si:1.11,P :1.06,S :1.02,Cl:0.99,Ar:0.97,
  K :1.96,Ca:1.74,Sc:1.44,Ti:1.56,V :1.25,/* Cr */Mn:1.39,Fe:1.25,Co:1.26,Ni:1.21,Cu:1.38,Zn:1.31,Ga:1.26,Ge:1.22,/* As */Se:1.16,Br:1.14,Kr:1.10,
  Rb:2.11,Sr:1.92,Y :1.62,Zr:1.48,Nb:1.37,Mo:1.45,Tc:1.56,Ru:1.26,Rh:1.35,Pd:1.31,Ag:1.53,Cd:1.48,In:1.44,Sn:1.41,Sb:1.38,Te:1.35,I :1.33,Xe:1.30,
  Cs:2.25,Ba:1.98,Lu:1.60,Hf:1.50,Ta:1.38,W :1.46,Re:1.59,Os:1.44,Ir:1.37,Pt:1.28,Au:1.44,Hg:1.49,Tl:1.48,Pb:1.47,Bi:1.46,/* Po *//* At */Rn:1.45,

  // None of the bottom row or any of the Lanthanides have bond lengths
};


/** Get the length used for bond identification for the specified element.
 * 
 */
export function bondLength(elem) {
  return bondTable[elem] || 1.6;
};

/** Set the length used for bond identification for the specified element.
 * 
 */
export function setBondLength(elem:string, radius:number) {
  if(radius < 0) radius = 0;
  bondTable[elem] = radius;
}

