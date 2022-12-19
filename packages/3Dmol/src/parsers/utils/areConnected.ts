import { bondLength } from "./bondLength";

/*
 * return true if atom1 and atom2 are probably bonded to each other based on distance alone
 */
export function areConnected(atom1, atom2) {
  var maxsq = bondLength(atom1.elem) + bondLength(atom2.elem);
  maxsq += 0.25; // fudge factor, especially important for md frames, also see 1i3d
  maxsq *= maxsq;

  var xdiff = atom1.x - atom2.x;
  xdiff *= xdiff;
  if (xdiff > maxsq) return false;
  var ydiff = atom1.y - atom2.y;
  ydiff *= ydiff;
  if (ydiff > maxsq) return false;
  var zdiff = atom1.z - atom2.z;
  zdiff *= zdiff;
  if (zdiff > maxsq) return false;

  var distSquared = xdiff + ydiff + zdiff;

  if (isNaN(distSquared)) return false;
  else if (distSquared < 0.5) return false; // maybe duplicate position.
  else if (distSquared > maxsq) return false;
  else if (
    atom1.altLoc != atom2.altLoc &&
    atom1.altLoc != " " &&
    atom2.altLoc != " "
  )
    return false; // don't connect across alternate locations
  else return true;
}
