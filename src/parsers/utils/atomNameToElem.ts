import { bondTable } from "./bondLength";

// Attempts to infer atomic element from an atom name
export function atomNameToElem(name: string, nothetero: boolean) {
  let elem = name.replace(/ /g, "");
  if (
    elem.length > 0 &&
    elem[0] === "H" &&
    elem !== "Hg" &&
    elem !== "He" &&
    elem !== "Hf" &&
    elem !== "Hs" &&
    elem !== "Ho"
  ) {
    elem = "H"; //workaround weird hydrogen names from MD, note mercury must use lowercase
  }
  if (elem.length > 1) {
    elem = elem[0].toUpperCase() + elem.substring(1).toLowerCase();
    if (bondTable[elem] === undefined) {
      //not a known element, probably should just use first letter
      elem = elem[0];
    } else if (nothetero) {
      if (elem === "Ca") {
        //alpha carbon, not calcium
        elem = "C";
      } else if (elem === "Cd") {
        elem = "C";
      }
    }
  }
  return elem;
}
