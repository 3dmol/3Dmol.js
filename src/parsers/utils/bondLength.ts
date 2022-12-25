import { bondTable } from "./bondTable";
export function bondLength(elem) {
  return bondTable[elem] || 1.6;
};