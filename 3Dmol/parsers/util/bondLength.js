import bondTable from "./bondTable";

export default function bondLength(elem) {
    return bondTable[elem] || 1.6;
  };