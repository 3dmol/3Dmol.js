export function isEmpty(obj: { [x: string]: { [x: string]: any; }; hasOwnProperty?: any; }) {
  var name;
  for (name in obj) {
    return false;
  }
  return true;
}
