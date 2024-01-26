export function isEmpty(obj: {
  [x: string]: { [x: string]: unknown };
  hasOwnProperty?: any;
}) {
  for (const _ in obj) {
    return false;
  }
  return true;
}
