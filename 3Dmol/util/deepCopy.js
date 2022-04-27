// deep copy, cannot deal with circular refs; undefined input becomes an empty object
// https://medium.com/javascript-in-plain-english/how-to-deep-copy-objects-and-arrays-in-javascript-7c911359b089
/**
 * @template T
 * @param {T} inObject 
 * @returns {T} 
 */
const deepCopy = inObject => {
  let outObject;
  let value;
  let key;

  if (inObject === undefined) {
    // @ts-ignore the generic type is not used
    return {};
  }
  if (typeof inObject != 'object' || inObject == null) {
    return inObject; // Return the value if inObject is not an object
  }

  // Create an array or object to hold the values
  // eslint-disable-next-line prefer-const
  outObject = Array.isArray(inObject) ? [] : {};

  for (key in inObject) {
    value = inObject[key];
    // Recursively (deep) copy for nested objects, including arrays
    outObject[key] = deepCopy(value);
  }

  // @ts-ignore
  return outObject;
};

export default deepCopy;
