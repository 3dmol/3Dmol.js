/**
 * Convert a base64 encoded string to a Uint8Array
 * @function $3Dmol.base64ToArray
 * @param {string} base64 encoded string
 */
export function base64ToArray(base64) {
  var binary_string =  window.atob(base64);
  var len = binary_string.length;
  var bytes = new Uint8Array( len );
  for (var i = 0; i < len; i++)        {
      bytes[i] = binary_string.charCodeAt(i);
  }
  return bytes;
};