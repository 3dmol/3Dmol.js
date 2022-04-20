/**
 * Convert a base64 encoded string to a Uint8Array
 * @function base64ToArray
 * @param {string} base64 encoded string
 */
 export default function base64ToArray(base64) {
   try{
    const binaryString = window.atob(base64);
    const len = binaryString.length;
    const bytes = new Uint8Array(len);
    for (let i = 0; i < len; i++) {
      bytes[i] = binaryString.charCodeAt(i);
    }
    return bytes;
   } catch(e) {
     console.error(`Failed to decode base64 string: ${base64}\n${e}`);
     throw e;
   }
  };