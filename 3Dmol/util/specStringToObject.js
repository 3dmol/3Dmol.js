// @ts-check
import isNumeric from './isNumeric';

/**
 * Parse a string that represents a style or atom selection and convert it
 * into an object.  The goal is to make it easier to write out these specifications
 * without resorting to json. Objects cannot be defined recursively.
 * ; - delineates fields of the object
 * : - if the field has a value other than an empty object, it comes after a colon
 * , - delineates key/value pairs of a value object
 *     If the value object consists of ONLY keys (no = present) the keys are
 *     converted to a list.  Otherwise a object of key/value pairs is created with
 *     any missing values set to null
 * = OR ~ - separates key/value pairs of a value object, if not provided value is null
 *     twiddle is supported since = has special meaning in URLs
 * @param {string} str
 * @returns {import('../specs').AtomStyleSpec}
 */
export default function specStringToObject(str) {
  if (typeof str == 'object') {
    return str; // not string, assume was converted already
  }
  if (typeof str == 'undefined' || str == null) {
    return str;
  }

  str = str.replace(/%7E/, '~'); // copy/pasting urls sometimes does this
  // convert things that look like numbers into numbers
  const massage = val => {
    if (isNumeric(val)) {
      /**  
       * hexadecimal does not parse as float 
       * disable the radix rules to abuse the 0x prefix parsing of hex numbers
       * */
      // eslint-disable-next-line radix
      if (Math.floor(Number.parseFloat(val)) === Number.parseInt(val)) {
        return Number.parseFloat(val);
      }
      if (val.indexOf('.') >= 0) {
        return Number.parseFloat(val); // ".7" for example, does not parseInt
      }

      return parseInt(val, 10);
    }
    // boolean conversions
    if (val === 'true') {
      return true;
    }
    if (val === 'false') {
      return false;
    }
    return val;
  };

  /** @type {import('../specs').AtomStyleSpec} */
  const ret = {};
  if (str === 'all') return ret;
  const fields = str.split(';');
  for (let i = 0; i < fields.length; i++) {
    const fv = fields[i].split(':');
    const f = fv[0];
    let val = {};
    let vstr = fv[1];
    if (vstr) {
      vstr = vstr.replace(/~/g, '=');
      if (vstr.indexOf('=') !== -1) {
        // has key=value pairs, must be object
        const kvs = vstr.split(',');
        for (let j = 0; j < kvs.length; j++) {
          const kv = kvs[j].split('=', 2);
          val[kv[0]] = massage(kv[1]);
        }
      } else if (vstr.indexOf(',') !== -1) {
        // has multiple values, must list
        val = vstr.split(',');
      } else {
        val = massage(vstr); // value itself
      }
    }
    ret[f] = val;
  }

  return ret;
}
