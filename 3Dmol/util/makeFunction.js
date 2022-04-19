/**
 * 
 * @param {string|Function} callback
 * @returns {Function|null}
 */
export default function makeFunction (callback) {
  // for py3dmol let users provide callback as string
  if (callback && typeof callback == 'string') {
    // eslint-disable-next-line no-eval
    callback = eval(`(${callback})`);
  }
  // report to console if callback is not a valid function
  if (callback && typeof callback === 'function') {
    return callback;
  }
  return null;
};


