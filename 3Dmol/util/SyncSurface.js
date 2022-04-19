/**
 * Render surface synchronously if true
 * @param {boolean} [SyncSurface=false]
 * @type {boolean}
 * Internet Explorer refuses to allow webworkers in data blobs.  I can find
 * no way of checking for this feature directly, so must do a browser check
 * */
const SyncSurface =  window.navigator.userAgent.indexOf('MSIE ') >= 0 ||
  window.navigator.userAgent.indexOf('Trident/') >= 0;

export default SyncSurface;