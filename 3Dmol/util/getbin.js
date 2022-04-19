import $ from "jquery"
/**
 * Download binary data (e.g. a gzipped file) into an array buffer and provide
 * arraybuffer to callback.
 * @function getbin
 * @param {string} uri - location of data
 * @param {Function} [callback] - Function to call with arraybuffer as argument.
 * @param {string} [request] - type of request
 * @param {any} [postdata] - data to send with request
 * @return {Promise}
 */
 export default function getbin(uri, callback, request, postdata) {
    const promise = new Promise((resolve, reject) => {
      $.ajax({
        url: uri,
        dataType: 'binary',
        method: request || 'GET',
        data: postdata,
        responseType: 'arraybuffer',
        processData: false,
      })
        .done(ret => {
          resolve(ret);
        })
        .fail((e, txt) => {
          console.log(txt);
          reject();
        });
    });
    // @ts-ignore
    if (callback) return promise.then(callback);
    return promise;
  };