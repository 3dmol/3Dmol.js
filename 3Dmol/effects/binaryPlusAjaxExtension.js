import $ from "jquery";

// use this transport for "binary" data type
if ($.ajaxTransport) $.ajaxTransport('+binary', (options, originalOptions, jqXHR) => {
    // check for conditions and support for blob / arraybuffer response type
    if (
      window.FormData &&
      ((options.dataType && options.dataType === 'binary') ||
        (options.data &&
          ((window.ArrayBuffer && options.data instanceof ArrayBuffer) ||
            (window.Blob && options.data instanceof Blob))))
    ) {
      return {
        // create new XMLHttpRequest
        send(headers, callback) {
          // setup all variables
          const xhr = new XMLHttpRequest();
          const {url} = options;
          const {type} = options;
          const async = options.async || true;
          // blob or arraybuffer. Default is blob
          const dataType = options.responseType || 'blob';
          let data = options.data || null;
          const username = options.username || null;
          const password = options.password || null;
  
          const xhrret = () => {
            data = {};
            data[options.dataType] = xhr.response;
            // make callback and send data
            callback(xhr.status, xhr.statusText, data, xhr.getAllResponseHeaders());
          };
  
          xhr.addEventListener('load', xhrret);
          xhr.addEventListener('error', xhrret);
          xhr.addEventListener('abort', xhrret);
  
          xhr.open(type, url, async, username, password);
  
          // setup custom headers
          for (const i in headers) {
            xhr.setRequestHeader(i, headers[i]);
          }
  
          xhr.responseType = dataType;
          xhr.send(data);
        },
        abort() {
          jqXHR.abort();
        },
      };
    }
    return null;
  });