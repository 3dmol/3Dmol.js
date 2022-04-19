import $ from 'jquery';
import getbin from './getbin';

/**
 * Load a PDB/PubChem structure into existing viewer. Automatically calls 'zoomTo' and 'render' on viewer after loading model
 * @function download
 * @param {string} query - String specifying pdb or pubchem id; must be prefaced with "pdb: " or "cid: ", respectively
 * @param {import("../GLViewer").default} viewer - Add new model to existing viewer
 * @param {Object} options - Specify additional options
 *                           format: file format to download, if multiple are available, default format is pdb
 *                           pdbUri: URI to retrieve PDB files, default URI is http://www.rcsb.org/pdb/files/
 * @param {Function} callback - Function to call with model as argument after data is loaded.
  
 * @return {Promise<import("../GLModel").default>| import("../GLModel").default | null} GLModel, Promise if callback is not provided
 * @example
 viewer.setBackgroundColor(0xffffffff);
       download('pdb:2nbd',viewer,{onemol: true,multimodel: true},function(m) {
        m.setStyle({'cartoon':{colorscheme:{prop:'ss',map:ssColors.Jmol}}});
       viewer.zoomTo();
       viewer.render(callback);
    });
 */
export default function download(query, viewer, options, callback) {
  let type = '';
  let pdbUri = '';
  let mmtfUri = '';
  let uri = '';
  let promise = null;
  const m = viewer.addModel();

  if (query.indexOf(':') < 0) {
    // no type specifier, guess
    if (query.length === 4) {
      query = `pdb:${query}`;
    } else if (!Number.isNaN(query)) {
      query = `cid:${query}`;
    } else {
      query = `url:${query}`;
    }
  }
  if (query.substring(0, 5) === 'mmtf:') {
    pdbUri = options && options.pdbUri ? options.pdbUri : 'https://mmtf.rcsb.org/v1.0/full/';
    query = query.substring(5).toUpperCase();
    uri = pdbUri + query;
    if (options && typeof options.noComputeSecondaryStructure == 'undefined') {
      // when fetch directly from pdb, trust structure annotations
      options.noComputeSecondaryStructure = true;
    }
    promise = new Promise(resolve => {
      getbin(uri).then(
        ret => {
          m.addMolData(ret, 'mmtf', options);
          viewer.zoomTo();
          viewer.render();
          resolve(m);
        },
        () => {
          console.log(`fetch of ${uri} failed.`);
        }
      );
    });
  } else {
    if (query.substr(0, 4) === 'pdb:') {
      type = 'mmtf';
      if (options && options.format) {
        type = options.format; // can override and require pdb
      }

      if (options && typeof options.noComputeSecondaryStructure == 'undefined') {
        // when fetch directly from pdb, trust structure annotations
        options.noComputeSecondaryStructure = true;
      }
      query = query.substring(4).toUpperCase();
      if (!query.match(/^[1-9][A-Za-z0-9]{3}$/)) {
        // eslint-disable-next-line no-alert
        alert('Wrong PDB ID');
        return null;
      }
      if (type === 'mmtf') {
        mmtfUri = options && options.mmtfUri ? options.mmtfUri : 'https://mmtf.rcsb.org/v1.0/full/';
        uri = mmtfUri + query.toUpperCase();
      } else {
        pdbUri = options && options.pdbUri ? options.pdbUri : 'https://files.rcsb.org/view/';
        uri = `${pdbUri + query}.${type}`;
      }
    } else if (query.substring(0, 4) === 'cid:') {
      type = 'sdf';
      query = query.substring(4);
      if (!query.match(/^[0-9]+$/)) {
        // eslint-disable-next-line no-alert
        alert('Wrong Compound ID');
        return null;
      }
      uri = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${query}/SDF?record_type=3d`;
    } else if (query.substring(0, 4) === 'url:') {
      uri = query.substring(4);
      type = uri;
    }

    const handler = ret => {
      m.addMolData(ret, type, options);
      viewer.zoomTo();
      viewer.render();
    };
    promise = new Promise(resolve => {
      if (type === 'mmtf') {
        // binary data
        getbin(uri)
          .then(ret => {
            handler(ret);
            resolve(m);
          })
          .catch(() => {
            // if mmtf server is being annoying, fallback to text
            pdbUri = options && options.pdbUri ? options.pdbUri : 'https://files.rcsb.org/view/';
            uri = `${pdbUri + query}.pdb`;
            console.log('falling back to pdb format');
            $.get(uri, ret => {
              handler(ret);
              resolve(m);
            }).fail(e => {
              handler('');
              resolve(m);
              console.log(`fetch of ${uri} failed: ${e.statusText}`);
            });
          }); // an error msg has already been printed
      } else {
        $.get(uri, ret => {
          handler(ret);
          resolve(m);
        }).fail(e => {
          handler('');
          resolve(m);
          console.log(`fetch of ${uri} failed: ${e.statusText}`);
        });
      }
    });
  }
  if (callback) {
    // @ts-ignore
    promise.then(callback);
    return m;
  }
  return promise;
}
