/*jshint esversion: 6 */
/*jshint unused: vars */

$3Dmol.Session = function (hostname) {

    var socket = null;
    this.name = null;
    this.state = null; //Stores the current state
    this.view = null; //Stores the current view
    this.viewer = null;

    var clicked_labels = []; // Contains the List of all the clicked labels
    var current_hover_label = null; // The current label displayed on hover


    var viewChangeCallback = function (new_view) {
        socket.emit('viewer view change event', {
            name: this.name,
            view: new_view
        });
    };

    var stateChangeCallback = function (new_state) {
        socket.emit('viewer state change event', {
            name: this.name,
            state: new_state
        });
    };

    var registerCallbacks = function (viewer) {
        viewer.setViewChangeCallback(viewChangeCallback);
        viewer.setStateChangeCallback(stateChangeCallback);
    };

    var unregisterCallbacks = function (viewer) {
        viewer.setViewChangeCallback(null);
        viewer.setViewChangeCallback(null);
    };
    //Clears the labels clicked or the one which is on hover (query labels)
    var clearQueryLabels = function (viewer) {
        for (let i = 0; i < clicked_labels.length; i++) {
            viewer.removeLabel(clicked_labels[i]);
        }
        clicked_labels = [];
        if (current_hover_label) {
            viewer.removeLabel(current_hover_label);
            current_hover_label = null;
        }
    };

    var result_labels = []; // array of labels (contains all query labels from all users)

    // remove result labels
    var clearResultLabels = function (viewer) {
        if (result_labels) {
            for (let i = 0; i < result_labels.length; i++) {
                viewer.removeLabel(result_labels[i]);
            }
            result_labels = [];
        }
    };

    var setResultLabels = function (viewer, query_result, labelStyle) {
        // first clear existing labels
        clearResultLabels(viewer);
        var max = 1;
        //calc max
        for (let i = 0; i < query_result.length; i++) {
            let val = parseInt(query_result[i][1]);
            if (val > max) max = val;
        }

        for (let i = 0; i < query_result.length; i++) {
            let pos = query_result[i][0];
            let p = {
                x: pos[0],
                y: pos[1],
                z: pos[2]
            };
            let val = parseInt(query_result[i][1]);
            if (!labelStyle) {
                labelStyle.backgroundOpacity = 0.5 + (val / max) * 0.5;
            }
            labelStyle.position = p;
            let l = viewer.addLabel(val, labelStyle);
            result_labels.push(l);
        }
    };

    //return information about what is currently clicked for the server
    var getClicked = function () {
        var ret = [];
        for (let i = 0; i < clicked_labels.length; i++) {
            let p = clicked_labels[i].stylespec.position;
            ret.push([p.x, p.y, p.z]);
        }
        return ret;
    };

    //how to label an atom
    var atomLabel = function (atom) {
        var l = '';
        if (atom.resn) l = atom.resn + ":" + atom.atom;
        else if (atom.atom) l = atom.atom;
        else l = atom.elem;
        return l;
    };

    var hoverCallback = function (atom, viewer, labelStyle) {
        if (!atom.hoverlabel && !atom.clicklabel) {
            labelStyle.position = atom;
            current_hover_label = atom.hoverlabel = viewer.addLabel(atomLabel(atom), labelStyle);
        }
    };
    var unhoverCallback = function (atom) {
        if (atom.hoverlabel) {
            current_hover_label = null;
            this.viewer.removeLabel(atom.hoverlabel);
            delete atom.hoverlabel;
        }
    };

    var clickCallback = function (atom, viewer, labelStyle) {
        if (atom.clicklabel) {
            //already clicked, deselect
            let idx = clicked_labels.indexOf(atom.clicklabel);
            if (idx >= 0) {
                clicked_labels.splice(idx, 1);
            }
            viewer.removeLabel(atom.clicklabel);
            delete atom.clicklabel;
        } else { //select
            labelStyle.position = atom;
            atom.clicklabel = viewer.addLabel(atomLabel(atom), labelStyle);
            clicked_labels.push(atom.clicklabel);
            if (atom.hoverlabel) {
                this.viewer.removeLabel(atom.hoverlabel);
                if (current_hover_label == atom.hoverlabel) current_hover_label = null;
                delete atom.hoverlabel;
            }

            //send clicked atom to server
            socket.emit('query update', {
                name: this.name,
                selected: getClicked()
            });
        }
    };


    this.connect = function () {
        // console.log(this.viewer);
        //webserver needs to have appropriate rules to forward to flask
        //https://stackoverflow.com/questions/36472920/apache-proxy-configuration-for-socket-io-project-not-in-root
        if (hostname == 'localhost') {
            //for debugging on localhost go straight to port to avoid having to setup webserver
            socket = io.connect(hostname + ":5000");
        } else {
            socket = io.connect(hostname);
        }
        return new Promise(function (resolve, reject) {
            socket.on('connect', () => {
                //By default, reconnect is true.
                resolve(socket);
            });
            socket.on('connect_error', (error) => {
                //Called after 20 failed connect attempts
                return reject(Error(error));
            });
            socket.on('connect_timeout', (timeout) => {
                //Called after 20 failed connect attempts
                return reject(Error(timeout));
            });

        });
    };
    this.checkName = function (name) {

        socket.emit('check session name event', {
            name: name
        });

        return new Promise(function (resolve, reject) {
            socket.on('check session name response', (msg) => {
                if (msg == 'exists') {
                    reject(Error(msg));
                } else {
                    resolve();
                }
            });
        });
    };

    this.create = function (name) {
        socket.emit('create session event', {
            name: name,
            state: this.state,
            view: this.view
        });
        this.name = name;
        registerCallbacks(this.viewer);
        return new Promise(function (resolve, reject) {
            socket.on('create session response', (msg) => {
                if (msg) {
                    resolve();
                } else {
                    reject(Error("Session name was already taken/ could not be created. Try Again"));
                }
            });
        });
    };

    this.join = function (name) {

        socket.emit('join session event', {
            name: name
        });
        this.name = name;
        registerCallbacks(this.viewer);
        return new Promise(function (resolve, reject) {
            socket.on('join session response', (msg) => {
                if (msg) {
                    resolve();
                } else {
                    reject(Error("Session does not exist/Could not be joined"));
                }
            });
        });
    };

    this.delete = function () {
        socket.emit('delete session event', {
            name: this.name
        });
        this.name = null;
        unregisterCallbacks(this.viewer);
        clearResultLabels(this.viewer);
        return new Promise(function (resolve) {
            socket.on('delete session response', () => {
                resolve();
            });
        });
    };

    this.leave = function () {
        socket.emit('leave session event', {
            name: this.name
        });
        this.name = null;
        unregisterCallbacks(this.viewer);
        clearQueryLabels(this.viewer);
        return new Promise(function (resolve) {
            socket.on('leave session response', () => {
                resolve();
            });
        });
    };

    this.ask = function (update_function) {

        // tell server to query clients
        socket.emit('query start', {
            name: this.name,
            kind: "atoms"
        });
        this.viewer.setHoverable({}, true, hoverCallback, unhoverCallback);
        this.viewer.setClickable({}, true, clickCallback);

        return new Promise(function (resolve) {
            socket.on('query start response', () => {
                socket.on('query update response', (result) => {
                    update_function(result);
                });
                resolve();
            });

        });

    };

    this.showResults = function (labelStyle) {
        socket.emit('query fetch', {
            name: this.name
        });
        setResultLabels(this.viewer, query_result, labelStyle);
        return new Promise(function (resolve) {
            socket.on('query fetch response', (query_result) => {
                resolve(result_labels);
            });
        });

    };

    this.clearResults = function () {
        socket.emit('query end', {
            name: this.name
        });
        clearResultLabels(this.viewer);
        return new Promise(function (resolve) {
            socket.on('query end response', () => {
                resolve();
            });
        });
    };

};