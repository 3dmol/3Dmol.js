/*jshint esversion: 6 */
/*jshint unused: vars */

$3Dmol.Session = function (hostname) {

    var socket = null;
    this.name = null;
    this.state = null; //Stores the current state
    this.view = null; //Stores the current view

    this.connect = function () {
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
            socket.on('check session name response', function (msg) {
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
        return new Promise(function (resolve, reject) {
            socket.on('create session response', function (msg) {
                if (msg) {
                    this.name = name;
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
        return new Promise(function (resolve, reject) {
            socket.on('join session response', function (msg) {
                if (msg) {
                    this.name = name;
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
        return new Promise(function (resolve) {
            socket.on('delete session response', function () {
                    this.name = null;
                    resolve();
            });
        });
    };
    this.leave = function () {
        socket.emit('leave session event', {
            name: this.name
        });
        return new Promise(function (resolve) {
            socket.on('leave session response', function () {
                    this.name = null;
                    resolve();
            });
        });
    };

};