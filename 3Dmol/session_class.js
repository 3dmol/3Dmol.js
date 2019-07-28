
$3Dmol.Session = function(json) {
    this.name = json.name;
    this.state = null;
    this.view = null;
    this.count = 0;
    this.socket = null;

    this.initSocket = function(){
    //webserver needs to have appropriate rules to forward to flask
    //https://stackoverflow.com/questions/36472920/apache-proxy-configuration-for-socket-io-project-not-in-root
    if(window.location.hostname == 'localhost') {
        //for debugging on localhost go straight to port to avoid having to setup webserver
        this.socket = io.connect(window.location.hostname+":5000");
      } else {
        this.socket = io.connect(window.location.hostname);
      }
    };

    this.connect = function(callback) {
        this.socket.on('connect', callback);
    };

    this.disconnect = function(callback){
        this.socket.on('disconnect', callback);
    };

    this.do = function(event_name, json, callback) {
        this.socket.emit(event_name, json);
        if(typeof(callback) === 'function') {
            callback();
        }
    };

    this.viewUpdateCallback = function(new_view) {
        this.socket.emit('viewer view change event', {           //Replace 'viewer view change event' with event_name
            name : this.name,
            view : new_view
        });
    };

    this.stateUpdateCallback = function(new_state) {        //Replace 'viewer state change event' with event_name
        this.socket.emit('viewer state change event', {
            name : this.name,
            state : new_state
        });
    };
   
};