let socket = null;
let sessionName = null;

const joinSession = function(name) {
    sessionName = name;
    socket.emit('join session event', {
        name : sessionName
    });
};

// setup socket handlers and session related event handlers
const initSessions = function() {
    let initiator = false;
    let joined = false;

    // webserver needs to have appropriate rules to forward to flask
    // https://stackoverflow.com/questions/36472920/apache-proxy-configuration-for-socket-io-project-not-in-root
    if(window.location.hostname === 'localhost') {
      // for debugging on localhost go straight to port to avoid having to setup webserver
      socket = io.connect(`${window.location.hostname}:5000`);
    } else {
      socket = io.connect(window.location.hostname);
    }
    socket.on('connect', () => {
        socket.send('User has connected!');
    });

    // show only the connect button, make sure sidebar is back
    const resetUpperRight = function() {
        $('#sessionbutton').show();
        $('#sessionconnect').hide();
        $('#sessioncontrol').hide();
        $('#sessionmonitor').hide();
        $("#menu").show();
        $("#sidenav").show();
    };

    // open session connection dialog
    $('#sessionbutton').on('click', () => {
        $('#sessionbutton').hide();
        $('#sessionconnect').show();
    });

    const setNumConnections = function(count) {
        $('#userinfo').html(`Users: ${  count}`)
    };
    setNumConnections(0);

    // close dialog and show sessions button
    $('#sessionconnectclose').on('click', () => {
        resetUpperRight();
    });

    // as the user types a session name, check for validity
    $('#session_name_input').on('input', () => {
        const name = $('#session_name_input').val();
        if (name.length > 0) {
            socket.emit('check session name event', {
                name
            });
        } else {
            $('#sessionconnectbutton').prop('disabled', true);
        }
    });

    // update button based on existance of name
    socket.on('check session name response', (msg) => {
        $('#sessionconnectbutton').prop('disabled', false);
        if (msg === 'exists') {
            $('#sessionconnectbutton').html("Join");
        } else {
            $('#sessionconnectbutton').html("Create");
        }
    });

    // handle connecting to a session
    const connectSession = function() {
        sessionName = $('#session_name_input').val();
        if ($('#sessionconnectbutton').text() === "Create") {
            socket.emit('create session event', {
                name : sessionName,
                state : glviewer.getInternalState(),
                view : glviewer.getView()
            });
        } else {
            socket.emit('join session event', {
                name : sessionName,
            });
        }
    };

    $('#sessionconnectbutton').on('click', connectSession);
    $('#session_name_input').on('keyup', (e) => {
        if (e.keyCode === 13) { // enter key
            connectSession();
        }
    });

    // register change callbacks
    const viewUpdateCallback = function(newView) {
        socket.emit('viewer view change event', {
            name : sessionName,
            view : newView
        });
    };

    const stateUpdateCallback = function(newState) {
        socket.emit('viewer state change event', {
            name : sessionName,
            state : newState
        });
    };

    socket.on('create session response',(msg) => {
        if (msg) {
            console.log("session created successfully")
            $('.sessionname').html(sessionName);

            // setup callbacks
            glviewer.setViewChangeCallback(viewUpdateCallback);
            glviewer.setStateChangeCallback(stateUpdateCallback);
            initiator = parseInt(msg);
    
            // close the connection create pane and open the
            // connection monitoring
            $('#sessionbutton').hide();
            $('#sessionconnect').hide();
            $('#sessioncontrol').show();
        } else {
            alert("Session name was already taken/ could not be created. Try Again");
            initiator = 0;
        }
    });

    socket.on('join session response', (msg) => {
        if (msg !== "0") {
            joined = msg;
            // close the connection create pane and open the connection
            // monitoring
            // disable the sidebar
            $('.sessionname').html(sessionName);
            $('#sessionbutton').hide();
            $('#sessionconnect').hide();
            $("#menu").hide();
            $("#sidenav").hide();
            $('#sessionmonitor').show();

        } else
            alert("Session Doesn't Exist")
    });

    const deleteSession = function() {
        if (sessionName)
            socket.emit('delete session event', {
                name : sessionName
            });
        clearResultLabels();
        sessionName = null;
    };

    const leaveSession = function() {
        if (sessionName) {
            socket.emit('leave session event', {
                name : sessionName,
            });
        }
        sessionName = null;
    };

    $('#sessiondestroy').on('click', deleteSession);
    $('#sessionleave').on('click', leaveSession);

    // make sure to unregister a session if user closes window
    $(window).on('beforeunload', () => {
        if (initiator)
            deleteSession();
        else
            leaveSession();
    });

    const sessionEnded = function(reason) {
        joined = false;
        resetUpperRight();
        clientFinishQuery();
    };

    socket.on('delete session response', () => {
        $('#createSession,#joinSession').prop('disabled', false);
        initiator = 0;
        // remove callbacks
        glviewer.setViewChangeCallback(null);
        glviewer.setStateChangeCallback(null);
        resetUpperRight();

    });

    socket.on('leave session response', sessionEnded);
    socket.on('disconnect', (reason) => {
      console.log(`disconnect: ${reason}`);
      // do NOT end session since we will attempt to reconnect
    });

    socket.on('error: restart connection', () => {
        console.log("error: restart connection");
    });
    
    socket.on('error', (error) => {
        console.log(`error: ${error}`);
    });
    
    socket.on("reconnect", () => {
      console.log("Reconnected");
      if(initiator) {
        socket.emit('reconnection', {name: sessionName, secret: initiator, sid: null});
      } else if(joined) {
        socket.emit('reconnection', {name: sessionName, secret: null, sid: joined});        
      }
    });    
    
    socket.on("connect_timeout", () => {
      console.log("Connect timeout");
    });    
    
    socket.on('session count', setNumConnections);

    socket.on('viewer view change response', (newView) => {
        if (!initiator) {
            glviewer.setView(newView);
        }
    });

    socket.on('viewer state change response', (newState) => {
        if (!initiator) {
            glviewer.setInternalState(newState);
        }
    });

    let resultLabel = []; // array of labels

    // remove result labels
    let clearResultLabels = function() {
        if (resultLabel) {
            for (let i = 0; i < resultLabel.length; i++) {
                glviewer.removeLabel(resultLabel[i]);
            }
            resultLabel = [];
        }
    };

    $('#refreshresults').on('click', () => {
        socket.emit('query fetch', {
            name : sessionName,
        });
     });

    $('#askbutton').on('click', () => {
        const val = $('#askbutton').text();
        if (val[0] === 'Q') { // start a query
            // tell server to query clients
            socket.emit('query start', {
                name : sessionName,
                kind : "atoms"
            });
        } else if (val[0] === 'S') { // show results
            $('#askbutton').html('Clear Labels');
            $('#refreshresults').show();            
            socket.emit('query fetch', {
                name : sessionName,
            });
        } else if (val[0] === 'C') { // clear results
            socket.emit('query end', {
                name : sessionName,
            });
        }
    });

    let clientFinishQuery = function () {
        $('#selectmessage').hide();
        glviewer.setHoverable({},false);
        glviewer.setClickable({},false);
        clearQueryLabels();    
    };
    
    let clickedLabels = [];
    let currentHoverLabel = null;
    
    let clearQueryLabels = function() {
      for(let i = 0; i < clickedLabels.length; i++) {
          glviewer.removeLabel(clickedLabels[i]);
      }
      clickedLabels = [];
      if(currentHoverLabel) {
          glviewer.removeLabel(currentHoverLabel);
          currentHoverLabel = null;
      }
    };
    
    // return information about what is currently clicked for the server
    const getClicked = function() {
        const ret = [];
        for(let i = 0; i < clickedLabels.length; i++) {
            const p = clickedLabels[i].stylespec.position;
            ret.push([p.x,p.y,p.z]);
        }
        return ret;
    };
    socket.on('query start response', () => {
        if (initiator) {
            $('#askbutton').html('Show Results');
            $('#responseinfo').html('Responses: 0');
        } else {
            // show message to select atoms, setup callbacks
            $('#selectmessage').show();
            
            // how to label an atom
            const atomLabel = function(atom) {
                let l = '';
                if(atom.resn) l = `${atom.resn}:${atom.atom}`;
                else if(atom.atom) l = atom.atom;
                else l = atom.elem;
                return l;
            };
            
            glviewer.setHoverable({},true,(atom,viewer,event,container) => {
                if(!atom.hoverlabel && !atom.clicklabel) {                    
                    currentHoverLabel = atom.hoverlabel = viewer.addLabel(atomLabel(atom),{position: atom, backgroundColor: 0x800080, backgroundOpacity: 0.7, fontColor:'white'});
                }
            },
            (atom) => { 
                if(atom.hoverlabel) {
                 currentHoverLabel = null;
                 glviewer.removeLabel(atom.hoverlabel);
                 delete atom.hoverlabel;
                }
             }
            );            
            
            glviewer.setClickable({},true, (atom, viewer) => {
               if(atom.clicklabel) {
                   // already clicked, deselect
                   const idx = clickedLabels.indexOf(atom.clicklabel);
                   if(idx >= 0) {
                       clickedLabels.splice(idx,1);
                   }
                   viewer.removeLabel(atom.clicklabel);
                   delete atom.clicklabel;
               } else { // select
                   atom.clicklabel = viewer.addLabel(atomLabel(atom), {position: atom, backgroundColor: 0x800080, backgroundOpacity: 1.0, fontColor: 'white'});
                   clickedLabels.push(atom.clicklabel);
                   if(atom.hoverlabel) {
                       glviewer.removeLabel(atom.hoverlabel);
                       if(currentHoverLabel === atom.hoverlabel) currentHoverLabel = null;
                       delete atom.hoverlabel;
                   }                  
               }
               // send clicked atoms to server
               socket.emit('query update', {
                   name : sessionName,
                   selected : getClicked()
               });               
            });
            
            glviewer.render(); // calling render is required to register clickables - should probably fix that
        }
    });

    socket.on('query end response', () => {
        if (initiator) {
            $('#askbutton').html('Query Atoms');
            $('#responseinfo').html('');
            $('#refreshresults').hide();

            clearResultLabels();
        } else {
            // remove message, callbacks
            clientFinishQuery();
        }
    });

    // update number of responses - result is list of {'cnt': N, 'position':
    // (x,y,z}
    socket.on('query update response', (result) => {
        if (initiator) {
            $('#responseinfo').html(`Responses: ${  result}`);            
        }
    });
    
    socket.on('query fetch response', (queryResult) => {
            // first clear existing labels
            clearResultLabels();
            let max = 1;
            // calc max
            for (let i = 0; i < queryResult.length; i++) {
              const val = parseInt(queryResult[i][1]);
              if(val > max) max = val;
            }

            for (let i = 0; i < queryResult.length; i++) {
                const pos = queryResult[i][0];
                const p = {x: pos[0], y: pos[1], z: pos[2]};
                const val = parseInt(queryResult[i][1]);                
                const l = glviewer.addLabel(val, {position : p,backgroundOpacity:0.5+(val/max)*.5});                
                resultLabel.push(l);                
            }
    });
    
};
