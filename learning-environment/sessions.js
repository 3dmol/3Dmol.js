var socket = null;
var session_name = null;

var joinSession = function(name) {
    session_name = name;
    socket.emit('join session event', {
        name : session_name
    });
};

// Setup socket handlers and session related event handlers
var initSessions = function() {
    var initiator = false;
    var joined = false;
    var pendingJoin = false;
    var pendingCreate = false;

    joinSession = function(name) {
        if (pendingJoin || pendingCreate) return;
        session_name = name;
        pendingJoin = true;
        socket.emit('join session event', { name: name });
    };

    var showViewerControls = function(visible) {
        var app = window.viewerApp;
        if (app && app.state && app.state.ui && app.state.ui.element) {
            app.state.ui.element.hidden = !visible;
        } else {
            $("#menu,#sidenav").toggle(visible);
        }
    };

    var sessionChanged = function(name, following) {
        if (window.viewerApp && typeof window.viewerApp.sessionChanged === 'function') {
            window.viewerApp.sessionChanged(name, following);
        }
    };

    // Local development can serve the viewer separately from Flask on port 5000.
    // Other hosts preserve the page origin unless an embedding supplies a server.
    var sessionServer = window.viewerSessionServer || (window.location.hostname === 'localhost'
        ? window.location.protocol + '//localhost:5000'
        : window.location.origin);
    socket = io.connect(sessionServer);
    socket.on('connect', function() {
        socket.send('User has connected!');
    });

    // Show only the connect button, make sure sidebar is back
    var resetUpperRight = function() {
        $('#sessionbutton').show();
        $('#sessionconnect').hide();
        $('#sessioncontrol').hide();
        $('#sessionmonitor').hide();
        showViewerControls(true);
        sessionChanged(null, false);
    };

    // Open session connection dialog
    $('#sessionbutton').on('click', function() {
        $('#sessionbutton').hide();
        $('#sessionconnect').show();
    });

    var setNumConnections = function(count) {
        $('#userinfo').html('Users: ' + count)
    };
    setNumConnections(0);

    // Close dialog and show sessions button
    $('#sessionconnectclose').on('click', function() {
        if (pendingJoin || pendingCreate) endSession(true);
        else resetUpperRight();
    });

    // as the user types a session name, check for validity
    $('#session_name_input').on('input', function() {
        var name = $('#session_name_input').val();
        if (name.length > 0) {
            socket.emit('check session name event', {
                name : name
            });
        } else {
            $('#sessionconnectbutton').prop('disabled', true);
        }
    });

    // Update button based on existance of name
    socket.on('check session name response', function(msg) {
        if (pendingJoin || pendingCreate || joined || initiator) return;
        $('#sessionconnectbutton').prop('disabled', false);
        if (msg == 'exists') {
            $('#sessionconnectbutton').html("Join");
        } else {
            $('#sessionconnectbutton').html("Create");
        }
    });

    // Handle connecting to a session
    var connectSession = function() {
        if (pendingJoin || pendingCreate) return;
        session_name = $('#session_name_input').val();
        if ($('#sessionconnectbutton').text() == "Create") {
            pendingCreate = true;
            socket.emit('create session event', {
                name : session_name,
                state : Object.assign({}, glviewer.getInternalState(), { frame: glviewer.getFrame() }),
                view : glviewer.getView()
            });
        } else {
            joinSession(session_name);
        }
    };

    $('#sessionconnectbutton').on('click', connectSession);
    $('#session_name_input').on('keyup', function(e) {
        if (e.keyCode == 13) { // enter key
            connectSession();
        }
    });

    // Register change callbacks
    var viewUpdateCallback = function(new_view) {
        if (!initiator || !session_name) return;
        socket.emit('viewer view change event', {
            name : session_name,
            view : new_view
        });
    };

    var stateUpdateCallback = function(new_state) {
        if (!initiator || !session_name) return;
        socket.emit('viewer state change event', {
            name : session_name,
            state : Object.assign({}, new_state, { frame: glviewer.getFrame() })
        });
    };

    socket.on('create session response',function(msg) {
        if (!pendingCreate || !session_name) return;
        pendingCreate = false;
        if (msg) {
            console.log("session created successfully")
            $('.sessionname').text(session_name);

            // Setup callbacks
            glviewer.setViewChangeCallback(viewUpdateCallback);
            glviewer.setStateChangeCallback(stateUpdateCallback);
            initiator = parseInt(msg);
    
            // Close the connection create pane and open the connection monitoring
            $('#sessionbutton').hide();
            $('#sessionconnect').hide();
            $('#sessioncontrol').show();
        } else {
            alert("Session name was already taken/ could not be created. Try Again");
            initiator = 0;
            session_name = null;
        }
    });

    socket.on('join session response', function(msg) {
        if (!pendingJoin || !session_name) return;
        pendingJoin = false;
        if (msg != "0") {
            joined = msg;
            // Close the connection create pane and open the connection monitoring
            $('.sessionname').text(session_name);
            $('#sessionbutton').hide();
            $('#sessionconnect').hide();
            // Disable the sidebar
            sessionChanged(session_name, true);
            showViewerControls(false);
            $('#sessionmonitor').show();

        } else {
            session_name = null;
            resetUpperRight();
            alert("Session Doesn't Exist");
        }
    });

    var endSession = function(notifyServer) {
        var name = session_name;
        var wasInitiator = initiator;
        var wasCreating = pendingCreate;
        initiator = false;
        joined = false;
        pendingJoin = false;
        pendingCreate = false;
        session_name = null;
        if (wasInitiator) {
            glviewer.setViewChangeCallback(null);
            glviewer.setStateChangeCallback(null);
        }
        if (notifyServer && name) {
            socket.emit(wasInitiator || wasCreating ? 'delete session event' : 'leave session event', { name: name });
        }
        clearResultLabels();
        clientFinishQuery();
        resetUpperRight();
    };

    // An embedding can leave a classroom before replacing its molecular scene.
    window.leaveViewerSession = function() {
        if (session_name || initiator || joined || pendingJoin || pendingCreate) endSession(true);
    };
    $('#sessiondestroy').on('click', window.leaveViewerSession);
    $('#sessionleave').on('click', window.leaveViewerSession);

    // Make sure to unregister a session if user closes window
    $(window).on('beforeunload', function() {
        if (session_name) endSession(true);
    });

    var sessionEnded = function(reason) {
        if (joined || initiator) endSession(false);
    };

    socket.on('delete session response', function() {
        if (!initiator) return;
        $('#createSession,#joinSession').prop('disabled', false);
        endSession(false);
    });

    socket.on('leave session response', sessionEnded);
    socket.on('disconnect', function(reason) {
      console.log("disconnect: "+reason);
      // Do NOT end session since we will attempt to reconnect
    });

    socket.on('error: restart connection', function() {
        console.log("error: restart connection");
    });
    
    socket.on('error', (error) => {
        console.log('error: '+error);
    });
    
    socket.on("reconnect", function() {
      console.log("Reconnected");
      if(initiator) {
        socket.emit('reconnection', {name: session_name, secret: initiator, sid: null});
      } else if(joined) {
        socket.emit('reconnection', {name: session_name, secret: null, sid: joined});        
      }
    });    
    
    socket.on("connect_timeout", function() {
      console.log("Connect timeout");
    });    
    
    socket.on('session count', function(count) {
        if (initiator) setNumConnections(count);
    });

    socket.on('viewer view change response', function(new_view) {
        if (joined && !initiator) {
            glviewer.setView(new_view);
        }
    });

    socket.on('viewer state change response', function(new_state) {
        if (joined && !initiator) {
            if (window.viewerApp && typeof window.viewerApp.acceptSessionState === 'function') {
                var applied = window.viewerApp.acceptSessionState(new_state);
                applied.catch(window.viewerApp.showError);
                return applied;
            } else {
                glviewer.setInternalState(new_state);
            }
        }
    });

    var result_labels = []; // array of labels

    // Remove result labels
    var clearResultLabels = function() {
        if (result_labels) {
            if (result_labels.length) glviewer.removeLabels(result_labels);
            result_labels = [];
        }
    };

    $('#refreshresults').on('click', function() {
        if (!initiator) return;
        socket.emit('query fetch', {
            name : session_name,
        });
     });

    $('#askbutton').on('click', function() {
        if (!initiator) return;
        var val = $('#askbutton').text();
        if (val[0] == 'Q') { // start a query
            // tell server to query clients
            socket.emit('query start', {
                name : session_name,
                kind : "atoms"
            });
        } else if (val[0] == 'S') { // show results
            $('#askbutton').html('Clear Labels');
            $('#refreshresults').show();            
            socket.emit('query fetch', {
                name : session_name,
            });
        } else if (val[0] == 'C') { // clear results
            socket.emit('query end', {
                name : session_name,
            });
        }
    });

    var clientFinishQuery = function () {
        $('#selectmessage').hide();
        glviewer.setHoverable({},false);
        glviewer.setClickable({},false);
        clearQueryLabels();    
    };
    
    var clicked_labels = [];
    var current_hover_label = null;
    
    var clearQueryLabels = function() {
      var labels = clicked_labels.slice();
      if (current_hover_label) labels.push(current_hover_label);
      if (labels.length) glviewer.removeLabels(labels);
      clicked_labels = [];
      current_hover_label = null;
    };
    
    // Return information about what is currently clicked for the server
    var getClicked = function() {
        const ret = Array.from(clicked_labels, ({ stylespec }) => [stylespec.position.x, stylespec.position.y, stylespec.position.z]);
        return ret;
    };
    socket.on('query start response', function() {
        if (!initiator && !joined) return;
        if (initiator) {
            $('#askbutton').html('Show Results');
            $('#responseinfo').html('Responses: 0');
        } else {
            // Show message to select atoms, setup callbacks
            $('#selectmessage').show();
            
            //how to label an atom
            var atomLabel = function(atom) {
                var l = '';
                if(atom.resn) l = atom.resn+":"+atom.atom;
                else if(atom.atom) l = atom.atom;
                else l = atom.elem;
                return l;
            };
            
            glviewer.setHoverable({},true,function(atom,viewer,event,container) {
                if(!atom.hoverlabel && !atom.clicklabel) {                    
                    current_hover_label = atom.hoverlabel = viewer.addLabel(atomLabel(atom),{position: atom, backgroundColor: 0x800080, backgroundOpacity: 0.7, fontColor:'white'});
                }
            },
            function(atom) { 
                if(atom.hoverlabel) {
                 current_hover_label = null;
                 glviewer.removeLabel(atom.hoverlabel);
                 delete atom.hoverlabel;
                }
             }
            );            
            
            glviewer.setClickable({},true, function(atom, viewer) {
               if(atom.clicklabel) {
                   // If already clicked, deselect
                   let idx = clicked_labels.indexOf(atom.clicklabel);
                   if(idx >= 0) {
                       clicked_labels.splice(idx,1);
                   }
                   viewer.removeLabel(atom.clicklabel);
                   delete atom.clicklabel;
               } else { 
                   // Select
                   atom.clicklabel = viewer.addLabel(atomLabel(atom), {position: atom, backgroundColor: 0x800080, backgroundOpacity: 1.0, fontColor: 'white'});
                   clicked_labels.push(atom.clicklabel);
                   if(atom.hoverlabel) {
                       glviewer.removeLabel(atom.hoverlabel);
                       if(current_hover_label == atom.hoverlabel) current_hover_label = null;
                       delete atom.hoverlabel;
                   }                  
               }
               // Send clicked atoms to server
               socket.emit('query update', {
                   name : session_name,
                   selected : getClicked()
               });               
            });
            
            glviewer.render(); //calling render is required to register clickables - should probably fix that
        }
    });

    socket.on('query end response', function() {
        if (!initiator && !joined) return;
        if (initiator) {
            $('#askbutton').html('Query Atoms');
            $('#responseinfo').html('');
            $('#refreshresults').hide();

            clearResultLabels();
        } else {
            // Remove message, callbacks
            clientFinishQuery();
        }
    });

    // Update number of responses - result is list of {'cnt': N, 'position': (x,y,z)}
    socket.on('query update response', function(result) {
        if (initiator) {
            $('#responseinfo').html('Responses: ' + result);            
        }
    });
    
    socket.on('query fetch response', function(query_result) {
        if (!initiator) return;
        
        // First clear existing labels
        clearResultLabels();
        
        let max = 1;
        let i = 0;
      
        // Calculate max
        query_result.forEach((result) => {
            const val = parseInt(result[1]);
            if (val > max) max = val;
          });
      
        i = 0;
        result_labels = Array.from(query_result, function(item) {
          const pos = item[0];
          const p = { x: pos[0], y: pos[1], z: pos[2] };
          const val = parseInt(item[1]);
          return glviewer.addLabel(val, { position: p, backgroundOpacity: 0.5 + (val / max) * 0.5 });
        });
      });
      
    
};
