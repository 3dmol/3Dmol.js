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

    // Webserver needs to have appropriate rules to forward to flask
    // https://stackoverflow.com/questions/36472920/apache-proxy-configuration-for-socket-io-project-not-in-root
    if(window.location.hostname == 'localhost') {
      // For debugging on localhost go straight to port to avoid having to setup webserver
      socket = io.connect(window.location.hostname+":5000");
    } else {
      socket = io.connect(window.location.hostname);
    }
    socket.on('connect', function() {
        socket.send('User has connected!');
    });

    // Show only the connect button, make sure sidebar is back
    var resetUpperRight = function() {
        $('#sessionbutton').show();
        $('#sessionconnect').hide();
        $('#sessioncontrol').hide();
        $('#sessionmonitor').hide();
        $("#menu").show();
        $("#sidenav").show();
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
        resetUpperRight();
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
        $('#sessionconnectbutton').prop('disabled', false);
        if (msg == 'exists') {
            $('#sessionconnectbutton').html("Join");
        } else {
            $('#sessionconnectbutton').html("Create");
        }
    });

    // Handle connecting to a session
    var connectSession = function() {
        session_name = $('#session_name_input').val();
        if ($('#sessionconnectbutton').text() == "Create") {
            socket.emit('create session event', {
                name : session_name,
                state : glviewer.getInternalState(),
                view : glviewer.getView()
            });
        } else {
            socket.emit('join session event', {
                name : session_name,
            });
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
        socket.emit('viewer view change event', {
            name : session_name,
            view : new_view
        });
    };

    var stateUpdateCallback = function(new_state) {
        socket.emit('viewer state change event', {
            name : session_name,
            state : new_state
        });
    };

    socket.on('create session response',function(msg) {
        if (msg) {
            console.log("session created successfully")
            $('.sessionname').html(session_name);

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
        }
    });

    socket.on('join session response', function(msg) {
        if (msg != "0") {
            joined = msg;
            // Close the connection create pane and open the connection monitoring
            $('.sessionname').html(session_name);
            $('#sessionbutton').hide();
            $('#sessionconnect').hide();
            $("#menu").hide();
            
            // Disable the sidebar
            $("#sidenav").hide();
            $('#sessionmonitor').show();

        } else
            alert("Session Doesn't Exist")
    });

    var deleteSession = function() {
        if (session_name)
            socket.emit('delete session event', {
                name : session_name
            });
        clearResultLabels();
        session_name = null;
    };

    var leaveSession = function() {
        if (session_name) {
            socket.emit('leave session event', {
                name : session_name,
            });
        }
        session_name = null;
    };

    $('#sessiondestroy').on('click', deleteSession);
    $('#sessionleave').on('click', leaveSession);

    // Make sure to unregister a session if user closes window
    $(window).on('beforeunload', function() {
        if (initiator)
            deleteSession();
        else
            leaveSession();
    });

    var sessionEnded = function(reason) {
        joined = false;
        resetUpperRight();
        clientFinishQuery();
    };

    socket.on('delete session response', function() {
        $('#createSession,#joinSession').prop('disabled', false);
        initiator = 0;
        // Remove callbacks
        glviewer.setViewChangeCallback(null);
        glviewer.setStateChangeCallback(null);
        resetUpperRight();

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
    
    socket.on('session count', setNumConnections);

    socket.on('viewer view change response', function(new_view) {
        if (!initiator) {
            glviewer.setView(new_view);
        }
    });

    socket.on('viewer state change response', function(new_state) {
        if (!initiator) {
            glviewer.setInternalState(new_state);
        }
    });

    var result_labels = []; // array of labels

    // Remove result labels
    var clearResultLabels = function() {
        if (result_labels) {
            for (let i = 0; i < result_labels.length; i++) {
                glviewer.removeLabel(result_labels[i]);
            }
            result_labels = [];
        }
    };

    $('#refreshresults').on('click', function() {
        socket.emit('query fetch', {
            name : session_name,
        });
     });

    $('#askbutton').on('click', function() {
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
      for(let i = 0; i < clicked_labels.length; i++) {
          glviewer.removeLabel(clicked_labels[i]);
      }
      clicked_labels = [];
      if(current_hover_label) {
          glviewer.removeLabel(current_hover_label);
          current_hover_label = null;
      }
    };
    
    // Return information about what is currently clicked for the server
    var getClicked = function() {
        const ret = Array.from(clicked_labels, ({ stylespec }) => [stylespec.position.x, stylespec.position.y, stylespec.position.z]);
        return ret;
    };
    socket.on('query start response', function() {
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