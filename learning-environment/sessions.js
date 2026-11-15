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
    // Keep one connection and one set of role-owning handlers for this window.
    if (socket) return socket;
    var initiator = false;
    var joined = false;
    var pendingJoin = false;
    var pendingCreate = false;
    var queryActive = false;
    var revision = 0, lastSentState = null, pendingSnapshot = null, receivedState = null;
    var stateTimer = null, stateDirty = false, forceFull = false, viewTimer = null, pendingView = null, lastView = '';
    var unsubscribeState = null, unsubscribeScene = null, reconnecting = false;
    var snapshot = function() {
        var app = window.viewerApp;
        var value = app && typeof app.getSessionState === 'function' ? app.getSessionState()
            : Object.assign({}, glviewer.getInternalState({ compact: true }), { frame: glviewer.getFrame() });
        if (value.labels) value = Object.assign({}, value, { labels: value.labels.filter(label => !label.style?._sessionTransient) });
        return value;
    };
    var scheduleState = function(full) {
        if (!initiator || !session_name) return;
        stateDirty = true; forceFull = forceFull || !!full;
        if (stateTimer === null) stateTimer = setTimeout(flushState, 100);
    };
    var flushState = function() {
        stateTimer = null;
        if (!initiator || !session_name || !stateDirty || pendingSnapshot || socket.connected === false) return;
        var next;
        try { next = JSON.parse(JSON.stringify(snapshot())); }
        catch (error) { window.viewerApp?.showError(error); return; }
        var patch = {};
        Object.keys(next).forEach(function(key) {
            if (!lastSentState || JSON.stringify(next[key]) !== JSON.stringify(lastSentState[key])) patch[key] = next[key];
        });
        stateDirty = false;
        if (!forceFull && lastSentState && !Object.keys(patch).length) return;
        var message = { name: session_name, baseRevision: revision, revision: revision + 1 };
        if (forceFull || !lastSentState) message.state = next;
        else message.patch = patch;
        forceFull = false; pendingSnapshot = next;
        socket.emit('viewer state change event', message);
    };
    var acknowledgeState = function(message) {
        if (!initiator || !Number.isInteger(message.revision) || message.revision < revision) return;
        revision = message.revision;
        if (pendingSnapshot) { lastSentState = pendingSnapshot; pendingSnapshot = null; }
        if (stateDirty) scheduleState(false);
    };
    var reconnect = function() {
        if (reconnecting || (!initiator && !joined) || !session_name) return;
        reconnecting = true;
        socket.emit('reconnection', { name: session_name, secret: initiator || null, sid: joined || null });
    };

    joinSession = function(name) {
        if (initiator || joined || pendingJoin || pendingCreate) return;
        session_name = name;
        pendingJoin = true;
        socket.emit('join session event', { name: name });
    };

    // All classroom dialogs share a floating position until the session closes.
    var dialogs = Array.from(window.document?.querySelectorAll('.rightdialog') || []);
    var dialogPosition = null, activeDialog = null, dialogDrag = null;
    var placeDialog = function(dialog, left, top) {
        var bounds = dialog.getBoundingClientRect();
        dialogPosition = {
            left: Math.max(12, Math.min(left, window.innerWidth - bounds.width - 12)),
            top: Math.max(12, Math.min(top, window.innerHeight - bounds.height - 12))
        };
        dialog.style.left = dialogPosition.left + 'px';
        dialog.style.top = dialogPosition.top + 'px';
        dialog.style.right = 'auto';
    };
    var endDialogDrag = function() {
        if (!dialogDrag) return;
        var previous = dialogDrag; dialogDrag = null;
        previous.dialog.classList.remove('is-dragging');
        window.document.removeEventListener('pointermove', moveDialog);
        window.document.removeEventListener('pointerup', finishDialogDrag);
        window.document.removeEventListener('pointercancel', finishDialogDrag);
        if (previous.handle.hasPointerCapture?.(previous.pointerId)) {
            previous.handle.releasePointerCapture(previous.pointerId);
        }
    };
    var moveDialog = function(event) {
        if (!dialogDrag || event.pointerId !== dialogDrag.pointerId) return;
        event.preventDefault(); event.stopPropagation();
        placeDialog(dialogDrag.dialog, event.clientX - dialogDrag.x, event.clientY - dialogDrag.y);
    };
    var finishDialogDrag = function(event) {
        if (!dialogDrag || event.pointerId !== dialogDrag.pointerId) return;
        if (event.type === 'lostpointercapture' && event.target !== dialogDrag.handle) return;
        event.stopPropagation(); endDialogDrag();
    };
    dialogs.forEach(function(dialog) {
        var handle = dialog.querySelector('[data-session-drag]');
        if (!handle) return;
        // Keep controls and dragging from reaching molecular viewer gestures.
        ['pointerdown', 'mousedown', 'touchstart', 'wheel', 'dblclick'].forEach(function(type) {
            dialog.addEventListener(type, function(event) { event.stopPropagation(); });
        });
        handle.addEventListener('pointerdown', function(event) {
            if (event.button !== 0 || event.isPrimary === false) return;
            endDialogDrag();
            var bounds = dialog.getBoundingClientRect();
            event.preventDefault(); event.stopPropagation();
            handle.focus({ preventScroll: true });
            dialogDrag = { dialog: dialog, handle: handle, pointerId: event.pointerId,
                x: event.clientX - bounds.left, y: event.clientY - bounds.top };
            dialog.classList.add('is-dragging');
            window.document.addEventListener('pointermove', moveDialog);
            window.document.addEventListener('pointerup', finishDialogDrag);
            window.document.addEventListener('pointercancel', finishDialogDrag);
            // Document listeners still handle drags when capture is unavailable.
            try { handle.setPointerCapture?.(event.pointerId); } catch (_) {}
        });
        handle.addEventListener('lostpointercapture', finishDialogDrag);
        handle.addEventListener('keydown', function(event) {
            var direction = { ArrowLeft: [-1, 0], ArrowRight: [1, 0], ArrowUp: [0, -1], ArrowDown: [0, 1] }[event.key];
            if (!direction || event.altKey || event.ctrlKey || event.metaKey) return;
            event.preventDefault(); event.stopPropagation();
            var bounds = dialog.getBoundingClientRect(), step = event.shiftKey ? 1 : 10;
            placeDialog(dialog, bounds.left + direction[0] * step, bounds.top + direction[1] * step);
        });
    });
    window.addEventListener?.('resize', function() {
        if (activeDialog && dialogPosition) placeDialog(activeDialog, dialogPosition.left, dialogPosition.top);
    });
    window.addEventListener?.('blur', endDialogDrag);
    var showDialog = function(selector) {
        endDialogDrag();
        $('#sessionbutton').hide();
        $('#sessionconnect,#sessioncontrol,#sessionmonitor').hide();
        $(selector).show();
        activeDialog = window.document?.querySelector(selector);
        if (activeDialog && dialogPosition) placeDialog(activeDialog, dialogPosition.left, dialogPosition.top);
    };
    var resetDialogs = function() {
        endDialogDrag();
        activeDialog = dialogPosition = null;
        dialogs.forEach(function(dialog) {
            dialog.style.left = dialog.style.top = dialog.style.right = '';
        });
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

    // Use the page's server (including any reverse proxy) unless explicitly overridden.
    var sessionServer = window.viewerSessionServer || window.location.origin;
    socket = io.connect(sessionServer);
    socket.on('connect', function() {
        socket.send('User has connected!');
        reconnect();
    });

    // Show only the connect button, make sure sidebar is back
    var resetUpperRight = function() {
        resetDialogs();
        $('#sessionbutton').show();
        $('#sessionconnect').hide();
        $('#sessioncontrol').hide();
        $('#sessionmonitor').hide();
        showViewerControls(true);
        sessionChanged(null, false);
    };

    // Open session connection dialog
    $('#sessionbutton').on('click', function() {
        showDialog('#sessionconnect');
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
        if (initiator || joined || pendingJoin || pendingCreate) return;
        session_name = $('#session_name_input').val();
        if ($('#sessionconnectbutton').text() == "Create") {
            pendingCreate = true;
            lastSentState = JSON.parse(JSON.stringify(snapshot()));
            socket.emit('create session event', {
                name : session_name,
                state : lastSentState,
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
        var key = JSON.stringify(new_view);
        if (key === lastView && viewTimer === null) return;
        pendingView = new_view.slice();
        if (viewTimer === null) viewTimer = setTimeout(function() {
            viewTimer = null;
            if (!initiator || !session_name || socket.connected === false) return;
            var key = JSON.stringify(pendingView);
            if (key !== lastView) {
                lastView = key;
                socket.emit('viewer view change event', { name: session_name, view: pendingView });
            }
        }, 33);
    };
    var stateUpdateCallback = function() { scheduleState(false); };

    socket.on('viewer state revision', acknowledgeState);
    socket.on('viewer state resync', function(message) {
        if (!initiator) return;
        revision = message.revision; pendingSnapshot = null;
        scheduleState(true);
    });
    socket.on('reconnection response', function(message) {
        reconnecting = false;
        revision = message.revision;
        if (message.controller && initiator) {
            pendingSnapshot = null;
            scheduleState(true);
            lastView = ''; viewUpdateCallback(glviewer.getView());
        } else if (joined) joined = message.sid;
    });

    socket.on('create session response',function(msg) {
        if (!pendingCreate || !session_name) return;
        pendingCreate = false;
        if (msg) {
            console.log("session created successfully")
            $('.sessionname').text(session_name);

            // Setup callbacks
            glviewer.setViewChangeCallback(viewUpdateCallback);
            initiator = parseInt(msg);
            revision = 0; pendingSnapshot = null; stateDirty = false;
            lastView = JSON.stringify(glviewer.getView());
            if (glviewer.subscribeStateChanges) unsubscribeScene = glviewer.subscribeStateChanges(stateUpdateCallback);
            else glviewer.setStateChangeCallback(stateUpdateCallback);
            if (window.viewerApp?.subscribeSessionChanges) unsubscribeState = window.viewerApp.subscribeSessionChanges(stateUpdateCallback);
            // Include edits made while the create request was in flight.
            scheduleState(false);
    
            // Close the connection create pane and open the connection monitoring
            showDialog('#sessioncontrol');
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
            joined = msg; revision = 0; receivedState = null;
            // Close the connection create pane and open the connection monitoring
            $('.sessionname').text(session_name);
            // Disable the sidebar
            sessionChanged(session_name, true);
            showViewerControls(false);
            showDialog('#sessionmonitor');

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
        if (stateTimer !== null) clearTimeout(stateTimer);
        if (viewTimer !== null) clearTimeout(viewTimer);
        stateTimer = viewTimer = null; stateDirty = forceFull = false;
        pendingSnapshot = lastSentState = receivedState = pendingView = null; revision = 0; lastView = ''; reconnecting = false;
        unsubscribeState?.(); unsubscribeScene?.(); unsubscribeState = unsubscribeScene = null;
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
      reconnecting = false;
      // Do NOT end session since we will attempt to reconnect
    });

    socket.on('error: restart connection', function() {
        console.log("error: restart connection");
    });
    
    socket.on('error', (error) => {
        console.log('error: '+error);
    });
    
    // Older Socket.IO clients emit reconnect on the socket; current clients
    // reach the same path through connect. The guard avoids duplicate requests.
    socket.on('reconnect', reconnect);

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

    socket.on('viewer state change response', function(message) {
        if (initiator) { acknowledgeState(message); return; }
        if (!joined) return;
        var sections;
        if (message.patch) {
            if (Number.isInteger(message.revision) && message.revision <= revision) return;
            if (!receivedState || message.baseRevision !== revision) {
                socket.emit('request state', { name: session_name });
                return;
            }
            sections = Object.keys(message.patch);
            receivedState = Object.assign({}, receivedState, message.patch);
            revision = message.revision;
        } else {
            if (Number.isInteger(message.revision) && message.revision < revision) return;
            receivedState = message.state || message;
            revision = Number.isInteger(message.revision) ? message.revision : 0;
        }
        if (window.viewerApp && typeof window.viewerApp.acceptSessionState === 'function') {
            var applied = window.viewerApp.acceptSessionState(receivedState, sections).then(restoreQueryInteractions);
            applied.catch(window.viewerApp.showError);
            return applied;
        }
        glviewer.setInternalState(sections ? message.patch : receivedState, { partial: !!sections });
        restoreQueryInteractions();
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
            $('#refreshresults').css('display', 'inline-flex');
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
        queryActive = false;
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
    var clientStartQuery = function() {
        queryActive = true;
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
            if (!joined || !queryActive) return;
            if(!atom.hoverlabel && !atom.clicklabel) {
                current_hover_label = atom.hoverlabel = viewer.addLabel(atomLabel(atom),{_sessionTransient: true, position: atom, backgroundColor: 0x800080, backgroundOpacity: 0.7, fontColor:'white'});
            }
        },
        function(atom) {
            if (!joined || !queryActive) return;
            if(atom.hoverlabel) {
             current_hover_label = null;
             glviewer.removeLabel(atom.hoverlabel);
             delete atom.hoverlabel;
            }
         }
        );

        glviewer.setClickable({},true, function(atom, viewer) {
           if (!joined || !queryActive) return;
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
               atom.clicklabel = viewer.addLabel(atomLabel(atom), {_sessionTransient: true, position: atom, backgroundColor: 0x800080, backgroundOpacity: 1.0, fontColor: 'white'});
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
    };

    // Instructor picking flags are scene metadata, not permission for students
    // to interact. Preserve the student's current query across scene updates.
    var restoreQueryInteractions = function() {
        if (!joined || initiator) return;
        if (queryActive) clientStartQuery();
        else {
            glviewer.setHoverable({}, false);
            glviewer.setClickable({}, false);
            glviewer.render();
        }
    };

    socket.on('query start response', function() {
        if (!initiator && !joined) return;
        if (initiator) {
            $('#askbutton').html('Show Results');
            $('#responseinfo').html('Responses: 0');
        } else clientStartQuery();
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
          return glviewer.addLabel(val, { _sessionTransient: true, position: p, backgroundOpacity: 0.5 + (val / max) * 0.5 });
        });
      });
      
    return socket;
};
