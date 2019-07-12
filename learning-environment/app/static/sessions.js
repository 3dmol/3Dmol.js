    
var socket = null;
var session_name = null;

var joinSession = function(name) {
    session_name = name;
    socket.emit( 'join session event', {
        name : name
    });
};



//setup socket handlers and session related event handlers
var initSessions = function () {
	var initiator = false;
	var joined = false;
	var old_url = null;
	var old_view = null;
	var fps = 20;
	socket = io.connect(window.location.hostname+":5000");
    socket.on('connect', function() {
    	socket.send('User has connected!');
    });
    
    //show only the connect button, make sure sidebar is back
    var resetUpperRight = function() {
        $('#sessionbutton').show();
        $('#sessionconnect').hide();
        $('#sessioncontrol').hide();
        $('#sessionmonitor').hide();
        $("#menu").show();
        $("#sidenav").show();
    };
    
    //open session connection dialog
    $('#sessionbutton').on('click', function() {
        $('#sessionbutton').hide();
        $('#sessionconnect').show();
    });
    
    //close dialog and show sessions button
    $('#sessionconnectclose').on('click',function() {
        resetUpperRight();
    });
    
    //as the user types a session name, check for validity
    $('#session_name_input').on('input', function() {
        var name = $('#session_name_input').val();
        if(name.length > 0) {
            socket.emit('check session name event', {
                name : name
            });     
        } else {
            $('#sessionconnectbutton').prop('disabled',true);
        }
    });
    
    //update button based on existance of name
    socket.on('check session name response', function(msg){
        $('#sessionconnectbutton').prop('disabled',false);
        if(msg == 'exists'){
            $('#sessionconnectbutton').html("Join");
        } else {
            $('#sessionconnectbutton').html("Create");
        }
    });
    
    //handle connecting to a session
    var connectSession = function() {
        session_name = $('#session_name_input').val();
        if($('#sessionconnectbutton').text() == "Create") {
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
    
    $('#sessionconnectbutton').on('click',connectSession);
    $('#session_name_input').on('keyup',function(e){
       if(e.keyCode == 13) { //enter key
           connectSession();
       } 
    });
    
    //register change callbacks
    var viewUpdateCallback = function(new_view){
        socket.emit('viewer view change event', {
            name : session_name,
            view : new_view
        });         
    };

    var stateUpdateCallback = function(new_state){
            socket.emit('viewer state change event', {
                name : session_name,
                state : new_state
            });         
    };
        
    socket.on('create session response', function(msg){
    	if(msg == 1){
    		console.log("session created successfully")
    		
    		//setup callbacks
    		glviewer.setViewChangeCallback(viewUpdateCallback);
    		glviewer.setStateChangeCallback(stateUpdateCallback);   
    		initiator = true;
    		
    		//close the connection create pane and open the connection monitoring
            $('#sessionbutton').hide();
            $('#sessionconnect').hide();
            $('#sessioncontrol').show();
    	}
    	else
    		alert("Session name was already taken/ could not be created. Try Again")
    });
    
    socket.on('join session response', function(msg){
    	if(msg==1){
    	    joined = true;
            //close the connection create pane and open the connection monitoring
    	    //disable the sidebar

            $('#sessionbutton').hide();
            $('#sessionconnect').hide();
            $("#menu").hide();
            $("#sidenav").hide();
            $('#sessionmonitor').show();

    	}
    	else
    		alert("Session Doesn't Exist")
    });
   
    
    var deleteSession = function() {
        if(session_name) socket.emit('delete session event', { name : session_name}); 
        session_name = null;
    };
    
    var leaveSession = function() {
        if(session_name) {
            socket.emit( 'leave session event', {
                name : session_name,
            });
        }
        session_name = null;
    };
    
    $('#sessiondestroy').on('click',deleteSession);
    $('#sessionleave').on('click',leaveSession);
    
    //make sure to unregister a session if user closes window
    $(window).on('beforeunload',function() { 
        if(initiator) deleteSession();
        else leaveSession();
    });    
    
    var sessionEnded = function() {
        joined = false;
        resetUpperRight();
    };

    socket.on('delete session response', function(){
        $('#createSession,#joinSession').prop('disabled', false);
        initiator = false;
        //remove callbacks
        glviewer.setViewChangeCallback(null);
        glviewer.setStateChangeCallback(null);
        resetUpperRight();

    });
    
    socket.on('leave session response', sessionEnded);
    socket.on('disconnect', sessionEnded);
       
    socket.on('error: restart connection', function(){
    	location.reload();
    })
    socket.on('session count', function(count){
    	$('#sessioncontrolinfo').html('Active Connections: ' + count)
    });
    
    
    socket.on('viewer view change response', function(new_view){
    	if(!initiator){
    		glviewer.setView(new_view);
    	}
    });
    
    socket.on('viewer state change response', function(new_state){
    		if(!initiator){
    			glviewer.setInternalState(new_state);
    		}
    	});
    		
    
 };
