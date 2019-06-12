
$(document).ready(function() {
	checkCookie();
	var initiator = false;
	var joined = false;
	var old_url = null;
	var old_view = null;
	var fps = 20;
	var socket = io.connect('http://localhost:5000');
	socket.on('connect', function() {
		socket.send('User has connected!');
	});
	
	function setCookie(exdays){
		var d = new Date();
		d.setTime(d.getTime() + (exdays*24*60*60*1000)) //Cookie stored for exDays
		var expires = "expires="+ d.toUTCString();
		var name = $('#session_name1').val();
		var password = $('#session_password1').val();
		document.cookie = "Session_Name=" + name +";" + expires +";path=/";
		document.cookie = "Session_Password=" + password + ";" + expires +";path=/";
	}
	function getCookie(cname) {
		var name = cname + "=";
		var decodedCookie = decodeURIComponent(document.cookie);
		var ca = decodedCookie.split(';');
		for(var i = 0; i <ca.length; i++) {
		  var c = ca[i];
		  while (c.charAt(0) == ' ') {
			c = c.substring(1);
		  }
		  if (c.indexOf(name) == 0) {
			return c.substring(name.length, c.length);
		  }
		}
		return "";
	  }
	function checkCookie(){
		var cookie_session_name = getCookie("Session_Name");
		if(cookie_session_name != ""){
			$('#session_name1').val(cookie_session_name);
			console.log("Cookie Session Name: " + cookie_session_name);
		}
		var cookie_session_password = getCookie("Session_Password");
		if(cookie_session_password != ""){
			$('#session_password1').val(cookie_session_password);
			console.log("Cookie Session Password:" + cookie_session_password);
		}
	} 

	$('input#session_name1').on('input', function(){
		socket.emit( 'check session name event', {
			name : $( 'input#session_name1' ).val()
		})
	})	
	socket.on('check session name response', function(msg){
		if(msg == 1)
			$('#check_name').html('&#10004;'); 	// Tick Mark
		else
			$('#check_name').html('&#10060;');	// Cross Mark
	});
	$( '#create_session_form' ).on( 'submit', function( e ) {
		e.preventDefault();
		
		socket.emit( 'create session event', {
			name : $( 'input#session_name1' ).val(),
			password : $( 'input#session_password1' ).val()
		});
		
		
	});
	socket.on('create session response', function(msg){
		if(msg == 1){
			console.log("session created successfully")
			setCookie(1);
			$('#session_list1').hide();
			$('#session_list3').show();
			$('#createSession,#joinSession').prop('disabled', true);
			initiator = true;
			socket.emit( 'URL state change event', {
				name : $( 'input#session_name1' ).val(),
				url : window.location.href.substring(window.location.href.indexOf("?")+1)
		});
		}
		else
			alert("Session name was already taken/ could not be created. Try Again")
	});
	
	$( '#join_session_form' ).on( 'submit', function( e ) {
		e.preventDefault();
		
		socket.emit( 'join session event', {
			name : $( 'input#session_name2' ).val(),
			password : $( 'input#session_password2' ).val()
		});
			
	});
	socket.on('join session response', function(msg){
		if(msg==1){
		$('#session_list2').hide();
		$('#session_list4').show();
		$('#createSession,#joinSession,#addStyle,#addSurface,#addLabelRes').prop('disabled', true);
		joined = true;
		// window.location.href.substring(window.location.href.indexOf("?"))='';
		}
		else
			alert("Session Doesn't Exist")
	});
	$( '#delete_session_form' ).on( 'submit', function( e ) {
		e.preventDefault();
		
		socket.emit( 'delete session event', {
			name : $( 'input#session_name1' ).val(),
			password : $( 'input#session_password1' ).val()
		});
		
	});
	socket.on('delete session response', function(){
		$('#session_list3').hide();
		$('#createSession,#joinSession').prop('disabled', false);
		initiator = false;
	});
	$( '#leave_session_form' ).on( 'submit', function( e ) {
		e.preventDefault();
		
		socket.emit( 'leave session event', {
			name : $( 'input#session_name1' ).val(),
			password : $( 'input#session_password1' ).val()
		});
		
	});
	socket.on('leave session response', function(){
		$('#session_list4').hide();
		$('#createSession,#joinSession,#addStyle,#addSurface,#addLabelRes').prop('disabled', false);
		joined = false;

	});
	socket.on('error: restart connection', function(){
		location.reload();
	})
	socket.on('session count', function(count){
		$('#people_joined').html('People Joined: ' + count)
	});
	window.setInterval(function (){
		if(initiator==true){
			var new_url=window.location.href.substring(window.location.href.indexOf("?")+1);
			if(old_url != new_url){
				old_url = new_url;
				console.log(new_url);
				socket.emit( 'URL state change event', {
					name : $( 'input#session_name1' ).val(),
					url : new_url
			});
		}
		}
	}, 500);

	socket.on('URL state change response', function(url){
		if(initiator == false){
		console.log("URL Response");
		setURL(url)
		var query = urlToQuery(window.location.search.substring(1));
		buildHTMLTree(query);
		render(true);
		// glviewer.translate(width/2,0,0,false);		
		}
	});
	window.setInterval(function(){
		if(initiator == true){
			var new_view = glviewer.getView();
				socket.emit('viewer state change event', {
					name : $( 'input#session_name1' ).val(),
					viewer : new_view
				});			
	}
	}, 1000/fps);

	socket.on('viewer state change response', function(new_view){
		if(initiator == false){
			glviewer.setView(new_view);
			render(true);
			// glviewer.render();
			// glviewer.translate(0,0,0,false);
		}
	})
});