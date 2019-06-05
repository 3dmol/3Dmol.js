
$(document).ready(function() {
	checkCookie();
	var socket = io.connect('http://localhost:5000');
	socket.on('connect', function() {
		socket.send('User has connected!');
	});

	function setCookie(exdays){
		var d = new Date();
		d.setTime(d.getTime() + (exdays*24*60*60*1000)) //Cookie stored for exDays
		var expires = "expires="+ d.toUTCString();
		var session_name = $('#session_name1').val();
		var session_password = $('#session_password1').val();
		document.cookie = "Session_Name=" + session_name +";" + expires +";path=/";
		document.cookie = "Session_Password=" + session_password + ";" + expires +";path=/";
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
			console.log("Session Name: " + cookie_session_name);
		}
		var cookie_session_password = getCookie("Session_Password");
		if(cookie_session_password != ""){
			$('#session_password1').val(cookie_session_password);
			console.log("Session Password:" + cookie_session_password);
		}
	} 

	$('input#session_name1').on('input', function(){
		socket.emit( 'check session name event', {
			session_name : $( 'input#session_name1' ).val()
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
			session_name : $( 'input#session_name1' ).val(),
			session_password : $( 'input#session_password1' ).val()
		});
		
	});
	socket.on('create session response', function(msg){
		if(msg == 1){
			console.log("session created successfully")
			setCookie(1);
			$('#session_list1').hide();
			$('#session_list3').show();
			$('#createSession,#joinSession').prop('disabled', true);
		}
		else
			alert("Session name was already taken/ could not be created. Try Again")
	});
	$()
	$( '#join_session_form' ).on( 'submit', function( e ) {
		e.preventDefault();
		
		socket.emit( 'join session event', {
			session_name : $( 'input#session_name2' ).val(),
			session_password : $( 'input#session_password2' ).val()
		});
		
	});
	$( '#delete_session_form' ).on( 'submit', function( e ) {
		e.preventDefault();
		
		socket.emit( 'delete session event', {
			session_name : $( 'input#session_name1' ).val(),
			session_password : $( 'input#session_password1' ).val()
		});
		
	});
});