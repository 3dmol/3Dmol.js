
$(document).ready(function() {

	var socket = io.connect('http://localhost:5000');
	socket.on('connect', function() {
		socket.send('User has connected!');
	});
	$('input#session_name1').on('input', function(){
		socket.emit( 'check session name event', {
			session_name : $( 'input#session_name1' ).val()
		})
	})	
	socket.on('check session name response', function(msg){
		console.log(msg);
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

	  $( '#join_session_form' ).on( 'submit', function( e ) {
		e.preventDefault();
		
		socket.emit( 'join session event', {
			session_name : $( 'input#session_name2' ).val(),
			session_password : $( 'input#session_password2' ).val()
		});
		
	  });

});