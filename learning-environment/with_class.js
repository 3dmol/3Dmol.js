$(document).ready(function () {
    var session = new $3Dmol.Session("localhost");
    session.viewer = glviewer;
    session.connect()
        .then(() => {
            // show only the connect button, make sure sidebar is back
            var resetUpperRight = function () {
                $('#sessionbutton').show();
                $('#sessionconnect').hide();
                $('#sessioncontrol').hide();
                $('#sessionmonitor').hide();
                $("#menu").show();
                $("#sidenav").show();
            };

            // open session connection dialog
            $('#sessionbutton').on('click', function () {
                $('#sessionbutton').hide();
                $('#sessionconnect').show();
            });

            var setNumConnections = function (count) {
                $('#userinfo').html('Users: ' + count);
            };
            setNumConnections(0);

            // close dialog and show sessions button
            $('#sessionconnectclose').on('click', function () {
                resetUpperRight();
            });

        })
        .catch((err) => console.log(err));

    $('#session_name_input').on('input', function () {
        var name = $('#session_name_input').val();
        $('#sessionconnectbutton').prop('disabled', false);
        session.checkName(name)
            .then(() => {
                $('#sessionconnectbutton').html("Create");
            })
            .catch((err) => {
                $('#sessionconnectbutton').html("Join");
            });
    });
    var connectSession = function () {
        var name = $('#session_name_input').val();
        console.log("ada")
        if ($('#sessionconnectbutton').text() == "Create") {
            session.create(name).then(() => {
                    console.log("session created successfully");
                    $('.sessionname').html(name);
                    // close the connection create pane and open the
                    // connection monitoring
                    $('#sessionbutton').hide();
                    $('#sessionconnect').hide();
                    $('#sessioncontrol').show();
                })
                .catch((err) => {
                    alert(err);
                });
        } else {
            session.join(name).then(() => {
                    $('.sessionname').html(name);
                    $('#sessionbutton').hide();
                    $('#sessionconnect').hide();
                    $("#menu").hide();
                    $("#sidenav").hide();
                    $('#sessionmonitor').show();

                })
                .catch((err) => {
                    alert(err);
                });
        }
    };
    // show only the connect button, make sure sidebar is back
    var resetUpperRight = function () {
        $('#sessionbutton').show();
        $('#sessionconnect').hide();
        $('#sessioncontrol').hide();
        $('#sessionmonitor').hide();
        $("#menu").show();
        $("#sidenav").show();
    };
    var deleteSession = function () {
        session.delete().then(() => {
            $('#createSession,#joinSession').prop('disabled', false);
            resetUpperRight();
            console.log("Session Deleted", session.name);
        });

    };

    var leaveSession = function () {
        session.leave().then(() => {
            console.log("Session Left", session.name);
            resetUpperRight();
        });
    };
    
    $('#sessiondestroy').on('click', deleteSession);
    $('#sessionleave').on('click', leaveSession);
    $('#sessionconnectbutton').on('click', connectSession);
    $('#session_name_input').on('keyup', function (e) {
        if (e.keyCode == 13) { // enter key
            connectSession();
        }
    });

});