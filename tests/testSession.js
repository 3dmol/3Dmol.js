var session = new $3Dmol.Session("localhost");
session.connect().then(function (response) { //Client is connected to the server
    console.log("Success!", response);
    var name = "Session#1";
    session.checkName(name)
    .then(() => { //Checks whether name exists 
        session.create(name)
         
            .then(() => { // If it does not exist then session is created
                console.log("Session Created");
                return session.delete();
            })
            .then((func) => {
                console.log(func);
                func().then(() => { //Once the session is created, it is deleted
                    console.log("Session Deleted");
                });
            })
            .catch(function (error) {
                console.log("Failed!", error);
            });
    }, (error) => {
        session.join(name)
            .then(() => { // Else the session is joined
                console.log("Session Joined");
                return session.leave();
            })
            .then((func) => {
                func().then(() => { //Once the session is joined, it is left
                    console.log("Session Left");
                });
            })
            .catch(function (error) {
                console.log("Failed!", error);
            });
    });

}).catch(function (error) {
    console.log("Failed!", error);
});
