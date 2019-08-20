var session = new $3Dmol.Session("localhost");

session.viewer = $3Dmol.createViewer(
    "gldiv", //id of div to create canvas in
    {
        defaultcolors: $3Dmol.rasmolElementColors
    }
);
var initiator;

var fun = function (viewer) {
    console.log('fun!');
    viewer.clear();
    viewer.addSphere({
        center: {
            x: 0,
            y: 0,
            z: 0
        },
        radius: 10.0,
        color: 'red'
    });
    viewer.addLine({
        color: 'blue',
        start: {
            x: 0,
            y: 0,
            z: 0
        },
        end: {
            x: 0,
            y: 5,
            z: 0
        }
    });
    viewer.addLine({
        color: 'green',
        start: {
            x: 0,
            y: 0,
            z: 0
        },
        end: {
            x: 0,
            y: 0,
            z: 5
        }
    });
    viewer.addBox({
        center: {
            x: 0,
            y: 0,
            z: 0
        },
        dimensions: {
            w: 3,
            h: 4,
            d: 2
        },
        color: 'magenta'
    });
    viewer.setStyle({}, {
        line: {}
    });
    viewer.render();
    viewer.zoomTo();
};

session.connect()
    .then(function (response) { //Client is connected to the server
        console.log("Success! Connected");
        var name = "Session#1";
        session.checkName(name)
            .then(() => { //Checks whether name exists 
                console.log("Name is not taken");
                session.create(name)
                    .then(() => {
                        console.log("Session Created");
                        initiator = true;
                        //Chain .then with functions of your choice to be called after session creation
                    })
                    .then(() => {
                        console.log("Write the next function here");
                        fun(session.viewer);
                    })
                    .then(() => {
                        console.log(session.name);
                        // session.delete().then(console.log("Session Deleted"));

                    })
                    .catch((error) => {
                        console.log("Failed!", error);
                    });
            })
            .catch((error) => {
                console.log("Name is taken");
                session.join(name)
                    .then(() => {
                        console.log("Session Joined", session.name);
                        //Chain .then with functions of your choice to be called after session joining.
                    })
                    .then(() => {
                        console.log("Write the next function here");
                        // fun(session.viewer);
                    })
                    .then(() => {
                        // session.leave().then(console.log("Session Left"));   
                    })
                    .catch((error) => {
                        console.log("Failed!", error);
                    });
            });


    }).catch(function (error) {
        console.log("Failed! Connection", error);
    });
$(window).on('beforeunload', function () {
    if (initiator)
        session.delete().then(console.log("session deleted"));
    else
        session.leave().then(console.log("session left"));
});