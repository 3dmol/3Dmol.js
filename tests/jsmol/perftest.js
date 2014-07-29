//Test rendering performance for different sized pdb's
jmol_isReady = function(applet) {
    start();
};
//var glviewer = null;
var Info = {
    width: '100%',
    height: '100%',
    debug: false,
    color: "0xFFFFFF",
    use: "HTML5",   // JAVA HTML5 WEBGL are all options
    j2sPath: "./jsmol/jsmol/j2s", // this needs to point to where the j2s directory is.
    jarPath: "./jsmol/jsmol/java",// this needs to point to where the java directory is.
    jarFile: "JmolAppletSigned.jar",
    isSigned: true,
    script: "load test_structs/3M8L.pdb",
    readyFunction: jmol_isReady,
    disableJ2SLoadMonitor: true,
    disableInitialConsole: true,
    allowJavaScript: true,
    //defaultModel: "$dopamine",
    console: "none" // default will be jmolApplet0_infodiv, but you can designate another div here or "none"
};

$(document).ready(function() {
    stop();
    $("#gldiv").html(Jmol.getAppletHtml("viewer", Info));
    start();

}); 

var profile = QUnit.urlParams.profilecheck;
var resultXML = null;
var resultStr = "";

//QUnit-reporter hook to output test results in XML format
QUnit.jUnitReport = function(data) {

    resultXML = $.parseXML(data.xml);

    //Wrap XML result in JQuery object; parse and setup output string
    $result = $(resultXML);

    resultStr = "JSmol Performance Tests\n";
    var runTime = $result.find("testsuites").attr("time");
    var runDate = $result.find("testsuites").attr("timestamp");
    resultStr += "Total Test Time: " + runTime + " s\n";
    resultStr += "Date: " + runDate + "\n\n";

    var test = $result.find("testsuite").first();
    var moduleName = test.attr("name");
    var moduleTime = test.attr("time");
    resultStr += "\n" + moduleName;
    //alert(moduleName);
    test.find("testcase").each(function() {
        var testName = $(this).attr("name");
        var testTime = $(this).attr("time");
        resultStr += "\n\t" + testName + ":   " + testTime + " s";
        //alert(testName);
    });
    resultStr += "\n\tTotal:         " + moduleTime + " s\n";

    //Set up a link to download test results
    $("#qunit-testresult").append("<br><a id='download'>Download</a>");
    var url = "data:text/plain;charset=utf-8," + encodeURIComponent(resultStr);

    $("#download").attr("download", "webgltest.log");
    $("#download").attr("href", url);
    //alert(resultStr);

};

var styleSpec = ["line", "stick", "sphere", "cartoon"];


//Generic style render testcase

var testcase = function(styleType, profile) {
    
    var testName = styleType + " render";
    var timeName = styleType + " render time: ";
    var testMsg = styleType + " style set correctly";
    var script = "select all; spacefill off;";
    if (styleType === "line") {
        script += "wireframe;";
    }
    else if (styleType === "stick") {
        script += "wireframe 100;";
    }
    else if (styleType === "sphere") {
        script += "spacefill;";
    }
    else if (styleType === "cartoon")
        script += "cartoon;";
    
    test(testName, function() {
        console.group(testName);
        console.time(timeName);
        
        if (profile)
            console.profile();
        
        var msg = Jmol.scriptEcho(viewer, script);
        
        if (profile)
            console.profileEnd();
        console.log(msg);
        console.timeEnd(timeName);
        console.groupEnd();
        
        ok(true, testMsg); 
        
    });
};


var runtests = (function(profile) {
    for (var style in styleSpec)
        new testcase(styleSpec[style], profile);
});


QUnit.module( "C. Calicivirus Capsid, 12,362 atoms (3M8L)", {

    setupOnce: function() {
        
        console.group("Capsid (12,362 atoms)");
    },

    teardownOnce: function() {
        console.groupEnd();
            //glviewer.removeAllModels();
    }
});

runtests(profile);






