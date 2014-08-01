//Test rendering performance for different sized pdb's 

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
    var script = "select *;";
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
        script += "cartoon only;";
    
    window["endTime"] = function(a, b, c, d) {
        console.log(a);
        console.log(b);
    };    
    
    QUnit.test(testName, function() {

        console.group(testName);
        console.time(timeName);
        
        if (profile)
            console.profile();
        
        var arr = Jmol.scriptWait(viewer, script);
        
        while ((arr.length < 3))
            continue;
        
        if (profile)
            console.profileEnd();
        
        console.timeEnd(timeName);
        console.groupEnd();
        
        ok(true, testMsg); 
        
    });
};

QUnit.module( "C. Calicivirus Capsid, 12,362 atoms (3M8L)", {

    setupOnce: function() {
        
        console.group("Capsid (12,362 atoms)");
        
    },

    teardownOnce: function() {
        console.groupEnd();
            //glviewer.removeAllModels();
    },
    
    setup: function() {
        
        Jmol.scriptWait(viewer, "wireframe -0.1; spacefill off; cartoon off; set cartoonFancy true;");
        
    }
});

for (var style in styleSpec)
    testcase(styleSpec[style], profile);






