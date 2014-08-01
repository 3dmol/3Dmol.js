/* 
 * QUnit benchmark tests for WebMol, GLmol, and JSmol
 */
//Test rendering performance for different sized pdb's 

var profile = QUnit.urlParams.profilecheck;

//QUnit-reporter hook to output test results in XML format
QUnit.jUnitReport = function(data) {
    
    var resultXML = $.parseXML(data.xml);
    //Wrap XML result in JQuery object; parse and setup output string
    var result = $(resultXML);

    var resultStr = "Viewer Performance Tests: Calicivirus Capsid, 12,362 atoms (3M8L)\n";
    var runTime = result.find("testsuites").attr("time");
    var runDate = result.find("testsuites").attr("timestamp");
    resultStr += "Total Test Time: " + runTime + " s\n";
    resultStr += "Date: " + runDate + "\n\n";

    result.find("testsuite").each(function(){
        var moduleName = $(this).attr("name");
        var moduleTime = $(this).attr("time");
        resultStr += "\n" + moduleName;
        //alert(moduleName);
        $(this).find("testcase").each(function() {
            var testName = $(this).attr("name");
            var testTime = $(this).attr("time");
            resultStr += "\n\t" + testName + ":   " + testTime + " s";
            //alert(testName);
        });      
        
        resultStr += "\n\tTotal:         " + moduleTime + " s\n";
    });

    //Set up a link to download test results
    $("#qunit-testresult").append("<br><a id='download'>Download</a>");
    var url = "data:text/plain;charset=utf-8," + encodeURIComponent(resultStr);

    $("#download").attr("download", "webgltest.log");
    $("#download").attr("href", url);
    //alert(resultStr);

};

// Style types to test
var styleSpec = ["line", "stick", "sphere"];

// WebMol testcase generator

var genWebMolTestCase = function(styleType, profile) {
    
    var testName = styleType + " render";
    var timeMsg = styleType + " render time: ";
    var testMsg = styleType + " style set correctly";
    var style = {}; style[styleType] = {};
    
    test(testName, function() {
        
        viewer.setStyle({}, {cross:{}});
        viewer.render();
        console.group(testName);   
        
        console.time(timeMsg);
        
        if (profile)
            console.profile();
        
        var start = new Date();
        viewer.setStyle({}, style);           
        viewer.render();
        var end = new Date();
        var testTime = end - start;
        console.timeEnd(timeMsg);
        
        if (profile)
            console.profileEnd();
        
        console.log(timeMsg + (testTime) + "ms");
        console.groupEnd();
        
        QUnit.ok(true, testMsg);
        
    });
};

//GLmol test generator

var genGLmolTestCase = function(styleType, profile) {
    
    var testName = styleType + " render";
    var timeName = styleType + " render time: ";
    var testMsg = styleType + " style set correctly";
    
    var defineRep = function(style){
    
        return function() {
            var all = this.getAllAtoms();
            var target = this.modelGroup;

            if (style === 'stick') 
                this.drawBondsAsStick(target, all, this.cylinderRadius, this.cylinderRadius, true);
            else if (style === 'line')
                this.drawBondsAsLine(target, all, this.lineWidth);
            else if (style === 'cross') 
                this.drawAsCross(target, all, 0.3, true);
            else if (style === 'sphere') 
                this.drawAtomsAsSphere(target, all, this.sphereRadius);
            else if (style === 'cartoon') {
                //this.colorChainbow(all);
                this.drawCartoon(target, all, false, this.thickness);
            }
        }; 
        
    };
    
    QUnit.test(testName, function() {
        
        console.group(testName);
        console.time(timeName);
               
        if (profile)
            console.profile();
            
        viewer.defineRepresentation = defineRep(styleType);
        viewer.rebuildScene();
        viewer.show();       
        
        if (profile)
            console.profileEnd();
        
        console.timeEnd(timeName);
        console.groupEnd();
             
        ok(true, testMsg); 
        
    });
};

//JSmol testcase generator

var genJSmolTestCase = function(styleType, profile) {
    
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
      
    //Create test case
    QUnit.test(testName, function() {

        console.group(testName);
        console.time(timeName);
        
        if (profile)
            console.profile();
        
        Jmol.scriptWait(viewer, script);
        
        if (profile)
            console.profileEnd();
        
        console.timeEnd(timeName);
        console.groupEnd();
        
        QUnit.ok(true, testMsg); 
        
    });
};


// Setup test modules


//WebMol tests

QUnit.module( "WebMol Tests", {
	
    setupOnce: function() {
        QUnit.stop();
        $("#viewerdiv").empty();
        viewer = WebMol.createViewer("viewerdiv");
        viewer.setBackgroundColor(0xffffff);
        $.get("test_structs/3M8L.pdb", function(data) {
                viewer.addModel(data, "pdb");
                viewer.zoomTo();
                QUnit.start();
        }, "text");

        console.group("WebMol");
    },
		
    teardownOnce: function() {
        console.groupEnd();
    }
    
});

// WebMol test cases
for (var style in styleSpec)
    genWebMolTestCase(styleSpec[style], profile);


//GLmol testing module

QUnit.module( "GLmol Tests", {
	
    setupOnce: function() {
        QUnit.stop();
        $("#viewerdiv").empty();
        
        viewer = new GLmol("viewerdiv", true);
        viewer.initializeScene();
        viewer.setBackground(0xffffff);
        
        $.get("test_structs/3M8L.pdb", function(data) {
            
            viewer.protein = {sheet: [], helix: [], biomtChains: '', biomtMatrices: [], symMat: [], pdbID: '', title: ''};
            viewer.atoms = [];
            viewer.parsePDB2(data);
            var all = viewer.getAllAtoms();
            viewer.colorByAtom(all, {});
            viewer.zoomInto(all);
            viewer.drawAsCross(viewer.modelGroup, all, 0.3, true);
            QUnit.start();
            
        }, "text");
        
        
        console.group("GLmol");
    },

    teardownOnce: function() {
        console.groupEnd();
        //glviewer.removeAllModels();
    }
});

for (var style in styleSpec)
    genGLmolTestCase(styleSpec[style], profile);

//JSmol Tests
QUnit.module( "JSmol Tests", {

    setupOnce: function() {    
        QUnit.stop();
        if (viewer !== undefined && viewer instanceof WebMol.GLViewer)
            viewer.removeAllModels();
        $("#viewerdiv").empty();
        console.group("JSmol");  
        console.log("starting JSmol tests");
        $("#viewerdiv").html(Jmol.getAppletHtml("viewer", Info));
        console.log("Filled html");
    },

    teardownOnce: function() {
        console.groupEnd();
        //$("#viewerdiv").empty();
    },
    
    setup: function() {
        console.log("setting up test...");
        //Jmol.scriptWait(viewer, "wireframe -0.1; spacefill off; cartoon off; set cartoonFancy true;");       
    }
    
});

for (var style in styleSpec)
    genJSmolTestCase(styleSpec[style], profile);






