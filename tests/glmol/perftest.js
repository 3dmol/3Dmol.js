//GLMol performance tests

$(document).ready(function() {

    glviewer = new GLmol("gldiv", true);
    //starts QUnit tests
    start();
}); 

var resultXML = null;
var resultStr = "";
var profile = QUnit.urlParams.profilecheck;

//QUnit-reporter hook to output test results in XML format
QUnit.jUnitReport = function(data) {

	resultXML = $.parseXML(data.xml);
	
	//Wrap XML result in JQuery object; parse and setup output string
	$result = $(resultXML);
	
	resultStr += "GLMol Performance Tests\n";
	var runTime = $result.find("testsuites").attr("time");
	var runDate = $result.find("testsuites").attr("timestamp");
	resultStr += "Total Test Time: " + runTime + " s\n";
	resultStr += "Date: " + runDate + "\n\n";
	
	$result.find("testsuite").each(function(){
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

var styleSpec = ["stick", "line", "cross", "sphere", "cartoon"];

var defineRep = function(style) {
    
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
            this.colorChainbow(all);
            this.drawCartoon(target, all, false, this.thickness);
        }
    };
    
};

//Generic style render testcase

var testcase = function(styleType, profile) {
    
    var testName = styleType + " render";
    var timeName = styleType + " render time: ";
    var testMsg = styleType + " style set correctly";
    
    
    test(testName, function() {
        
        console.group(testName);
        console.time(timeName);
        
        
        if (profile)
            console.profile();
            
        glviewer.defineRepresentation = defineRep(styleType);
        glviewer.rebuildScene();
        glviewer.show();       
        
        if (profile)
            console.profileEnd();
        
        console.timeEnd(timeName);
        console.groupEnd();
        
        
        ok(true, testMsg); 
        
    });
};


var runtests = (function(profile) {
	for (var style in styleSpec)
		new testcase(styleSpec[style], profile);
});
//test cases

//TESTS


//moldata 3

QUnit.module( "C. Calicivirus Capsid, 12,362 atoms (3M8L)", {
	
    setupOnce: function() {

        stop();
        $.get("test_structs/3M8L.pdb", function(data) {
                //glviewer.loadMoleculeStr(false, data);
            glviewer.protein = {sheet: [], helix: [], biomtChains: '', biomtMatrices: [], symMat: [], pdbID: '', title: ''};
            glviewer.atoms = [];
            glviewer.parsePDB2(data);
            var all = glviewer.getAllAtoms();
            glviewer.colorByAtom(all, {});
            glviewer.initializeScene();
            glviewer.setBackground(0xffffff);
            glviewer.zoomInto(glviewer.getAllAtoms());

                start();
        }, "text");
        console.groupEnd();
        console.group("Capsid (12,362 atoms)");
    },

    teardownOnce: function() {
        console.groupEnd();
        //glviewer.removeAllModels();
    }
});

runtests(profile);






