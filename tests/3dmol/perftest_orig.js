// Test rendering performance for different sized pdb's

// var glviewer = null;

$(document).ready(() => {

    glviewer = $3Dmol.createViewer("gldiv", {defaultcolors: $3Dmol.rasmolElementColors, callback(viewer) {
        viewer.setBackgroundColor(0xffffff);
    }});

    // starts QUnit tests
    start();
}); 

const profile = QUnit.urlParams.profilecheck;
let resultXML = null;
let resultStr = "";

// QUnit-reporter hook to output test results in XML format
QUnit.jUnitReport = function(data) {

	resultXML = $.parseXML(data.xml);
	
	// Wrap XML result in JQuery object; parse and setup output string
	$result = $(resultXML);
	
	resultStr += "$3Dmol Performance Tests\n";
	const runTime = $result.find("testsuites").attr("time");
	const runDate = $result.find("testsuites").attr("timestamp");
	resultStr += `Total Test Time: ${  runTime  } s\n`;
	resultStr += `Date: ${  runDate  }\n\n`;
	
	$result.find("testsuite").each(function(){
		const moduleName = $(this).attr("name");
		const moduleTime = $(this).attr("time");
		resultStr += `\n${  moduleName}`;
		// alert(moduleName);
		$(this).find("testcase").each(function() {
			const testName = $(this).attr("name");
			const testTime = $(this).attr("time");
			resultStr += `\n\t${  testName  }:   ${  testTime  } s`;
			// alert(testName);
		});
		resultStr += `\n\tTotal:         ${  moduleTime  } s\n`;
	});
	
	// Set up a link to download test results
	$("#qunit-testresult").append("<br><a id='download'>Download</a>");
	const url = `data:text/plain;charset=utf-8,${  encodeURIComponent(resultStr)}`;
	
	$("#download").attr("download", "webgltest.log");
	$("#download").attr("href", url);
	// alert(resultStr);

};

const styleSpec = {"stick":{stick:{}}, "line":{line:{}}, "cross":{cross:{}}, "sphere":{sphere:{}}, "cartoon":{cartoon:{color:"spectrum"}}};


// Generic style render testcase

const testcase = function(styleType, profile) {
    
    const testName = `${styleType  } render`;
    const timeName = `${styleType  } render time: `;
    const testMsg = `${styleType  } style set correctly`;
    const styleExpected = styleSpec[styleType];
    
    test(testName, () => {
        const m = glviewer.getModel(0);
        console.group(testName);
        console.time(timeName);
        
        if (profile)
            console.profile();
            
        glviewer.setStyle({}, styleExpected);
        glviewer.render();
        
        if (profile)
            console.profileEnd();
        
        console.timeEnd(timeName);
        console.groupEnd();
        
        const styleActual = m.selectedAtoms()[0].style;
        equal(JSON.stringify(styleActual), JSON.stringify(styleExpected), testMsg); 
        
    });
};


const runtests = (function(profile) {
	for (const style in styleSpec){
		const tc = new testcase(style, profile);
	}
});
// test cases

// TESTS

// moldata 1

QUnit.module( "A. Bovine Calbindin, 637 atoms (1YCR)", {
	
	setupOnce() {
		
		glviewer.removeAllModels();
		stop();
   		$.get("test_structs/1YCR.pdb", (data) => {				
	      		glviewer.addModel(data, "pdb");
	      		glviewer.zoomTo();
	      		glviewer.render();
	      		start();
   		}, "text");
		
		console.group("Calbindin (637 atoms)");
	},
		
	teardownOnce() {		
		console.groupEnd();
	}
	
});

// new tester("sphere", profile);
runtests(profile);

// moldata 2

QUnit.module( "B. Cathodic Hemoglobin, 5,085 atoms (2AA1)", {
	
	setupOnce() {
		glviewer.removeAllModels();
		stop();
   		$.get("test_structs/2AA1.pdb", (data) => {				
	      		glviewer.addModel(data, "pdb");
	      		glviewer.zoomTo();
	      		glviewer.render();
	      		start();
   		}, "text");
   		console.groupEnd();
   		console.group("Hemoglobin (5,085 atoms)");
	},
	
	teardownOnce() {
		console.groupEnd();
	}
});

runtests(profile);

// moldata 3

QUnit.module( "C. Calicivirus Capsid, 12,362 atoms (3M8L)", {
	
	setupOnce() {
		console.log("Testing third molecule");
		glviewer.removeAllModels();
		stop();
   		$.get("test_structs/3M8L.pdb", (data) => {
	      		glviewer.addModel(data, "pdb");
	      		glviewer.zoomTo();
	      		glviewer.render();
	      		start();
   		}, "text");
   		console.groupEnd();
   		console.group("Capsid (12,362 atoms)");
	},
		
	teardownOnce() {
		console.groupEnd();
		// glviewer.removeAllModels();
	}
});

runtests(profile);


QUnit.module( "D. Dihydrolipoyl Transacetylase Biological Assembly,  71,820 atoms (1B5S)", {
    
    setupOnce() {
        
        glviewer.removeAllModels();
        stop();
        $.get("../../test_structs/1B5S.pdb", (data) => {             
                glviewer.addModel(data, "pdb");
                glviewer.zoomTo();
                glviewer.render();
                start();
        }, "text");
        
        console.group("Calbindin (71,820 atoms)");
    },
        
    teardownOnce() {      
        console.groupEnd();
    }
    
});

// runtests(profile);





