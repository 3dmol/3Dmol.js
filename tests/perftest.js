//Test rendering performance for different sized pdb's

QUnit.config.autostart = false;

var styleSpec = {"stick":{stick:{}}, "line":{line:{}}, "cross":{cross:{}}, "sphere":{sphere:{}}, "cartoon":{cartoon:{color:0x0000ff}}};

//test cases
var runtests = function() {

//setup new model

	

//Render stick - use some console logging
test("Stick Render", function() {
	var m = glviewer.getModel(0);
	var styleExpected = styleSpec["stick"];
	
	console.group("Stick Render");
	console.time("Stick render time: ");
	
	glviewer.setStyle({},styleExpected);
	glviewer.render();
	
	console.timeEnd("Stick render time: ");
	console.groupEnd();
	
	//Should return first atom's style from our model
	var styleActual = m.selectedAtoms()[0].style;
	
	equal(JSON.stringify(styleActual), JSON.stringify(styleExpected), "Stick style set correctly");	
});

//Render cross - use some console logging
test("Cross Render", function() {
	var m = glviewer.getModel(0);
	var styleExpected = styleSpec["cross"];
	
	console.group("Cross Render");
	console.time("Cross render time: ");
	
	glviewer.setStyle({},styleExpected);
	glviewer.render();
	
	console.timeEnd("Cross render time: ");
	console.groupEnd();
	
	//Should return first atom's style from our model
	var styleActual = m.selectedAtoms()[0].style;
	
	equal(JSON.stringify(styleActual), JSON.stringify(styleExpected), "Cross style set correctly");	
});

//Render sphere - use some console logging
test("Sphere Render", function() {
	var m = glviewer.getModel(0);
	var styleExpected = styleSpec["sphere"];
	
	console.group("Sphere Render");
	console.time("Sphere render time: ");
	
	glviewer.setStyle({},styleExpected);
	glviewer.render();
	
	console.timeEnd("Sphere render time: ");
	console.groupEnd();
	
	//Should return first atom's style from our model
	var styleActual = m.selectedAtoms()[0].style;
	
	equal(JSON.stringify(styleActual), JSON.stringify(styleExpected), "Sphere style set correctly");	
});

//Render cartoon - use some console logging
test("Cartoon Render", function() {
	var m = glviewer.getModel(0);
	var styleExpected = styleSpec["cartoon"];
	
	console.group("Cartoon Render");
	console.time("Cartoon render time: ");
	
	glviewer.setStyle({},styleExpected);
	glviewer.render();
	
	console.timeEnd("Cartoon render time: ");
	console.groupEnd();
	
	//Should return first atom's style from our model
	var styleActual = m.selectedAtoms()[0].style;
	
	equal(JSON.stringify(styleActual), JSON.stringify(styleExpected), "Cartoon style set correctly");	
});

};
//TESTS



//moldata 1

QUnit.module( "Bovine Calbindin, 76 res (1YCR)", {
	
	setupOnce: function() {
		console.group("MOL 1");
		console.log("Testing first molecule");
		glviewer.removeAllModels();
		var moldata = $("#moldata_1").val();
		console.log("moldata length: " + moldata.length);
		var m = glviewer.addModel(moldata, "pdb");
		glviewer.mapAtomProperties(WebMol.partialCharges);
		glviewer.zoomTo();
		glviewer.render();
	},
	
	teardownOnce: function() {
		console.groupEnd();
		glviewer.removeAllModels();
	}
	
});

runtests();

//moldata 2

QUnit.module( "Cathodic Hemoglobin, 143 res (2AA1)", {
	
	setupOnce: function() {
		console.group("MOL 2");
		console.log("Testing second molecule");
		glviewer.removeAllModels();
	},
	
	teardownOnce: function() {
		console.groupEnd();
		glviewer.removeAllModels();
	}
});

//I'm loading in a molecule from PDB because it's too large to add to html file directly
// Have to wait to make sure it's ajax request is finished before resuming tests
//TODO: See if there's a way to include these files and not have to load them
asyncTest("Load molecule", function() {
	WebMol.download("pdb:2AA1", glviewer);
	setTimeout(function() {
		ok(true, "resuming");
		start();
	}, 1000);
		
});

runtests();

//moldata 3

QUnit.module( "Calicivirus Capsid, 534 res (3M8L)", {
	
	setupOnce: function() {
		console.group("MOL 3");
		console.log("Testing second molecule");
		glviewer.removeAllModels();
	},
	
	teardownOnce: function() {
		console.groupEnd();
		glviewer.removeAllModels();
	}
});

asyncTest("Load molecule", function() {
	WebMol.download("pdb:3M8L", glviewer);
	setTimeout(function() {
		ok(true, "resuming");
		start();
	}, 1000);
		
});

runtests();





