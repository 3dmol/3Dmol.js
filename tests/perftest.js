//Test rendering performance for different sized pdb's

QUnit.config.autostart = false;

var styleSpec = {"stick":{stick:{}}, "line":{line:{}}, "cross":{cross:{}}, "sphere":{sphere:{}}, "cartoon":{cartoon:{color:0x0000ff}}};
var glviewer = null;


var runtests = function() {
	
test("Initial render", function(){ 
 	ok(glviewer);
});	



test("Model set correctly", function() {
	var m = glviewer.getModel(0);
	ok(m);
});


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
module( "First", {
	setup: function() {
		glviewer.removeAllModels();
		var moldata = $("#moldata_1").val();
        var m = glviewer.addModel(moldata,"pdb");
        glviewer.mapAtomProperties(WebMol.partialCharges);
        glviewer.zoomTo();
        glviewer.render();
	}
});

runtests();

module( "First", {
	setup: function() {
		glviewer.removeAllModels();
		var moldata = $("#moldata_12").val();
        var m = glviewer.addModel(moldata,"pdb");
        glviewer.mapAtomProperties(WebMol.partialCharges);
        glviewer.zoomTo();
        glviewer.render();
	}
});

runtests();


//Look over getting the closure variables set correctly...
/*
for (var style in styleSpec) (function(style)
{
	
	var testName = style + " render";
	var consoleMsg = style + " render time: ";
	test(testName, function(){
		var styleExpected = styleSpec[style];
		
		var m = glviewer.getModel(0);
		//var consoleMsg = style + " render time: ";
		
		console.group(testName);
		console.time(consoleMsg);
		
		glviewer.setStyle(styleExpected);
		glviewer.render();
		
		console.timeEnd(consoleMsg);
		console.groupEnd();

		//Should return first atom's style from our model
		var styleActual = m.selectedAtoms()[0].style;	
		//alert(JSON.stringify(styleActual));
		equal(JSON.stringify(styleActual), JSON.stringify(styleExpected));			
	});
});
*/
