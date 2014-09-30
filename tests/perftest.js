/* 
 * QUnit benchmark tests for WebMol, GLmol, and JSmol
 */
//Test rendering performance for different sized pdb's 

var profile = QUnit.urlParams.profilecheck;

var testSuite = testSuite || "webmol";

// Style types to test
var styleSpec = ["line", "stick", "sphere", "cartoon"];

// WebMol testcase generator

var genWebMolTestCase = function(styleType) {
    
    var testName = styleType + " render";
    var timeMsg = styleType + " render time: ";
    var testMsg = styleType + " style set correctly";
    var style = {}; style[styleType] = {};
    
    QUnit.test(testName, function() {
        
        viewer.setStyle({}, {cross:{}});
        viewer.render();
        console.group(testName);          
        
        var start = new Date();
        
        viewer.setStyle({}, style);           
        viewer.render();
        
        var end = new Date();
        var testTime = end - start;
        
        resultTimes[testName] = testTime;
        console.log(timeMsg + (testTime) + "ms");
        console.groupEnd();
        
        QUnit.ok(true, testMsg);
        
    });
};

//GLmol test generator

var genGLmolTestCase = function(styleType, profile) {
    
    var testName = styleType + " render";
    var timeMsg = styleType + " render time: ";
    var testMsg = styleType + " style set correctly";
    
    var defineRep = function(style){
        var all = viewer.getAllAtoms();
        
        return function() {
            //var all = this.getAllAtoms();
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
        
        var view = viewer.getView();
        viewer.initializeScene();
        var all = viewer.getAllAtoms();
        viewer.zoomInto(all);
        viewer.defineRepresentation = defineRep(styleType);
        
        var start = new Date();
        //Draw appropriate style        
        viewer.defineRepresentation();
        viewer.setView(view);        
        viewer.show();     
        
        var end = new Date();
        var testTime = end - start;
        
        resultTimes[testName] = testTime;
        console.log(timeMsg + (testTime) + "ms");        
        console.groupEnd();
             
        ok(true, testMsg); 
        
    });
};

//JSmol testcase generator

var genJSmolTestCase = function(styleType, profile) {
    
    var testName = styleType + " render";
    var timeMsg = styleType + " render time: ";
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
        script += "set cartoonFancy true; cartoon only;";
      
    //Create test case
    QUnit.test(testName, function() {
        
        Jmol.scriptWait(viewer, "select *; wireframe off; cartoon off; spacefill off;");
        
        console.group(testName);
        
        var start = new Date();
        
        Jmol.scriptWait(viewer, script);
        
        var end = new Date();      
        
        var testTime = end - start;
        
        resultTimes[testName] = testTime;
        console.log(timeMsg + (testTime) + "ms");          
        console.groupEnd();
        
        QUnit.ok(true, testMsg); 
        
    });
};


// Setup test modules


//WebMol tests
if (testSuite === 'webmol') {
    QUnit.module( "WebMol Tests", {

        setupOnce: function() {
            viewer.zoomTo();
            console.group("WebMol");
        },

        teardownOnce: function() {
            console.groupEnd();
        }

    });

    // WebMol test cases
    for (var style in styleSpec)
        genWebMolTestCase(styleSpec[style], profile);    

    //surface test

    var testName = "SURF render";
    var timeMsg = "surface render time: ";
    var testMsg = "surface style set correctly";
    
    QUnit.test(testName, function() {
        
        viewer.setStyle({}, {});
        viewer.render();
        console.group(testName);          
        
        var start = new Date();
        
        viewer.addSurface(WebMol.SurfaceType.VDW, {}, {}, {}, {});           
        viewer.render();
        
        var end = new Date();
        var testTime = end - start;
        
        resultTimes[testName] = testTime;
        console.log(timeMsg + (testTime) + "ms");
        console.groupEnd();
        
        QUnit.ok(true, testMsg);
        
    });

}


else if (testSuite === "glmol") {
//GLmol testing module

    QUnit.module( "GLmol Tests", {

        setupOnce: function() {
            var all = viewer.getAllAtoms();
            viewer.zoomInto(all);
            console.group("GLmol");
        },

        teardownOnce: function() {
            console.groupEnd();
            //glviewer.removeAllModels();
        }
    });

    for (var style in styleSpec)
        genGLmolTestCase(styleSpec[style], profile);  

    
}

else if (testSuite === "glmol_surf") {

    QUnit.module( "GLmol Surface Test", {

        setupOnce: function() {
            var all = viewer.getAllAtoms();
            viewer.zoomInto(all);
            console.group("GLmol");
        },

        teardownOnce: function() {
            console.groupEnd();
        }
    });


    var testName = "surf render";
    var timeMsg = "surf render time: ";
    var testMsg = "surf style set correctly";
    
    var defineRep = function(){
        var all = viewer.getAllAtoms();
        
        return function() {
            //var all = this.getAllAtoms();
            var target = this.modelGroup;          
            this.drawAsCross(target, all, 0.3, true);

            
        }; 
        
    };
    
    QUnit.test(testName, function() {
        
        console.group(testName);
        
        var view = viewer.getView();
        viewer.initializeScene();
        var all = viewer.getAllAtoms();
        viewer.zoomInto(all);
        viewer.defineRepresentation = defineRep();
        
        
        //Draw appropriate style        
        //viewer.defineRepresentation();
        //viewer.drawAsCross(viewer.modelGroup, all, 0.3, true);
        var start = new Date();
        viewer.generateMesh(viewer.modelGroup, all, 1, false);
        
        viewer.setView(view);        
        viewer.show();     
        
        var end = new Date();
        var testTime = end - start;
        
        resultTimes[testName] = testTime;
        console.log(timeMsg + (testTime) + "ms");        
        console.groupEnd();
             
        ok(true, testMsg); 
        
    });

}

else if (testSuite === "jmol") {
    
    //JSmol Tests
    QUnit.module( "Jmol Tests", {

        setupOnce: function() {   
            Jmol.scriptWait(viewer, "load test_structs/3M8L.pdb");
            var end = new Date();
            
            resultTimes['initialization'] = end - start;
            console.group("JSmol");  
            console.log("starting JSmol tests");
        },

        teardownOnce: function() {
            console.groupEnd();
            //$("#viewerdiv").empty();
        }

    });

    for (var style in styleSpec)
        genJSmolTestCase(styleSpec[style], profile);


    // Surf test
    var testName = "surface render";
    var timeMsg = "surface render time: ";
    var testMsg = "surface style set correctly";
      
    //Create test case
    QUnit.test(testName, function() {
        
        Jmol.scriptWait(viewer, "select *; wireframe off; cartoon off; spacefill off;");
        
        console.group(testName);
        
        var start = new Date();
        
        Jmol.scriptWait(viewer, "select *; isosurface vdw;");
        
        var end = new Date();      
        
        var testTime = end - start;
        
        resultTimes[testName] = testTime;
        console.log(timeMsg + (testTime) + "ms");          
        console.groupEnd();
        
        QUnit.ok(true, testMsg); 
        
    });
    
}







