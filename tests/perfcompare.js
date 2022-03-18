/* 
 * QUnit benchmark tests for $3Dmol, GLmol, and JSmol
 */
// Test rendering performance for different sized pdb's 

const profile = QUnit.urlParams.profilecheck;

// QUnit-reporter hook to output test results in XML format
QUnit.jUnitReport = function(data) {
    
    const resultXML = $.parseXML(data.xml);
    // Wrap XML result in JQuery object; parse and setup output string
    const result = $(resultXML);

    let resultStr = "Viewer Performance Tests: Calicivirus Capsid, 12,362 atoms (3M8L)\n";
    const runTime = result.find("testsuites").attr("time");
    const runDate = result.find("testsuites").attr("timestamp");
    resultStr += `Total Test Time: ${  runTime  } s\n`;
    resultStr += `Date: ${  runDate  }\n\n`;

    result.find("testsuite").each(function(){
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

// Style types to test
const styleSpec = ["line", "stick", "sphere"];

// $3Dmol testcase generator

const gen$3DmolTestCase = function(styleType, profile) {
    
    const testName = `${styleType  } render`;
    const timeMsg = `${styleType  } render time: `;
    const testMsg = `${styleType  } style set correctly`;
    const style = {}; style[styleType] = {};
    
    test(testName, () => {
        
        viewer.setStyle({}, {cross:{}});
        viewer.render();
        console.group(testName);   
        
        console.time(timeMsg);
        
        if (profile)
            console.profile();
        
        const start = new Date();
        viewer.setStyle({}, style);           
        viewer.render();
        const end = new Date();
        const testTime = end - start;
        console.timeEnd(timeMsg);
        
        if (profile)
            console.profileEnd();
        
        console.log(`${timeMsg + (testTime)  }ms`);
        console.groupEnd();
        
        QUnit.ok(true, testMsg);
        
    });
};

// GLmol test generator

const genGLmolTestCase = function(styleType, profile) {
    
    const testName = `${styleType  } render`;
    const timeName = `${styleType  } render time: `;
    const testMsg = `${styleType  } style set correctly`;
    
    const defineRep = function(style){
    
        return function() {
            const all = this.getAllAtoms();
            const target = this.modelGroup;

            if (style === 'stick') 
                this.drawBondsAsStick(target, all, this.cylinderRadius, this.cylinderRadius, true);
            else if (style === 'line')
                this.drawBondsAsLine(target, all, this.lineWidth);
            else if (style === 'cross') 
                this.drawAsCross(target, all, 0.3, true);
            else if (style === 'sphere') 
                this.drawAtomsAsSphere(target, all, this.sphereRadius);
            else if (style === 'cartoon') {
                // this.colorChainbow(all);
                this.drawCartoon(target, all, false, this.thickness);
            }
        }; 
        
    };
    
    QUnit.test(testName, () => {
        
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

// JSmol testcase generator

const genJSmolTestCase = function(styleType, profile) {
    
    const testName = `${styleType  } render`;
    const timeName = `${styleType  } render time: `;
    const testMsg = `${styleType  } style set correctly`;
    let script = "select *;";
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
      
    // Create test case
    QUnit.test(testName, () => {

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


// $3Dmol tests

QUnit.module( "$3Dmol Tests", {
	
    setupOnce() {
        QUnit.stop();
        $("#viewerdiv").empty();
        viewer = $3Dmol.createViewer("viewerdiv");
        viewer.setBackgroundColor(0xffffff);
        $.get("test_structs/3M8L.pdb", (data) => {
                viewer.addModel(data, "pdb");
                viewer.zoomTo();
                QUnit.start();
        }, "text");

        console.group("$3Dmol");
    },
		
    teardownOnce() {
        console.groupEnd();
    }
    
});

// $3Dmol test cases
for (const style in styleSpec)
    gen$3DmolTestCase(styleSpec[style], profile);


// GLmol testing module

QUnit.module( "GLmol Tests", {
	
    setupOnce() {
        QUnit.stop();
        $("#viewerdiv").empty();
        
        viewer = new GLmol("viewerdiv", true);
        viewer.initializeScene();
        viewer.setBackground(0xffffff);
        
        $.get("test_structs/3M8L.pdb", (data) => {
            
            viewer.protein = {sheet: [], helix: [], biomtChains: '', biomtMatrices: [], symMat: [], pdbID: '', title: ''};
            viewer.atoms = [];
            viewer.parsePDB2(data);
            const all = viewer.getAllAtoms();
            viewer.colorByAtom(all, {});
            viewer.zoomInto(all);
            viewer.drawAsCross(viewer.modelGroup, all, 0.3, true);
            QUnit.start();
            
        }, "text");
        
        
        console.group("GLmol");
    },

    teardownOnce() {
        console.groupEnd();
        // glviewer.removeAllModels();
    }
});

for (const style in styleSpec)
    genGLmolTestCase(styleSpec[style], profile);

// JSmol Tests
QUnit.module( "JSmol Tests", {

    setupOnce() {    
        QUnit.stop();
        if (viewer !== undefined && viewer instanceof $3Dmol.GLViewer)
            viewer.removeAllModels();
        $("#viewerdiv").empty();
        console.group("JSmol");  
        console.log("starting JSmol tests");
        $("#viewerdiv").html(Jmol.getAppletHtml("viewer", Info));
        console.log("Filled html");
    },

    teardownOnce() {
        console.groupEnd();
        // $("#viewerdiv").empty();
    },
    
    setup() {
        console.log("setting up test...");
        // Jmol.scriptWait(viewer, "wireframe -0.1; spacefill off; cartoon off; set cartoonFancy true;");       
    }
    
});

for (const style in styleSpec)
    genJSmolTestCase(styleSpec[style], profile);






