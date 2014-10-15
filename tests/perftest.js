/* 
 * QUnit benchmark tests for $3Dmol, GLmol, and JSmol
 */
//Test rendering performance for different sized pdb's 

var profile = QUnit.urlParams.profilecheck;

var testSuite = testSuite || "3Dmol";

// Style types to test
var styleSpec = ["line", "stick", "sphere", "cartoon"];

// $3Dmol testcase generator

var gen$3DmolTestCase = function(styleType) {
    
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
        
        viewer.rotate(10);
        var time2 = new Date() - end;
        
        resultTimes[testName] = testTime;
        console.log(timeMsg + (testTime) + "ms");
        console.log(timeMsg + (testTime) + "ms; rotate "+time2+"ms");          
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
    var script = "";
    if (styleType === "line") {
        script += "wireframe only;";
    }
    else if (styleType === "stick") {
        script += "wireframe 100 only;";
    }
    else if (styleType === "sphere") {
        script += "spacefill only;";
    }
    else if (styleType === "cartoon")
        script += "set cartoonFancy true; cartoon only;";
      
    //Create test case
    QUnit.test(testName, function() {
        
        console.group(testName);
        var start = new Date();
        Jmol.scriptWait(viewer, script);
        
        Jmol.scriptWaitOutput(viewer, "refresh;");
        
        var end = new Date();      
        
        var testTime = end - start;
        
        Jmol.scriptWaitOutput(viewer,"rotate 10; refresh;");
        
        var time2 = new Date() - end;
        resultTimes[testName] = testTime;
        console.log(timeMsg + (testTime) + "ms; rotate "+time2+"ms");          
        console.groupEnd();
        
        QUnit.ok(true, testMsg); 
        
    });
};


// Setup test modules


//$3Dmol tests
if (testSuite === '3Dmol') {
    QUnit.module( "$3Dmol Tests", {

        setupOnce: function() {
            viewer.zoomTo();
            console.group("$3Dmol");
        },

        teardownOnce: function() {
            console.groupEnd();
        }

    });

    // $3Dmol test cases 
    for (var style in styleSpec)
        gen$3DmolTestCase(styleSpec[style], profile);    
    //surface test
/*
    var testName = "SURF render";
    var timeMsg = "surface render time: ";
    var testMsg = "surface style set correctly";
    
    QUnit.test(testName, function() {
        
        viewer.setStyle({}, {});
        viewer.render();
        console.group(testName);          
        
        var start = new Date();
        
        viewer.addSurface($3Dmol.SurfaceType.VDW, {}, {}, {}, {});           
        viewer.render();
        
        var end = new Date();
        var testTime = end - start;
        
        resultTimes[testName] = testTime;
        console.log(timeMsg + (testTime) + "ms");
        console.groupEnd();
        
        QUnit.ok(true, testMsg);
        
    });
 */   
        
    //Combo testcase
	var comboName = "combined";
    QUnit.test(comboName, function() {
        
        viewer.setStyle({}, {});
        viewer.removeAllSurfaces();
        viewer.render();        
        console.group(comboName);
        var start = new Date();

        viewer.setStyle({chain: 'A'}, {sphere: {} });
        viewer.setStyle({chain: 'B'}, {cartoon: {color: 'spectrum'} });
        viewer.setStyle({chain: 'C'}, {stick: {} });
        viewer.render();
        //Jmol.scriptWait(viewer, "select chain=A; spacefill; select chain=B; cartoon; color group; wireframe; select chain=C; wireframe 80;");
        var end = new Date();      
        
        var testTime = end - start;
        
        resultTimes[comboName] = testTime;
        console.log(timeMsg + (testTime) + "ms");          
        console.groupEnd();
        
        QUnit.ok(true, testMsg); 
    });
    
    var rotateName = "rotate";
    QUnit.test(rotateName, function() {
        

        console.group(rotateName);
        var start = new Date();
        
        viewer.rotate(10);

        var end = new Date();      
        
        var testTime = end - start;
        
        resultTimes[rotateName] = testTime;
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
    var combinedName = "combined";
    QUnit.test(combinedName, function() {
        
        console.group(combinedName);
        
        var view = viewer.getView();
        viewer.initializeScene();
        var all = viewer.getAllAtoms();
		var chainA = [];
		var chainB = [];
		var chainC = [];
		
		for(var i = 0; i < all.length; i++) {
			if(viewer.atoms[all[i]].chain == 'A') chainA.push(all[i]);
			else if(viewer.atoms[all[i]].chain == 'B') chainB.push(all[i]);
			else if(viewer.atoms[all[i]].chain == 'C') chainC.push(all[i]);
		}
        
        var start = new Date();
        //Draw appropriate style        
          //      viewer.setStyle({chain: 'A'}, {sphere: {} });
        //viewer.setStyle({chain: 'B'}, {cartoon: {color: 'spectrum'} });
        viewer.drawAtomsAsSphere(viewer.modelGroup, chainA, viewer.sphereRadius);
        viewer.colorChainbow(chainB);
        viewer.drawCartoon(viewer.modelGroup, chainB, false, viewer.thickness);
        viewer.drawBondsAsStick(viewer.modelGroup, chainC, viewer.cylinderRadius, viewer.cylinderRadius, true);
        
        viewer.zoomInto(all);
        viewer.show();     
        
        var end = new Date();
        var testTime = end - start;
        
        resultTimes[combinedName] = testTime;
        console.log(timeMsg + (testTime) + "ms");        
        console.groupEnd();
             
        QUnit.ok(true, testMsg); 
        
    });
    
        var rotateNAme = "rotate";
    QUnit.test(rotateName, function() {
        
        console.group(rotateName);
        
        var view = viewer.getView();
        var dx = .1;
      var dy = 0;
      var r = Math.sqrt(dx * dx + dy * dy);

           var rs = Math.sin(r * Math.PI) / r;


        var start = new Date();
        //Draw appropriate style        
          //      viewer.setStyle({chain: 'A'}, {sphere: {} });
        //viewer.setStyle({chain: 'B'}, {cartoon: {color: 'spectrum'} });
                 viewer.dq.x = Math.cos(r * Math.PI); 
         viewer.dq.y = 0;
         viewer.dq.z =  rs * dx; 
         viewer.dq.w =  rs * dy;
         viewer.rotationGroup.quaternion = new THREE.Quaternion(1, 0, 0, 0); 
         viewer.rotationGroup.quaternion.multiplySelf(viewer.dq);
         viewer.rotationGroup.quaternion.multiplySelf(viewer.cq);

        viewer.show();     
        
        var end = new Date();
        var testTime = end - start;
        
        resultTimes[rotateName] = testTime;
        console.log(timeMsg + (testTime) + "ms");        
        console.groupEnd();
             
        QUnit.ok(true, testMsg); 
        
    });
    
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
        
        Jmol.scriptWait(viewer, "select *; wireframe off; cartoon off; spacefill off; refresh;");
        
        console.group(testName);
        
        var start = new Date();
        
        Jmol.scriptWait(viewer, "select *; isosurface vdw; refresh;");
        
        var end = new Date();      
        
        var testTime = end - start;
        
        resultTimes[testName] = testTime;
        console.log(timeMsg + (testTime) + "ms");          
        console.groupEnd();
        
        QUnit.ok(true, testMsg); 
        
    });
    
    //Combo testcase
	var comboName = "combined";
    QUnit.test(comboName, function() {
        
        Jmol.scriptWait(viewer, "select *; wireframe off; cartoon off; spacefill off; isosurface off; refresh;");
        
        console.group(comboName);
        var start = new Date();

        
        Jmol.scriptWait(viewer, "select chain=A; spacefill; select chain=B; cartoon; color group; wireframe; select chain=C; wireframe 80;");
        Jmol.scriptWaitOutput(viewer,"refresh;");
        var end = new Date();      
        
        var testTime = end - start;
        
        resultTimes[comboName] = testTime;
        console.log(timeMsg + (testTime) + "ms");          
        console.groupEnd();
        
        QUnit.ok(true, testMsg); 
    });
    
    var rotateName = "rotate";
    QUnit.test(rotateName, function() {
        

        console.group(rotateName);
        var start = new Date();
        
        Jmol.scriptWaitOutput(viewer,"rotate 10; refresh;");
        var end = new Date();      
        
        var testTime = end - start;
        
        resultTimes[rotateName] = testTime;
        console.log(timeMsg + (testTime) + "ms");          
        console.groupEnd();
        
        QUnit.ok(true, testMsg); 
    });
    
}







