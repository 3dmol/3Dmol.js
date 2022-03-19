/* 
 * QUnit benchmark tests for $3Dmol, GLmol, and JSmol
 */
// Test rendering performance for different sized pdb's 

const profile = QUnit.urlParams.profilecheck;

const testSuite = testSuite || "3Dmol";

// Style types to test
const styleSpec = ["line", "stick", "sphere", "cartoon"];

// $3Dmol testcase generator

const gen$3DmolTestCase = function(styleType) {
    
    const testName = `${styleType  } render`;
    const timeMsg = `${styleType  } render time: `;
    const testMsg = `${styleType  } style set correctly`;
    const style = {}; style[styleType] = {};
    
    QUnit.test(testName, () => {
        
        viewer.setStyle({}, {cross:{}});
        viewer.render();
        console.group(testName);          
        
        const start = new Date();
        
        viewer.setStyle({}, style);           
        viewer.render();
        
        const end = new Date();
        const testTime = end - start;
        
        viewer.rotate(10);
        const time2 = new Date() - end;
        
        resultTimes[testName] = testTime;
        console.log(`${timeMsg + (testTime)  }ms`);
        console.log(`${timeMsg + (testTime)  }ms; rotate ${time2}ms`);          
        console.groupEnd();
        
        QUnit.ok(true, testMsg);
        
    });
};

// GLmol test generator

const genGLmolTestCase = function(styleType, profile) {
    
    const testName = `${styleType  } render`;
    const timeMsg = `${styleType  } render time: `;
    const testMsg = `${styleType  } style set correctly`;
    
    const defineRep = function(style){
        const all = viewer.getAllAtoms();
        
        return function() {
            // var all = this.getAllAtoms();
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
        
        const view = viewer.getView();
        viewer.initializeScene();
        const all = viewer.getAllAtoms();
        viewer.zoomInto(all);
        viewer.defineRepresentation = defineRep(styleType);
        
        const start = new Date();
        // Draw appropriate style        
        viewer.defineRepresentation();
        viewer.setView(view);        
        viewer.show();     
        
        const end = new Date();
        const testTime = end - start;
        
        resultTimes[testName] = testTime;
        console.log(`${timeMsg + (testTime)  }ms`);        
        console.groupEnd();
             
        ok(true, testMsg); 
        
    });
    

};

// JSmol testcase generator

const genJSmolTestCase = function(styleType, profile) {
    
    const testName = `${styleType  } render`;
    const timeMsg = `${styleType  } render time: `;
    const testMsg = `${styleType  } style set correctly`;
    let script = "";
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
      
    // Create test case
    QUnit.test(testName, () => {
        
        console.group(testName);
        const start = new Date();
        Jmol.scriptWait(viewer, script);
        
        Jmol.scriptWaitOutput(viewer, "refresh;");
        
        const end = new Date();      
        
        const testTime = end - start;
        
        Jmol.scriptWaitOutput(viewer,"rotate 10; refresh;");
        
        const time2 = new Date() - end;
        resultTimes[testName] = testTime;
        console.log(`${timeMsg + (testTime)  }ms; rotate ${time2}ms`);          
        console.groupEnd();
        
        QUnit.ok(true, testMsg); 
        
    });
};


// Setup test modules


// $3Dmol tests
if (testSuite === '3Dmol') {
    QUnit.module( "$3Dmol Tests", {

        setupOnce() {
            viewer.zoomTo();
            console.group("$3Dmol");
        },

        teardownOnce() {
            console.groupEnd();
        }

    });

    // $3Dmol test cases 
    for (const style in styleSpec)
        gen$3DmolTestCase(styleSpec[style], profile);    
    // surface test
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
        
    // Combo testcase
	const comboName = "combined";
    QUnit.test(comboName, () => {
        
        viewer.setStyle({}, {});
        viewer.removeAllSurfaces();
        viewer.render();        
        console.group(comboName);
        const start = new Date();

        viewer.setStyle({chain: 'A'}, {sphere: {} });
        viewer.setStyle({chain: 'B'}, {cartoon: {color: 'spectrum'} });
        viewer.setStyle({chain: 'C'}, {stick: {} });
        viewer.render();
        // Jmol.scriptWait(viewer, "select chain=A; spacefill; select chain=B; cartoon; color group; wireframe; select chain=C; wireframe 80;");
        const end = new Date();      
        
        const testTime = end - start;
        
        resultTimes[comboName] = testTime;
        console.log(`${timeMsg + (testTime)  }ms`);          
        console.groupEnd();
        
        QUnit.ok(true, testMsg); 
    });
    
    const rotateName = "rotate";
    QUnit.test(rotateName, () => {
        

        console.group(rotateName);
        const start = new Date();
        
        viewer.rotate(10);

        const end = new Date();      
        
        const testTime = end - start;
        
        resultTimes[rotateName] = testTime;
        console.log(`${timeMsg + (testTime)  }ms`);          
        console.groupEnd();
        
        QUnit.ok(true, testMsg); 
    });

}


else if (testSuite === "glmol") {
// GLmol testing module

    QUnit.module( "GLmol Tests", {

        setupOnce() {
            const all = viewer.getAllAtoms();
            viewer.zoomInto(all);
            console.group("GLmol");
        },

        teardownOnce() {
            console.groupEnd();
            // glviewer.removeAllModels();
        }
    });

//    for (var style in styleSpec)
//        genGLmolTestCase(styleSpec[style], profile);  
    const combinedName = "combined";
    QUnit.test(combinedName, () => {
        
        console.group(combinedName);
        
        const view = viewer.getView();
        viewer.initializeScene();
        const all = viewer.getAllAtoms();
		const chainA = [];
		const chainB = [];
		const chainC = [];
		
		for(let i = 0; i < all.length; i++) {
			if(viewer.atoms[all[i]].chain === 'A') chainA.push(all[i]);
			else if(viewer.atoms[all[i]].chain === 'B') chainB.push(all[i]);
			else if(viewer.atoms[all[i]].chain === 'C') chainC.push(all[i]);
		}
        
        const start = new Date();
        // Draw appropriate style        
          //      viewer.setStyle({chain: 'A'}, {sphere: {} });
        // viewer.setStyle({chain: 'B'}, {cartoon: {color: 'spectrum'} });
        viewer.drawAtomsAsSphere(viewer.modelGroup, chainA, viewer.sphereRadius);
        viewer.colorChainbow(chainB);
        viewer.drawCartoon(viewer.modelGroup, chainB, false, viewer.thickness);
        viewer.drawBondsAsStick(viewer.modelGroup, chainC, viewer.cylinderRadius, viewer.cylinderRadius, true);
        
        viewer.zoomInto(all);
        viewer.show();     
        
        const end = new Date();
        const testTime = end - start;
        
        resultTimes[combinedName] = testTime;
        console.log(`${timeMsg + (testTime)  }ms`);        
        console.groupEnd();
             
        QUnit.ok(true, testMsg); 
        
    });
    
        const rotateName = "rotate";
    QUnit.test(rotateName, () => {
        
        console.group(rotateName);
        
        const view = viewer.getView();
        const dx = .1;
      const dy = 0;
      const r = Math.sqrt(dx * dx + dy * dy);

           const rs = Math.sin(r * Math.PI) / r;


        const start = new Date();
        // Draw appropriate style        
          //      viewer.setStyle({chain: 'A'}, {sphere: {} });
        // viewer.setStyle({chain: 'B'}, {cartoon: {color: 'spectrum'} });
                 viewer.dq.x = Math.cos(r * Math.PI); 
         viewer.dq.y = 0;
         viewer.dq.z =  rs * dx; 
         viewer.dq.w =  rs * dy;
         viewer.rotationGroup.quaternion = new THREE.Quaternion(1, 0, 0, 0); 
         viewer.rotationGroup.quaternion.multiplySelf(viewer.dq);
         viewer.rotationGroup.quaternion.multiplySelf(viewer.cq);

        viewer.show();     
        
        const end = new Date();
        const testTime = end - start;
        
        resultTimes[rotateName] = testTime;
        console.log(`${timeMsg + (testTime)  }ms`);        
        console.groupEnd();
             
        QUnit.ok(true, testMsg); 
        
    });
    
}

else if (testSuite === "glmol_surf") {

    QUnit.module( "GLmol Surface Test", {

        setupOnce() {
            const all = viewer.getAllAtoms();
            viewer.zoomInto(all);
            console.group("GLmol");
        },

        teardownOnce() {
            console.groupEnd();
        }
    });


    const testName = "surf render";
    const timeMsg = "surf render time: ";
    const testMsg = "surf style set correctly";
    
    const defineRep = function(){
        const all = viewer.getAllAtoms();
        
        return function() {
            // var all = this.getAllAtoms();
            const target = this.modelGroup;          
            this.drawAsCross(target, all, 0.3, true);

            
        }; 
        
    };
    
    QUnit.test(testName, () => {
        
        console.group(testName);
        
        const view = viewer.getView();
        viewer.initializeScene();
        const all = viewer.getAllAtoms();
        viewer.zoomInto(all);
        viewer.defineRepresentation = defineRep();
        
        
        // Draw appropriate style        
        // viewer.defineRepresentation();
        // viewer.drawAsCross(viewer.modelGroup, all, 0.3, true);
        const start = new Date();
        viewer.generateMesh(viewer.modelGroup, all, 1, false);
        
        viewer.setView(view);        
        viewer.show();     
        
        const end = new Date();
        const testTime = end - start;
        
        resultTimes[testName] = testTime;
        console.log(`${timeMsg + (testTime)  }ms`);        
        console.groupEnd();
             
        ok(true, testMsg); 
        
    });

}

else if (testSuite === "jmol") {
    
    // JSmol Tests
    QUnit.module( "Jmol Tests", {

        setupOnce() {   
            Jmol.scriptWait(viewer, "load test_structs/3M8L.pdb");
            const end = new Date();
            
            resultTimes.initialization = end - start;
            console.group("JSmol");  
            console.log("starting JSmol tests");
        },

        teardownOnce() {
            console.groupEnd();
            // $("#viewerdiv").empty();
        }

    });
/*
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
 */   
    // Combo testcase
	const comboName = "combined";
    QUnit.test(comboName, () => {
        
        Jmol.scriptWait(viewer, "select *; wireframe off; cartoon off; spacefill off; isosurface off; refresh;");
        
        console.group(comboName);
        const start = new Date();

        
        Jmol.scriptWait(viewer, "select chain=A; spacefill; select chain=B; cartoon; color group; wireframe; select chain=C; wireframe 80;");
        Jmol.scriptWaitOutput(viewer,"refresh;");
        const end = new Date();      
        
        const testTime = end - start;
        
        resultTimes[comboName] = testTime;
        console.log(`${timeMsg + (testTime)  }ms`);          
        console.groupEnd();
        
        QUnit.ok(true, testMsg); 
    });
    
    const rotateName = "rotate";
    QUnit.test(rotateName, () => {
        

        console.group(rotateName);
        const start = new Date();
        
        Jmol.scriptWaitOutput(viewer,"rotate 10; refresh;");
        const end = new Date();      
        
        const testTime = end - start;
        
        resultTimes[rotateName] = testTime;
        console.log(`${timeMsg + (testTime)  }ms`);          
        console.groupEnd();
        
        QUnit.ok(true, testMsg); 
    });
    
}







