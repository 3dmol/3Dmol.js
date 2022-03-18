$(document).ready(()=> {
    window.devicePixelRatio = 2.001; // for compat with headless
    const getKeys = function(obj){
       const keys = [];
       for(const key in obj){
          keys.push(key);
       }
       return keys;
    }

    const canvasCount= function(){

        return $('#gldiv').children().length;
    }

    /**
    Tests for global vairiables 
    */
    const GlobalTester = (function(){
        const fields={}
        const before = function(w){
            for(const field in w){
                fields[field] = true;
            }
            return fields;
        }

        const after = function(w){
            const errors=[];
            for(const field in w){
                if(!fields[field]){
                    delete window[field];
                    errors.push(field);
                }            
            }
            return errors;

        }
        return {before, after};
    }());

    function waitfor(test, expectedValue, msec, count, source, callback) {
        // Check if condition met. If not, re-check later (msec).
        if (test() !== expectedValue) {
            count+=1;
            setTimeout(() => {
                waitfor(test, expectedValue, msec, count, source, callback);
            }, msec);
            return;
        }
        // Condition finally met. callback() can be executed.
        callback();
    }

    function createRow(key){

        const tr = $('<tr class="examplerow">').prop('id',`${key}_row`);
        $('<td class="label">').appendTo(tr).append(`<a class="testname">${key}</a>`);
        $('<td class="rendered">').appendTo(tr);
        $('<td class="reference">').appendTo(tr).append(`<img src="../glcheck/reference-images/tests-glcheck-render-tests-${key}.html.png">`);
        $('<td class="difference">').appendTo(tr);
        $('#tests_table').append(tr);
        return tr;
    }

    $('#sorttable').click(() => {
        // sort the table with the worse examples first
        // detach rows from tbody
        const rows = $('#tests_table tr').detach().get();
        // sort
        rows.sort((a, b) => {

            let A = $(a).find('.percentage').html();
            let B = $(b).find('.percentage').html();
            A = parseFloat(A);
            B = parseFloat(B);
            if(isNaN(A)) A = 10000; // error messages go to front
            if(isNaN(B)) B = 10000;
            return B-A; // large to small

        });

        // reattach
        $('#tests_table').append(rows);
       
    });
    
    const keys=getKeys(system)
    keys.sort(); // lexographic ordering

    const par=$(`<p>   0/${keys.length}</p>`).css('display','inline');
    $("#summary_list").append(par);

    let beforeGlobals;
    const i=0;
    $('#gldiv').hide();


    function imageFromWebGlCanvas(canvas) {
        return canvas.toDataURL("image/png");
        /* the below has trouble with transparent backgrounds, it was needed
         * for older versions of safari.
        var w = canvas.width;
        var h = canvas.height;
        var c = $("<canvas>").get(0);
        c.width = w;
        c.height = h;
        var drawctx = c.getContext("2d");
        var ctx3d = canvas.getContext("webgl2");
        if (!ctx3d) ctx3d = canvas.getContext("webgl");
        const webglPixels = new Uint8Array(4 * w * h);

        //in what was probably a horrible design decision,
        //the ordering of pixels extracted from a webgl context is
        //flipped vertically from that of a 2d
        ctx3d.readPixels(0,0,w,h,ctx3d.RGBA,ctx3d.UNSIGNED_BYTE,webglPixels);

        var idata = drawctx.getImageData(0, 0, w, h);

        idata.data.set(webglPixels);
        //so we flip
        drawctx.putImageData(idata, 0, 0);  
        drawctx.scale(1, -1);
        drawctx.drawImage(c, 0, -h);
        return c.toDataURL("image/png"); // standard in 2d
        */
    }

    const beginTime = Date.now();
    let failures = 0;
    resemble.outputSettings({
        useCrossOrigin: false
    });
    
    // https://stackoverflow.com/questions/19327749/javascript-blob-filename-without-link
    function saveFile(blob, filename) {
      if (window.navigator.msSaveOrOpenBlob) {
        window.navigator.msSaveOrOpenBlob(blob, filename);
      } else {
        const a = document.createElement('a');
        document.body.appendChild(a);
        const url = window.URL.createObjectURL(blob);
        a.href = url;
        a.download = filename;
        a.click();
        setTimeout(() => {
          window.URL.revokeObjectURL(url);
          document.body.removeChild(a);
        }, 0)
      }
    }
    
    function runTest(i){
        if(i === keys.length) {
            // finished
            const endTime = Date.now();
            $('#summary_list').append(`<li class="totaltime">Total time: ${(endTime-beginTime)/1000}s</p>`);
            $('#summary_list').append(`<li class="failures">Total failures: ${failures}</p>`);
            
            zip.generateAsync({type:'blob'})
                .then((content) => {
                    $('#download').prop('disabled',false);
                    $('#download').click(() => {
                        saveFile(content,'images.zip');
                    });
                });
            return;
        }
        console.log(`%c-------------------------- ${i}: ${keys[i]} -----------------------------`,'background: green; color: white; display: block;')
        const before=Date.now();
        const key=keys[i];
        
        const viewer=$3Dmol.createViewer($("#gldiv"),{id:key});
        let afterGlobals;
        const nexti = i+1; // just in case exception happens after i increment

        // create the table row
        const tableRow=createRow(key);
        const percentage=$('<p class="percentage">');
        tableRow.find(".label").append(percentage);
        $('#tests').find('tbody').append(tableRow);
        
        par.html(`   ${i+1}/${keys.length}`);

        // setup error handling
        const setError = function(msg) {
            const listElement=$('<li class="erroritem">').append(`<a href="#${key}_row">${key} ${msg}</a>`);
            tableRow.find('.label').css('backgroundColor','red');
            $('#summary_list').append(listElement)
            percentage.html(msg);
            failures+=1;
        }
        
        // override default error handling (try/catch isn't sufficient due to callbacks)
        window.onerror = function(message, file, lineNumber) {
            setError(message);
            failures+=1;
            setTimeout(() => {runTest(nexti);}, 1); // let page update  
            return true; 
          };
          
        system[key](viewer,(viewer)=> {
            waitfor(() => viewer.surfacesFinished() && !viewer.isAnimated() , true , 1000 , 0 , "" , ()=> {
                const after=Date.now();
                // gets the canvas
                const canvas=viewer.getCanvas(); // $("canvas#"+key).get(0);
                // creates an image for the canvas
                const canvasImageData = imageFromWebGlCanvas(canvas);
                zip.file(`${key}.png`,canvasImageData.split('base64,')[1], {base64: true});
                const canvasImage=$(`<img name=${key} class='renderedImage'>`).attr('src',canvasImageData);

                // click event for canvas
                canvasImage.click(()=> {
                    const win = window.open();
                    win.location=`generate_test.cgi?test=${key}`;
                });
                tableRow.find('.rendered').append(canvasImage);
                tableRow.find('a.testname').attr('href',canvasImageData);
                tableRow.find('a.testname').attr('download',`${key}.png`);

                const differenceImage=$('<img>');
                let differ=0;

                const diff = resemble(canvasImageData).compareTo(`../glcheck/reference-images/tests-glcheck-render-tests-${key}.html.png`).set3DmolTolerances().scaleToSameSize().onComplete((data)=> {
                    // scaletosamesize is necessary for retina displays
                    if(data.error) {
                        throw data.error; // probably need to create an image file
                    }
                    differ=data.rawMisMatchPercentage;// (100-blankDiff);
                    percentage.html(`${differ}<br>${after-before}ms`);
                    differenceImage.attr('src',data.getImageDataUrl());
                    tableRow.find('.difference').append(differenceImage);
                    
                    // compare globals before and after
                    afterGlobals=GlobalTester.after(window);
                    
                    if(afterGlobals.length > 0) {
                        setError(`Globals added: ${afterGlobals}`);
                    } else if(differ>5){
                        setError(differ);
                    } else{
                        tableRow.find(".label").css('backgroundColor',"green");
                    }
                    
                    // remove possible div
                    $(".viewer_3Dmoljs").remove();
                    // remove canvas
                    $(canvas).remove();
                    
                    setTimeout(() => {runTest(nexti);}, 1); // let page update             
                }); // end onComplete
            });
            
        });
            
    }    
    
    // initialize a viewer since jquery adds some event handling stuff to the window
    // that we don't want to caught by the global tester
    $3Dmol.createViewer($("#gldiv"));    
    zip = new JSZip();
    GlobalTester.before(window);     
    
    runTest(0);
});
