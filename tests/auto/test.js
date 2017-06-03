$(document).ready(function(){
var getKeys = function(obj){
   var keys = [];
   for(var key in obj){
      keys.push(key);
   }
   return keys;
}

var canvasCount= function(){

    return $('#gldiv').children().length;
}

/**
Tests for global vairiables 
*/
var GlobalTester = (function(){
    var fields={}
    var before = function(w){
        for(var field in w){
            fields[field] = true;
        }
        return fields;
    }

    var after = function(w){
        var errors={};
        for(var field in w){
            if(!fields[field]){
                delete window[field];
                errors[field]=window[field];
            }            
        }
        return errors;

    }
    return {before: before, after:after};
}());

function waitfor(test, expectedValue, msec, count, source, callback) {
    // Check if condition met. If not, re-check later (msec).
    while (test() !== expectedValue) {
        count++;
        setTimeout(function() {
            waitfor(test, expectedValue, msec, count, source, callback);
        }, msec);
        return;
    }
    // Condition finally met. callback() can be executed.
    callback();
}

function createRow(key){

    var tr = $('<tr>').prop('id',key);
    $('<td class="label">').appendTo(tr).append($('<p>'+key+'</p>'));
    $('<td class="rendered">').appendTo(tr);
    $('<td class="reference">').appendTo(tr).append($('<img src="imgs/'+key+'.png">'));
    $('<td class="difference">').appendTo(tr);

    return tr;
}

var h1=$('<h3>Tests</h3>').css('display','inline');
$("#summary_scroll").append(h1);
imgs="imgs";
var keys=getKeys(system)
keys.sort();
var tests=[];
var copy=keys.slice(0);
for(var i=0;i<keys.length;i++){
    if(keys[i].substring(0,4)=="test"){
        tests.push(i);
        copy.splice(i);
    }
}
//sort the tests
for(var i=1;i<tests.length;i++){
    var j=i;
    while(j>0 && parseInt(keys[tests[j-1]].substring(4)) > parseInt(keys[tests[j]].substring(4))){
        var hold=tests[j];
        tests[j]=tests[j-1];
        tests[j-1]=hold;
        j--;
    }
}

var par=$('<p>   0/'+keys.length+'</p>').css('display','inline');
$("#summary_scroll").append(par);

var new_arr=[]
for(var test=0;test<tests.length;test++){
    new_arr.push(keys[tests[test]]);
}
keys=new_arr;
keys=keys.concat(copy);
var beforeGlobals;
var i=0;
$('#gldiv').hide();
beforeGlobals=GlobalTester.before(window);

// apparently toDataURL isn't technically a standard for webgl canvases
// and (some versions of) safari return the image flipped vertically
// this returns an image uri in a (hopefully) portable way
function imageFromWebGlCanvas(canvas) {
    var w = canvas.width;
    var h = canvas.height;
    var c = $("<canvas>").get(0);
    c.width = w;
    c.height = h;
    var drawctx = c.getContext("2d");
    var ctx3d = canvas.getContext("webgl");
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
}

function runTest(i){
    console.log("%c-------------------------- "+keys[i]+" -----------------------------",'background: green; color: white; display: block;')
    var before=Date.now();
    var key=keys[i];
    
    var viewer=$3Dmol.createViewer($("#gldiv"),{id:key});
    var afterGlobals;

    try{
        system[key](viewer,function(){
            waitfor(function() { return viewer.surfacesFinished() && !viewer.isAnimated() } , true , 250 , 0 , "" , function(){
                //create the table row
                var tableRow=createRow(key);
                var percentage=$('<p>');
                tableRow.find(".label").append(percentage);

                var setError = function(msg) {
                    var listElement=$('<li>').append($('<a href="#'+key+'">'+key+' '+msg+'</a>'));
                    tableRow.find('.label').css('backgroundColor','red');
                    $('#summary_scroll').append(listElement)
                    percentage.html(msg);
                }
                
                try {
                    var after=Date.now();
                    //gets the canvas
                    var canvas=document.getElementById(key);
                    //creates an image for the canvas
                    var canvasImageData = imageFromWebGlCanvas(canvas);
                    var canvasImage=$("<img class='referenceImage'>").attr('src',canvasImageData);

                    //click event for canvas
                    canvasImage.click(function(){
                        var win = window.open();
                        win.location="generate_test.cgi?test="+key;
                    });
                    tableRow.find('.rendered').append(canvasImage);
                    $('#tests').find('tbody').append(tableRow);

                    var differenceImage=$('<img>');
                    var differ=0;

                    par.html("   "+i+"/"+keys.length);

                    resemble.outputSettings({
                        useCrossOrigin: false
                    });

                    var diff = resemble(canvasImageData).compareTo("imgs/"+key+".png").ignoreAntialiasing().scaleToSameSize().onComplete(function(data){
                        //ignoreantialiasing provides some flex - scaletosamesize is necessary for retina displays
                        try {
                            differ=data.rawMisMatchPercentage;//(100-blankDiff);
                            percentage.html(differ+'<br>'+(after-before)+'ms');
                            differenceImage.attr('src',data.getImageDataUrl());
                            tableRow.find('.difference').append(differenceImage);
                            if(differ>5){
                                setError(differ);
                            }else{
                                tableRow.find(".label").css('backgroundColor',"green");
                            }
                            
                            //remove possible div
                            $(".viewer_3Dmoljs").remove();
                            //remove canvas
                            $(canvas).remove();
                            //compare globals before and after
                            afterGlobals=GlobalTester.after(window);

                        } catch(e) {
                            setError("Img Error "+e);
                        }
                        
                        i++;
                        if(i<keys.length){
                            //run the next test
                            runTest(i);
                        }                   
                    }); //end onComplete
                }catch(e) {
                    setError("Error "+e);
                }
            });
            
        });

    }catch(e){
        console.log("caught");
        console.log(e);
        console.log("caught");
        i+=1
        console.log(i);
        runTest(i);
    }
        
}    
 
runTest(70);
});
