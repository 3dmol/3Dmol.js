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

	var tr=document.createElement('tr');
		tr.id=key;

	var column1=document.createElement('th');//label
		column1.className="label";
	var column2=document.createElement('th');//rendered image
		column2.className="rendered";

	var column3=document.createElement('th');//reference image
		column3.className="reference";

	var column4=document.createElement('th');//reference image
		column4.className="difference";

	var label=document.createElement('p');
	label.innerHTML=key;
	column1.appendChild(label);

	var reference=document.createElement('img');
	reference.src="imgs/"+key+".png";
	column3.appendChild(reference);

	tr.appendChild(column1);
	tr.appendChild(column2);
	tr.appendChild(column3);
	tr.appendChild(column4);

	return tr;
}

var h1=document.createElement('h3');
h1.innerHTML="Tests";
h1.style="display:inline";
document.getElementById("summary_scroll").appendChild(h1);
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

var par=document.createElement('p');
par.innerHTML="   0/"+keys.length;
document.getElementById("summary_scroll").appendChild(par);

par.style="display:inline";
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

var getUrl=function(){
var url=window.location.search.substring(1);
url=url+".html";
url=url.substring(0,url.indexOf(".html"));
if(url!==""){
document.documentElement.innerHTML = '';
	for(var key in keys){
		if(keys[key]===url){
			var div_test=document.createElement("div");
			div_test.id="div_"+keys[key];
			div_test.style="height:0px;width:0px;"
			document.body.appendChild(div_test);

			document.body.style="overflow:hidden";
			var div=document.createElement('div');
			div.id="gldiv";
			div.style="width: 100vw; height: 100vh; position: relative;";

			document.body.appendChild(div);

			var script=document.createElement('script');
			script.innerHTML+=system[keys[key]].toString();
			script.innerHTML="var callback=function(){};\nvar viewer=$3Dmol.createViewer($(\"#gldiv\"));\n"+script.innerHTML.substring(script.innerHTML.indexOf("try{")+4,script.innerHTML.indexOf("}catch(err)"));
			
			console.log(script.innerHTML);
			
			document.body.appendChild(script);
			
			}
		}

}else{

	runTest(i);   
}
}

function runTest(i){
	console.log("%c-------------------------- "+keys[i]+" -----------------------------",'background: green; color: white; display: block;')
	var before=Date.now();
	var key=keys[i];
	
	var viewer=$3Dmol.createViewer($("#gldiv"),{id:key});
	var afterGlobals;

	try{
	system[key](viewer,function(){
		waitfor(viewer.surfacesFinished , true , 100 , 0 , "" , function(){
			var after=Date.now();
			//gets the canvas
			var canvas=document.getElementById(key);
			//creates an image for the canvas
			var canvasImage=document.createElement("img");
				canvasImage.class="referenceImage";

			canvasImage.src=canvas.toDataURL("image/png");
			//click event for canvas
			canvasImage.onclick=function(){
				var win = window.open();
				win.location="generate_test.cgi?test="+key;
			};
			//create the table row
			var tableRow=createRow(key);
			tableRow.getElementsByClassName('rendered')[0].appendChild(canvasImage);
			document.getElementById("tests").getElementsByTagName('tbody')[0].appendChild(tableRow);
			var percentage=document.createElement('p');
			var differenceImage=document.createElement('img');
			var differ=0;
			var listElement=document.createElement('li');
			var anchor=document.createElement('a');
			anchor.href="#"+key;
			listElement.appendChild(anchor);
			anchor.innerHTML=key;
			par.innerHTML="   "+i+"/"+keys.length;

			resemble.outputSettings({
  				useCrossOrigin: false
			});

    		var diff = resemble(canvasImage.src).compareTo("imgs/"+key+".png").ignoreAntialiasing().scaleToSameSize().onComplete(function(data){
    		    //ignoreantialiasing provides some flex - scaletosamesize is necessary for retina displays
    			differ=data.rawMisMatchPercentage;//(100-blankDiff);
    			percentage.innerHTML=differ;
    			differenceImage.src=data.getImageDataUrl();
    			tableRow.getElementsByClassName('difference')[0].appendChild(differenceImage);
   				if(differ>5){
					tableRow.getElementsByClassName("label")[0].style.backgroundColor="red";
					document.getElementById('summary_scroll').appendChild(listElement)
					anchor.innerHTML+=" "+percentage.innerHTML;
				}else{
					tableRow.getElementsByClassName("label")[0].style.backgroundColor="green";
				}
   				
   	            tableRow.getElementsByClassName("label")[0].appendChild(percentage);
   	            //remove possible div
   	            $(".viewer_3Dmoljs").remove();
   	            //remove canvas
   	            $(canvas).remove();
   	            //compare globals before and after
   	            afterGlobals=GlobalTester.after(window);
   	            var str="";
   	            for(var field in afterGlobals){
   	                str+=field+"\n";
   	            }
   	            if(i<keys.length-1 ){
   	                //run the next test
   	                i+=1;
   	                runTest(i);
   	            }   				
    		}); //end onComplete

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
 
getUrl();
});
