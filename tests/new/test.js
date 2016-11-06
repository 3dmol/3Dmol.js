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
if(url!==""){

	document.body.innerHTML="";
	for(var key in keys){
		if(keys[key]===url){
			var script=document.createElement('script');
			document.body.write("<div id=\"gldiv\" style=\"width: 100vw; height: 100vh; position: relative;\"></div>");

			document.body.appendChild(script);
			}
		}

}else{

	runTest(i);   
}
}

function runTest(i){

	var before=Date.now();
	var key=keys[i];
	console.log("running test "+key+":");
	//first put the reference images up
	
	var right_head=document.createElement('h4');
	right_head.innerHTML=key;
	document.getElementById("right").appendChild(right_head);
	var img=document.createElement('img');
	img.style.width="400px";
	img.style.height="400px";
	img.src=imgs+"/"+key+".png";
	document.getElementById("right").appendChild(img);
	//then load the other image
	var div=document.createElement('div');
	div.id="div_"+key;
	document.getElementById("left").appendChild(div);
	
	
	var textarea=document.createElement('textarea');
	textarea.style.height="100px";
	textarea.style.width="80%"
	textarea.id="right_"+key;
	textarea.style.visibility="hidden";
	document.getElementById("right").appendChild(textarea);
	var texta=document.createElement('textarea');
	texta.style.height="100px";
	texta.style.width="80%"
	texta.id="left_"+key;
	texta.value="";
	document.getElementById("left").appendChild(texta);
	var viewer=$3Dmol.createViewer($("#gldiv"),{id:key});
	var afterGlobals;

	system[key](viewer,function(canvasid){
		waitfor(viewer.surfacesFinished,true,100,0,"",function(){
		var after=Date.now();
		
		var textarea=document.getElementById("left_"+key);
		var left_head=document.createElement("h4");
		left_head.innerHTML=key;
		var canvas=document.getElementById(key);
		var image=document.createElement("img");
		image.style.width="400px";
		image.style.height="400px";

		image.src=canvas.toDataURL("image/png");
		image.onclick=function(){
			if(key.substring(0,4)=="test"){
				var win = window.open();
				win.document.write("<!DOCTYPE html><html><head><script src=\"../../build/3Dmol.js\"></script></head><body style=\"overflow-y:hidden;overflow-x:hidden;\"><div id=\"gldiv\" style=\"width: 100vw; height: 100vh; position: relative;\"></div><script>var viewer = $3Dmol.createViewer($(\"#gldiv\"));</script>");
				win.document.write("<script src=\""+key+".js\"></script></body></html>");
		}else{
			var win = window.open();
			win.document.write("<!DOCTYPE html><html><head><script src=\"../../build/3Dmol.js\"></script></head><body style=\"overflow-y:hidden;overflow-x:hidden;\"><div id=\"gldiv\" style=\"width: 100vw; height: 100vh; position: relative;\"></div><script>var viewer = $3Dmol.createViewer($(\"#gldiv\"));</script>");
			win.document.write("<script>var sys={func:"+system[key].toString()+"};sys.func(viewer,function(){});</script>");
		}
		};
		document.getElementById("div_"+key).appendChild(left_head);
		document.getElementById("div_"+key).appendChild(image);
		//$("#undefined").remove();
		$(".viewer_3Dmoljs").remove();
			$(canvas).remove();
		//$("#gldiv").children()[0].remove();
		textarea.value+="\ntimestamp : "+(after-before)+" ms\n";

		afterGlobals=GlobalTester.after(window);
		var str="";
		for(var field in afterGlobals){
			str+=field+"\n";
		}
		if(afterGlobals!==undefined){
			
			textarea.value+="\n"+str;
			
		}

		if(i<keys.length-1 ){
			i+=1;
			runTest(i);
		}

		});
		
		});
		
}    
 
getUrl();
});
