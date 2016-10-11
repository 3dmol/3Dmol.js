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
    var fields = {};
    var before = function(w){
        for(var field in w){
            fields[field] = true;
        };
    };

    var after = function(w){
        for(var field in w){
            if(!fields[field]){
                 console.log("%c "+field + " has been added","color:#00cc00");
            }            
        };

    };
    return {before: before, after:after};
}());

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
	var par=document.createElement('p');
	par.innerHTML="placeholder";
	par.style.visibility="hidden";
	document.getElementById("right").appendChild(par);
	var viewer=$3Dmol.createViewer($("#gldiv"),{id:key});
	viewer.render_callback=function(){
		var left_head=document.createElement('h4');
		left_head.innerHTML=key;
		var canvas=document.getElementById(key);
		var image=document.createElement('img');
		image.style.width="400px";
		image.style.height="400px";

		image.src=canvas.toDataURL('image/png');
		image.onclick=function(){
			if(key.substring(0,4)=="test"){
				var win = window.open();
				win.document.write(`<!DOCTYPE html><html><head><script src="../../build/3Dmol.js"></script></head><body><div id="gldiv" style="width: 100vw; height: 100vh; position: relative;"></div><script>var viewer = $3Dmol.createViewer($("#gldiv"));</script>`);
				win.document.write('<script src="'+key+`.js"></script></body></html>`);
		}else{
			var win = window.open();
			win.document.write(`<!DOCTYPE html><html><head><script src="../../build/3Dmol.js"></script></head><body><div id="gldiv" style="width: 100vw; height: 100vh; position: relative;"></div><script>var viewer = $3Dmol.createViewer($("#gldiv"));</script>`);
			win.document.write("<script>var sys={func:"+system[key].toString()+"};sys.func();</script>");
		}
		};
		document.getElementById("div_"+key).appendChild(left_head);
		document.getElementById("div_"+key).appendChild(image);
		$("#undefined").remove();
		$(canvas).remove();
		//$("#gldiv").children()[0].remove();
		var after=Date.now();
		var p=document.createElement('p');
		p.innerHTML="timestamp : "+(after-before)+" ms";
		document.getElementById("div_"+key).appendChild(p);
		if(i<keys.length-1 ){
			i+=1;
			afterGlobals=GlobalTester.after(window);
			runTest(i);
		}
	};
	beforeGlobals=GlobalTester.before(window);
	var afterGlobals;
	
	system[key](viewer);
}
runTest(i);        
});
