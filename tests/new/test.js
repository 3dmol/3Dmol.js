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
	system[key](function(){
		var left_head=document.createElement('h4');
		left_head.innerHTML=key;
		var canvas=document.getElementById(key);
		var image=document.createElement('img');
		image.style.width="400px";
		image.style.height="400px";

		image.src=canvas.toDataURL('image/png');
		image.onclick=function(){
			var win = window.open();
				win.document.write(`<!DOCTYPE html>
<html>
  <head>
    <script src="../../build/3Dmol.js"></script>
  </head>

  <body>

  
          <div id="test10" style="width: 400px; height: 400px; position: relative;"></div>
        <script>
          var viewer = $3Dmol.createViewer($("#gldiv"));</script>`);
				win.document.write('<div id="gldiv" style="width: 100vw; height: 100vh"></div><script src="'+key+`.js"></script>
  </body>
</html>
`);
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
		console.log(after-before);
		if(i<keys.length-1 ){
			i+=1;
			runTest(i);
		}

	});

}


runTest(i);
//install a right-click handler on every canvas to export png
        
});
