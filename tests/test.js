$(document).ready(function(){
	console.log(system);

var getKeys = function(obj){
   var keys = [];
   for(var key in obj){
      keys.push(key);
   }
   return keys;
}

imgs="check/imgs";
var keys=getKeys(system)
var i=0;

function runTest(i){
	var key=keys[i];
	console.log("------------------------------");
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
	console.log(key);
	var div=document.createElement('div');
	div.id="div_"+key;
	document.getElementById("left").appendChild(div);
	system[key](function(){
		var left_head=document.createElement('h4');
		left_head.innerHTML=key;
		var canvas=document.getElementById(key);
		var image=document.createElement('img');
		document.getElementById("div_"+key).appendChild(left_head);
		document.getElementById("div_"+key).appendChild(image);
		console.log(canvas);
		image.src=canvas.toDataURL('image/png');
		$("#undefined").remove();
		$(canvas).remove();
		//$("#gldiv").children()[0].remove();
		if(i<keys.length-1){
			i+=1;
			console.log(i);
			runTest(i);
		}
	});


}

runTest(i);
//install a right-click handler on every canvas to export png
        
});
