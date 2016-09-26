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
	
	var right_head=document.createElement('h1');
	right_head.innerHTML=key;
	document.getElementById("right").appendChild(right_head);
	var img=document.createElement('img');
	img.style.width="400px";
	img.style.height="400px";
	img.src=imgs+"/"+key+".png";
	document.getElementById("right").appendChild(img);
	
	//then load the other image
	var image=document.createElement('img');
	console.log(key);
	var left_head=document.createElement('h1');
	left_head.innerHTML=key;
	document.getElementById("left").appendChild(left_head);
	system[key](function(){
		$('#gldiv').find('div').first().remove();
		if(i<keys.length){
			i+=1;
			runTest(i);
		}
	});


}
runTest(i);