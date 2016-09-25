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
/*
keys=[keys[0]];
function runTest(i){
	var key=keys[i];
	console.log("------------------------------");
	//first put the reference images up
	
	var img=document.createElement('img');
	img.style.width="400px";
	img.style.height="400px";
	img.src=imgs+"/"+key+".png";
	document.getElementById("right").appendChild(img);

	//then load the other image
	var image=document.createElement('img');
	console.log(key);
	try{
	system[key]();
	}catch(err){
		if(i<keys.length){
			i+=1;
			runTest(i);
		}
	}

}
runTest(i);

*/
system[keys[28]]();