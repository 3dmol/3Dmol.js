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
for(var key in keys){
	//first put the reference images up
	var img=document.createElement('img');
	img.src=imgs+"/"+keys[key]+".png";
	document.getElementById("right").appendChild(img);

	//then load the other image
	var image=document.createElement('img');

	system[keys[key]](function(){
		image.src=document.getElementById("gldiv").firstChild.toDataURL('image/png');
		document.getElementById("left").appendChild(image);
		console.log("finished");
	});

}

