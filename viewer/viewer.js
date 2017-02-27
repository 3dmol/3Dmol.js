
var Viewer = function(){
	
	var entries=[];
	var selections=document.getElementById("selections");
	this.init=function(){

	};

	/*
	Adds a new selection object and initializes it 
	*/
	this.add_entry = function(type, parent){	

		var li = document.createElement('li');
		if(type == "Model"){
			var model = new Model();
			entries.push(model);

			var model_string = document.createElement('input');


			var label = document.createElement('p');
			label.className="label"
			label.innerHTML="M";
			li.appendChild(label);

			li.appendChild(model_string);
			li.className="entry"
			li.style.float="left";


		}else if(type == "Selection"){
			var selection = new Selection();
			entries.push(selection);

			var selection_string
		}else if(type == "Label"){

		}

		var selections=document.getElementById("selections");
		selections.appendChild(li);
		
		
	};

	this.read_url = function(){
		var url = window.location.search.substr(1);
	}

	this.update_url = function(){

	};



};

function openNav() {
    document.getElementById("sidenav").style.width = "250px";
    glviewer.translate(-125,0,500);
}

/* Set the width of the side navigation to 0 */
function closeNav() {
    document.getElementById("sidenav").style.width = "0";
    glviewer.translate(125,0,500);
}



var Model = function(){
	this.type="PDB";
	this.selections=[];
	this.model="";

	this.createModel = function(){
		console.log("model");
	}

	this.addChild= function(){

	}
}


var Selection = function(){
	this.type=null;
	this.query="";
	this.style;
	this.labels=[];

	this.addChild= function(){
		
	}
};

var Style = function(){

}

var Label = function(){
	this.text="";
	
};


var Menu = function(){
	this.items = [];
	this.actions = [];

}




var state = new Viewer();

state.init();