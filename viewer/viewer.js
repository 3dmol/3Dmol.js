
var Viewer = function(){
	
	this.entries=[];
	this.selections=document.getElementById("selections");
	this.init=function(){

	};


	/*
	Adds a new selection object and initializes it 
	*/
	this.add_entry = function(){
		this.entries.push(new Model());

		//add model to the page
		var model=document.createElement("li");
		model.className="model";
		
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
	this.type=null;
	this.selections=[];
	this.model=null;
}

var PDBModel = function(){
	
}

var Selection = function(){
	this.type=null;
	this.query="";
	this.style;
	this.labels=[];


};

var Label = function(){
	this.text="";
	
};


var Menu = function(){
	this.items = [];
	this.actions = [];

}




var state = new Viewer();

state.init();