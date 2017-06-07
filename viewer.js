    

var width=300;
var urlObject={};

/*
builds an html tree that goes inside of the selection portion of the viewer page
*/
var buildHTMLTree = function(query){
    var parent = document.getElementById('selection_list');

    //list file type and path
    document.getElementById("model_type").value=query.file.type;
    document.getElementById("model_input").value=query.file.path;

    //loops through selections and creates a selection tree
    var selections = query.selections;
    for(var i=0;i<selections.length; i++){
        var selection = document.createElement('li');
        selection.className =  "selection ui-sortable-handle";
        
        //add together sub selections
        var modifier = "";
        for(var j = 0; j<selections[i].subselections.length;j++){
            if(j!=0)
                modifier+=";"
            var subselections_keys=Object.keys(selections[i].subselections[j]);
            var subselections_values=Object.values(selections[i].subselections[j]);

            modifier+=subselections_keys[j]+":"+subselections_values[j]
        }
        var selection_div = document.createElement('div');
        selection_div.onClick="this.contentEditable='true';"
        selection_div.innerHTML=modifier
        selection_div.className='selection_div'
        selection.appendChild(selection_div)
        // onClick="this.contentEditable='true';"
        parent.appendChild(selection);
        //creates a style tree

    }
}

/*
takes the query object and updates the url based on its contents
*/
var unpackQuery = function(query){
    var url= "";
}

var updateHTMLTree = function(){

}

//style object for styling the selected model portion
function Style(modifier){
    this.attributes={};
    var sc_split=modifier.split(";")
    for(var i =0; i<sc_split.length;i++){//types such as line,cartoon
        var colon_split=sc_split[i].split(":")
        this.attributes[colon_split[0]]={}
        if(colon_split[1]!=undefined){
            var comma_split=colon_split[1].split(",")
            for(var j=0;j<comma_split.length;j++){
                var equals_split=comma_split[j].split("~");
                this.attributes[colon_split[0]][equals_split[0]]=equals_split[1]
            }
        }
    }
}
//selection object for displaying
function Selection(selection){
    this.subselections=[];
    var sc_split=selection.split(";")
    for(var i =0; i<sc_split.length;i++){//types such as line,cartoon
        var colon_split=sc_split[i].split(":")
        var obj={}
        obj[colon_split[0]]=colon_split[1]

        this.subselections.push(obj)
    }

    this.style=null;
}

function File(type,path){
    var fileTypes = {
        PDB:0,
        CID:1,
        URL:2,
    };

    this.path=path;
    this.type=type;
}

var Query = function(){
    this.globalStyle = null;

    this.selections = [];

    this.file = new File();

    //update function for query object
    this.update = function(){
        //this is called on every click of render
        //this should read the forms and change query according to the changes
    }

}


function setURL(urlPath){
    window.history.pushState({"html":"test","pageTitle":"test"},"", "viewer.html?"+urlPath);
}


var parseURL = function(url){
    var query = new Query();
    var tokens=url.split("&");
    console.log(tokens)

    function strType(str){
        if(str.indexOf("select")==0)
            return "select"
        else if(str.indexOf("pdb=")==0 || str.indexOf("cid=")==0 || str.indexOf("url=")==0)
            return "file"
        else if(str.indexOf("style")==0)
            return "style"
        return null;
    }
    var currentSelection=null;
    for(var i=0;i<tokens.length;i++){
        if(strType(tokens[i])=="file"){
            var split=tokens[i].split("=")
            query.file=new File(split[0],split[1])
        }else if(strType(tokens[i])=="style"){
            var split=tokens[i].split("=")
            var style=new Style(split[1])
            if(currentSelection==null){
                query.globalStyle=style;
            }else{
                currentSelection.style=style;
                currentSelection=null
            }

        }else if(strType(tokens[i])=="select"){
            var split=tokens[i].split("=")
            currentSelection=new Selection(split[1])
            query.selections.push(currentSelection)
        }
    }
    console.log(query)
    return query;
}


var query = parseURL(window.location.search.substring(1));

//initializes the sidebar based on the given url
var initSide = function(url){
    buildHTMLTree(query);
    //model type value

    var list=document.getElementById("selection_list");

    for(var i=0;i<query.selections.length;i++){
        list.appendChild(createSelection(query.selections[i]));
    }
}


//opens up the side bar
var openSide= function(){
    console.log("open")
    document.getElementById("sidenav").style.width = width+"px";
    document.getElementById("menu").style.visibility="hidden";
    //document.getElementById("url").value=window.location.href.substring(window.location.href.indexOf("?")+1);
    glviewer.translate(width/2,0);
    glviewer.zoomTo();
    glviewer.render();
}
//closes the side bar
var closeSide= function(){

    document.getElementById("menu").style.visibility="visible";
    document.getElementById("sidenav").style.width = "0";

    glviewer.translate(-width/2,0);
}

