    

var width=300;
var urlObject={};


var createSelection = function(selection_object){
    var li=document.createElement('li');

    li.className="selection";

    var delete_button=document.createElement('a');
    delete_button.className="delete_button";
    delete_button.innerHTML="&#x2715;";
    $(delete_button).click(function(event){
        removeSelection(selection);
    });
    delete_button.style.float="right";
    delete_button.style.color="white";
    var selection_div=document.createElement('div');

    selection_div.innerHTML=selection_object.selection;

    selection_div.className="selection_div";
    selection_div.appendChild(delete_button);

    li.appendChild(selection_div);
    //sub selections(list of the style,surface, or labelres)
    var div=document.createElement('div');
    div.className="subSelection";

    var sublist=document.createElement('ul');

    sublist.style.visibility="hidden";

    for(var i=0;i<selection_object.list.length;i++){
        var lst=createSubSelection(selection_object.list[i]);

        sublist.appendChild(lst);
    }

    div.appendChild(sublist);

    li.appendChild(div);

    //sub sub selections

    var list=document.createElement('ul');


    for(var i=0;i<selection_object.list.length;i++){
        var sublist_object=createSubSelection(selection_object.list[i]);
        list.appendChild(sublist_object);
    }

    list.style.visibility="hidden";

    li.appendChild(list);

    return li;
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

    console.log(query);
    return query;
}


var query= parseURL(window.location.href.substring(29));




//initializes the sidebar based on the given url
var initSide = function(url){

    console.log(query.file);
    //model type value
    document.getElementById("model_type").value=query.file.fileType;
    //query value
    document.getElementById("model_input").value=query.file.fileValue;

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

