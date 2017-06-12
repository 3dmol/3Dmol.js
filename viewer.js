var width=400;
var urlObject={};

/*
builds an html tree that goes inside of the selection portion of the viewer page
*/
var buildHTMLTree = function(query){
    var parent = document.getElementById('selection_list');
    parent.innerHTML="";
    //list file type and path
    document.getElementById("model_type").value=query.file.type;
    document.getElementById("model_input").value=query.file.path;

    //loops through selections and creates a selection tree
    var selections = query.selections;
    for(var i=0;i<selections.length; i++){

        var createSelection = function(){
            //creates container
            var selection = document.createElement('li');
            selection.className =  "selection";// ui-sortable-handle
            
            //add together sub selections
            var modifier = "";
            for(var j = 0; j<selections[i].subselections.length;j++){
                if(j!=0)
                    modifier+=";"
                modifier+=Object.keys(selections[i].subselections[j])[0]+":"+selections[i].subselections[j][Object.keys(selections[i].subselections[j])[0]]
            }
            
            var selection_spec = document.createElement('div');
            selection_spec.innerHTML=modifier
            var att= document.createAttribute('contenteditable');
            att.value='true';
            selection_spec.setAttributeNode(att);
            selection_spec.className='selection_spec'
            selection.appendChild(selection_spec)

            //creates style if it exists
            var model_specifications=document.createElement('ul');
            model_specifications.className='model_specifications';

            var createModelSpecification = function(type,spec){

                var model_specification = null;
                var createStyle = function(spec){

                    var style = document.createElement('li');
                    style.innerHTML='Style';
                    style.className='style';

                    var style_specs= document.createElement('ul');
                    style_specs.className='style_specs';
                    style.appendChild(style_specs)

                    var createStyleSpec= function(specification,typ,type){
                        var style_spec = document.createElement('li');
                        style_spec.className = 'style_spec';

                        var style_spec_name = document.createElement('div');
                        style_spec_name.innerHTML= typ;

                        style_spec.appendChild(style_spec_name);

                        var style_spec_attributes = document.createElement('ul');
                        style_spec_attributes.className="style_spec_attributes";

                        style_spec.appendChild(style_spec_attributes);

                        var createAttribute = function(name,value){
                            var attr = document.createElement('li');
                            attr.className = 'attribute';

                            var name_div = document.createElement('span');
                            name_div.className = "attribute_name";
                            name_div.innerHTML = name;
                            var att_name= document.createAttribute('contenteditable');
                            att.value='true';
                            name_div.setAttributeNode(att_name);

                            
                            var value_div = document.createElement('span');
                            value_div.className = "attribute_value";
                            value_div.innerHTML= value;

                            var att_value= document.createAttribute('contenteditable');
                            att_value.value='true';
                            value_div.setAttributeNode(att_value);

                            attr.appendChild(name_div);
                            attr.appendChild(value_div);


                            return attr;
                        }

                        var attrs = Object.keys(specification)

                        for(var attri=0; attri<attrs.length;attri++){
                            console.log(attrs[attri] +" "+ specification[attrs[attri]])
                            style_spec_attributes.appendChild(createAttribute(attrs[attri],specification[attrs[attri]]))
                        }

                        var add_attribute= document.createElement('button');
                        add_attribute.className="add_attribute";
                        add_attribute.innerHTML="Add Attribute"

                        var atr = document.createAttribute('obj');
                        console.log(i,type,typ)
                        atr.value=i+","+type+","+typ;
                        add_attribute.setAttributeNode(atr);

                        add_attribute.onclick= function(){addAttribute(this)};
                        style_spec.appendChild(add_attribute);
                        
                        return style_spec;
                    }

                    var specAttrKeys = Object.keys(spec.attributes);
                    for(var index = 0; index<specAttrKeys.length;index++){
                        //loop through style_specs and create style specs for them
                        var specif = spec.attributes[specAttrKeys[index]];
                        var typ = specAttrKeys[index]
                        style_specs.appendChild(createStyleSpec(specif,typ,type));                      

                    }                    
                    return style; 
                }

                //check for type
                if(type=="style"){
                   model_specification = createStyle(spec)
                }
                var add_style_spec= document.createElement('button');
                add_style_spec.className="add_style_spec";
                add_style_spec.innerHTML="Add Style Spec";

                var obj = document.createAttribute('obj');
                obj.value=[i,type];
                add_style_spec.setAttributeNode(obj);

                add_style_spec.onclick = function(){addStyleSpec(this)};

                model_specification.appendChild(add_style_spec);
                return model_specification;

            }

            //check if style exists and if so create it
            if(selections[i].style !=undefined){
                var style = createModelSpecification("style",selections[i].style);
                model_specifications.appendChild(style);
            }

            selection.appendChild(model_specifications);


            //this is where surfaces and other thigns will be created
            var add_model_spec = document.createElement('button');
            add_model_spec.className="add_model_spec";
            add_model_spec.innerHTML="Add Model Spec";
            var obj = document.createAttribute('obj');
            obj.value=i;
            add_model_spec.setAttributeNode(obj);
            add_model_spec.onclick=function(){addModelSpec(this)};
            selection.appendChild(add_model_spec);
            return selection;
        }

        var selection=createSelection()
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

    var unpackStyle = function(style){
        var attr = style.attributes;
        var str="style=";
        for(var attribute in attr){
            var index = Object.keys(attr).indexOf(attribute)
            if(index!=0)
                str+=";";
            str+=attribute;
            if(Object.keys(attr[attribute])!=undefined){
                if(Object.keys(attr[attribute]).length!=0 )
                    str+=":";
                for(var style in attr[attribute]){
                    if(Object.keys(attr[attribute]).indexOf(style))
                        str+=",";
                    str+=style+"~"+attr[attribute][style];
                }
            }
        }
        return str;
    }

    var unpackSelection = function(selection){
        var sels=selection.subselections;
        var str = "select=";
        for(var sel in sels){
            if(Object.keys(sels).indexOf(sel)!=0)
                str+=";";
            if(Object.keys(sels[sel])!=undefined){
                for(var style in sels[sel]){

                    str+=style+":"+sels[sel][style];
                }
            }
        }
        str+="&"+unpackStyle(selection.style);
        for(var i in selection.dumps){
            if(Object.keys(sels).length!=0)
                str+="&"
            
            str+=selection.dumps[i]
        }
        return str;
    }
    //unpack file type and name
    url+=query.file.type+"="+query.file.path+"&";
    //unpack global style if it exists
    
    if(query.globalStyle!=null){
        url+=unpackStyle(query.globalStyle)+"&";
    }

    //unpack other selections and styles
    var selections= query.selections;
    for(var sel in selections){
        url+=unpackSelection(selections[sel])+"&";
    }
    console.log(url)
    return url;
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
    this.dumps=[];
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
            }
        }else if(strType(tokens[i])=="select"){
            var split=tokens[i].split("=")
            currentSelection=new Selection(split[1])
            query.selections.push(currentSelection)
        }else{
            currentSelection.dumps.push(tokens[i]);
        }

        console.log(tokens[i])
    }
    console.log(query);
    unpackQuery(query);
    return query;
}

//these functions all edit the query object 
var addSelection = function(){
    query.selections.push(new Selection(""))
    buildHTMLTree(query);
}
var addModelSpec = function(selection){
    if(query.selections[selection.getAttribute("obj")].style==null)
        query.selections[selection.getAttribute("obj")].style=new Style("");
    buildHTMLTree(query);
}

var addStyleSpec = function(model_spec){
    var str=model_spec.getAttribute("obj")
    var i=str.split(",")[0];
    var type=str.split(",")[1];
    console.log(query.selections[i][type])
    query.selections[i][type].attributes[query.selections[i][type].attributes.length]="";
    buildHTMLTree(query);
}

var addAttribute = function(style_spec){
    var list=style_spec.getAttribute("obj").split(",");
    console.log(query.selections[list[0]][list[1]].attributes[list[2]])
    query.selections[list[0]][list[1]].attributes[list[2]][""]="";
    buildHTMLTree(query);
}
//this function reads the form changes and upates the query accordingly
var updateQuery = function(){
    //console.log(document.getElementById("model_input").value)
    query.file.path=document.getElementById("model_input").value;
    query.file.type=document.getElementById("model_type").value;
}


var center = function(){
    glviewer.center({},1000,false);
}

var query = parseURL(window.location.search.substring(1));
//this function compresses the html object back into a url
var render = function(){
    //calls update query
    updateQuery();
    setURL(unpackQuery(query));
}
//initializes the sidebar based on the given url
var initSide = function(url){
    var list = document.createElement('ul')
    document.getElementById('container').appendChild(list);
    buildHTMLTree(query);
    //model type value
    
}

//opens up the side bar
var openSide= function(){
    console.log("open")
    document.getElementById("sidenav").style.width = width+"px";
    document.getElementById("menu").style.visibility="hidden";
    //document.getElementById("url").value=window.location.href.substring(window.location.href.indexOf("?")+1);
    
    glviewer.zoomTo();
    glviewer.render();
}
//closes the side bar
var closeSide= function(){

    document.getElementById("menu").style.visibility="visible";
    document.getElementById("sidenav").style.width = "0";

    glviewer.translate(-width/2,0);
}