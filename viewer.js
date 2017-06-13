/*
builds an html tree that goes inside of the selection portion of the viewer page
*/
var buildHTMLTree = function(query){
    //get parent object for the html tree
    var parent = $('#selection_list');
    parent.text("");
    //list file type and path
    $("#model_type").attr("value",query.file.type);
    $("#model_input").attr("value",query.file.path);

    //loops through selections and creates a selection tree
    for(var selection_index in query.selections){
        var selection_object = query.selections[selection_index];

        //this function creates the selection object
        var createSelection = function(){
            //creates container
            var selection = $("<li/>",{
                class:"selection"
            });
            
            //add together sub selections
            // i think using object.keys is valid here
            var attribute_pairs =[];
            for(var subselection in selection_object.subselections){
                var obj=selection_object.subselections[subselection];
                attribute_pairs.push(Object.keys(obj)[0]+":"+obj[Object.keys(obj)[0]]);
            }
            var modifier=attribute_pairs.join(";");

            var selection_spec=$('<div/>', {
                class:'selection_spec',
                text:modifier,
                contenteditable:'true',
            }).appendTo(selection);        

            //creates style if it exists
            var model_specifications=$('<ul/>',{
                class:'model_specifications'
            });

            var createModelSpecification = function(model_spec_type,model_spec_object){
                var model_specification = null;
                var createStyle = function(spec){

                    var style=$('<li/>',{
                        text:"Style",
                        "class":"style",
                    });

                    var style_specs = $('<ul/>',{
                        "class":'style_specs',
                    }).appendTo(style);

                    var createStyleSpec = function(style_spec_object,style_spec_type,model_spec_type){
                        var style_spec=$('<li/>',{
                            "class":"style_spec",
                        });

                        var style_spec_name=$('<div/>',{
                            text:style_spec_type,
                            class:"style_spec_name",
                        }).appendTo(style_spec);

                        var style_spec_attributes = $('<ul/>',{
                            class:'style_spec_attributes',
                        }).appendTo(style_spec);

                        var createAttribute = function(name,value){

                            var attribute = $('<li/>',{
                                class:'attribute'
                            });
                            
                            var attribute_name = $('<span/>',{
                                class:'attribute_name',
                                text:name,
                                contenteditable:'true',
                            }).appendTo(attribute);

                            var attribute_value = $('<span/>',{
                                class:'attribute_value',
                                text:value,
                                contenteditable:'true',
                            }).appendTo(attribute);
                                                     
                            return attribute;
                        }

                        for(var attribute_index in style_spec_object){
                            createAttribute(attribute_index,style_spec_object[attribute_index]).appendTo(style_spec_attributes);
                        }

                        var add_attribute = $('<button/>',{
                            "class":"add_attribute",
                            "text":"Add Attribute",
                            "data-index":selection_index,
                            "data-type":model_spec_type,
                            "data-spectype":style_spec_type,
                            "click":function(){addAttribute(this)},
                        }).appendTo(style_spec);
                        
                        return style_spec;
                    }       
                    for(var attribute_index in model_spec_object.attributes){
                        createStyleSpec(model_spec_object.attributes[attribute_index],attribute_index,model_spec_type).appendTo(style_specs);
                    }

                    return style; 
                }

                //check for type
                if(model_spec_type=="style"){
                   model_specification = createStyle(model_spec_object.attributes)
                }

                var add_style_spec = $('<button/>',{
                    "class":"add_style_spec",
                    "text":"Add Style Spec",
                    "data-index":selection_index,
                    "data-type":model_spec_type,
                    "click":function(){addStyleSpec(this)},
                }).appendTo(model_specification);

                return model_specification;
            }

            //check if style exists and if so create the object
            if(selection_object.style !=undefined){
                var style = createModelSpecification("style",selection_object.style);
                //add style to model_specifications
                style.appendTo(model_specifications);
            }
            //add model_specifications to selection
            model_specifications.appendTo(selection);

            var add_model_spec = $('<button/>',{
                "class":"add_model_spec",
                "text":"Add Model Spec",
                "data-index":selection_index,
                "click":function(){addModelSpec(this)},
            }).appendTo(selection);

            return selection;
        }

        var selection=createSelection()
        // onClick="this.contentEditable='true';"
        selection.appendTo(parent);
        //creates a style tree
    }
}
/*
takes the query object and updates the url based on its contents
*/
var unpackQuery = function(query){
    var url= "";

    var unpackStyle = function(style){
        var style_attributes = style.attributes;
        var string="style=";

        var subStyles = [];
        $.each(style_attributes, function(key,value){

            var unpackSubStyle = function(val){
                var assignments = [];
                $.each(val, function(key, value){
                   assignments.push(key+"~"+value);
                });
                return assignments.join(",");
            };
            var sub_style_unpacked=unpackSubStyle(value);

            if(sub_style_unpacked=="")
                subStyles.push(key);
            else
                subStyles.push(key+":"+sub_style_unpacked);

        });
        string+=subStyles.join(";");
        return string;
    }

    var unpackSelection = function(selection){
        var subselections=selection.subselections;

        var subSelections = []
        $.each(subselections, function(index, value){
            $.each(value, function(key, value){
                subSelections.push(key+":"+value);
            });
        });

        var subselections_string = "select="+subSelections.join(";");

        var statements = [];
        statements=statements.concat(subselections_string);
        statements=statements.concat(unpackStyle(selection.style))
        statements=statements.concat(selection.dumps);
        
        string = statements.join("&");
        console.log(string)
        return string
    }
    //unpack file type and name
    url+=query.file.type+"="+query.file.path+"&";
    //unpack global style if it exists
    
    if(query.globalStyle!=null){
        url+=unpackStyle(query.globalStyle)+"&";
    }

    //unpack other selections and styles
    var unpacked =[]
    for(var sel in query.selections){
        unpacked.push(unpackSelection(query.selections[sel]));
    }
    url+=unpacked.join("&");
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
    glviewer.center({},1000,true);
}

var query = parseURL(window.location.search.substring(1));
console.log($3Dmol.specStringToObject(window.location.search.substring(1)))
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
    var width=400;
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

    
}

