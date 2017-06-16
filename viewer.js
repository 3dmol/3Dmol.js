//removes style,labelres, and surface from a copy of the selection object and returns it
var augmentSelection = function(selection){
    var copiedObject = jQuery.extend(true,{}, selection);//deep copy
        
    if(copiedObject.style!=undefined){
        delete copiedObject.style;
    }if(copiedObject.labelres!=undefined){
        delete copiedObject.labelres;
    }if(copiedObject.surface!=undefined){
        delete copiedObject.surface;
    }
    return copiedObject;
}

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

var createOtherModelSpec = function(spec,type){
    var other=$('<li/>',{
        text:type,
        "class":type.toLowerCase(),
    });

    var attributes = $('<ul/>',{
        "class":type.toLowerCase()+'_attributes',
    }).appendTo(other);

    for(var attribute_index in spec){
        createAttribute(attribute_index,spec[attribute_index]).appendTo(attributes);
    }
    return other;
}

var createStyleSpec = function(style_spec_object,style_spec_type,model_spec_type,selection_index){
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

    for(var attribute_index in style_spec_object){
        createAttribute(attribute_index,style_spec_object[attribute_index]).appendTo(style_spec_attributes);
    }
    var add_attribute = $('<button/>',{
        "class":"add_attribute",
        "text":"Add Attribute",
        "data-index":selection_index,
        "data-type":model_spec_type,
        "data-styletype":style_spec_type,
        "click":function(){addAttribute(this)},
    }).appendTo(style_spec);
                        
    return style_spec;
}

var createStyle = function(model_spec_object,model_spec_type,selection_index){
    var style=$('<span/>',{
        text:"Style",
        "class":"style",
    });

    var style_specs = $('<ul/>',{
        "class":'style_specs',
    }).appendTo(style);
                        
    for(var attribute_index in model_spec_object){
        createStyleSpec(model_spec_object[attribute_index],attribute_index,model_spec_type,selection_index).appendTo(style_specs);
    }

    return style; 
}

var createModelSpecification = function(model_spec_type,model_spec_object,selection_index){
    var model_specification = null;
    //check for type
    if(model_spec_type=="style"){
        model_specification = createStyle(model_spec_object,model_spec_type,selection_index)

        var add_style_spec = $('<button/>',{
            "class":"add_style_spec",
            "text":"Add Style Spec",
            "data-index":selection_index,
            "data-type":model_spec_type,
            "click":function(){addStyleSpec(this)},
        }).appendTo(model_specification);
         
    }else if(model_spec_type=="surface"){
        model_specification = createOtherModelSpec(model_spec_object,"Surface")
    }else if(model_spec_type=="labelres"){
        model_specification = createOtherModelSpec(model_spec_object,"LabelRes")
    }             

    return model_specification;
}

var createModelSpecButtons = function(selection_index,selection){

    //add model spec
    var add_style = $('<button/>',{
        "class":"add_style",
        "text":"Add Style",
        "data-index":selection_index,
        "click":function(){addModelSpec("style",this)},
    }).appendTo(selection);

    var add_surface = $('<button/>',{
        "class":"add_surface",
        "text":"Add Surface",
        "data-index":selection_index,
        "click":function(){addModelSpec("surface",this)},
    }).appendTo(selection);

    var add_labelres = $('<button/>',{
        "class":"add_labelres",
        "text":"Add LabelRes",
        "data-index":selection_index,
        "click":function(){addModelSpec("labelres",this)},
    }).appendTo(selection);
}
//this function creates the selection object
var createSelection = function(selection_object, selection_index,selection_booleans){

    //creates container
    var selection = $("<li/>",{
        class:"selection"
    });
            
    //add together sub selections
    var attribute_pairs =[];
    var sel = augmentSelection(selection_object)
    for(var subselection in sel){
        var obj=sel[subselection];
        attribute_pairs.push(subselection+":"+obj);
    }

    var modifier=attribute_pairs.join(";");

    var selection_spec=$('<div/>', {
        class:'selection_spec',
        text:modifier,
        contenteditable:'true',
     }).appendTo(selection);        


    //check if style exists and if so create the object
    if(selection_object.style !=null && !selection_booleans.style){
        var style = createModelSpecification("style",selection_object.style, selection_index);
        style.appendTo(selection);
        selection_booleans.style=true;
    }else if(selection_object.surface !=null && !selection_booleans.surface){
        var surface = createModelSpecification("surface", selection_object.surface, selection_index);
        surface.appendTo(selection);
        selection_booleans.surface=true;
    }else if(selection_object.labelres != null && !selection_booleans.labelres){
        var labelres= createModelSpecification("labelres", selection_object.labelres, selection_index);
        labelres.appendTo(selection);
        selection_booleans.labelres=true;
    }

    createModelSpecButtons(selection_index,selection);

    return selection;
}
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

        var selection_count = 0;

        var selection_booleans = {
            surface:false,
            style:false,
            labelres:false,
        }

        if(selection_object.surface !=undefined)
            selection_count++;
        if(selection_object.style != undefined)
            selection_count++;
        if(selection_object.labelres !=undefined)
            selection_count++;
        if(selection_count==0)//empty selection
            selection_count++;
        //creates individual selections for each surface, style and labelres
        for(var i=0;i<selection_count;i++){
            var selection=createSelection(selection_object,selection_index,selection_booleans)
            selection.appendTo(parent);
        }
    }
}

//reads the html tree and redefines the query object based on its contents
var readHTMLTree = function(tree){

}

Object.size = function(obj) {
    var size = 0, key;
    for (key in obj) {
        if (obj.hasOwnProperty(key)) size++;
    }
    return size;
};

//takes the queyr object and creates a url for it
var unpackQuery = function(query){
    var url = "";
    //unpacks everything except for style which has multiple layers 
    var unpackObject = function (object){
        var objs =[]
        $.each(object, function(key,value){
            //array values 
            if(Array.isArray(value)){
                //sperate by commas
                objs.push(key+":"+value.join(","));
            }else{
                objs.push(key+":"+value);
            }
        });
        return objs.join(";");
    }

    var unpackStyle = function(object){
        var subStyles=[]
        $.each(object, function(sub_style,sub_style_object){
            var string="";
            string+=sub_style;
            if(Object.size(sub_style_object)!=0)
                string+=":";
            var assignments =[]
            $.each(sub_style_object, function(key,value){
                assignments.push(key+"~"+value);
            });
            string+=assignments.join(",");
            subStyles.push(string)
        });

        return subStyles.join(";");
    }
    //global 

    var unpackSelection = function(object){
        var copiedObject = jQuery.extend(true,{}, object)
        var objs=[];
        var string="";
        
        if(copiedObject.style!=undefined){
            objs.push("style="+unpackStyle(copiedObject.style));
            delete copiedObject.style;
        }if(copiedObject.labelres!=undefined){
            objs.push("labelres="+unpackObject(copiedObject.labelres));
            delete copiedObject.labelres;
        }
        if(copiedObject.surface!=undefined){
            objs.push("surface="+unpackObject(copiedObject.surface));
            delete copiedObject.surface;
        }
        var select="select="+unpackObject(copiedObject);
        objs.unshift(select);//prepend
        return objs.join("&");
    }

    var objects = [];

    objects.push(query.file.type+"="+query.file.path);

    if(query.labelres!=null){
        objects.push("labelres="+unpackObject(query.labelres));
    }else if(query.surface!=null){
        objects.push("surface="+unpackObject(query.surface));
    }else if(query.style!=null){
        objects.push("style="+unpackObject(query.style));
    }
    for(var selection in query.selections){
        objects.push(unpackSelection(query.selections[selection]))
    }

    return objects.join("&");
}
function File(path,type){
    this.path=path;
    this.type=type;
}

var Query = function(){
    this.style = null;
    this.selections = [];
    this.file = new File();
    this.surface=null;
    this.labelres=null;
}

function setURL(urlPath){
    window.history.pushState({"html":"test","pageTitle":"test"},"", "viewer.html?"+urlPath);
}


//takes the search url string and makes a query object for it 
var parseURL = function(url){
    var query = new Query();
    var tokens = url.split("&");
    
    //still using indexOf because otherwise i would need to check to see if the first substring in the string is "select" and check to see if the string isnt to small
    function stringType(string){
        if(string.indexOf("select")==0)
            return "select"
        else if(string.indexOf("pdb=")==0 || string.indexOf("cid=")==0 || string.indexOf("url=")==0)
            return "file"
        else if(string.indexOf("style")==0)
            return "style"
        else if(string.indexOf("surface")==0)
            return "surface"
        else if(string.indexOf("labelres")==0)
            return "labelres"
        return null;
    }

    var currentSelection = null;
    for(var token in tokens){
        var strings = tokens[token].split("=");
        var type = stringType(tokens[token]);
        var string = strings[1];
        var object = $3Dmol.specStringToObject(string);

        if(type == "file"){
            query.file = new File(string,strings[0]);
        }else if(type == "select"){
            var selection = object;
            query.selections.push(selection);
            currentSelection = selection;
        }else if(type == "style"){
            if(currentSelection==null)
                query.style = object;
            else
                currentSelection.style = object;
        }else if(type == "surface"){
            if(currentSelection==null)
                query.surface = object;
            else
                currentSelection.surface = object;
        }else if(type == "labelres"){
             if(currentSelection==null)
                query.labelres = object;
            else
                currentSelection.labelres = object;
        }
    }
    return query;
}

//these functions all edit the query object 
var addSelection = function(){
    query.selections.push({})
    buildHTMLTree(query);
}

var addModelSpec = function(type,selection){
    var current_selection = query.selections[selection.dataset.index]
    if(type == "style" || type == "surface" || type == "labelres"){
        if(current_selection[type]==null)
            current_selection[type]={};
        else
            console.log(type+" already defined for selection");//TODO error handling
    }
    
    buildHTMLTree(query);
}

var addStyleSpec = function(model_spec){
    var str=model_spec.getAttribute("obj")
    var i=str.split(",")[0];
    var type=str.split(",")[1];
    query.selections[i][type].attributes[query.selections[i][type].attributes.length]="";
    buildHTMLTree(query);
}

var addAttribute = function(style_spec){
    var list=style_spec.getAttribute("obj").split(",");
    query.selections[list[0]][list[1]].attributes[list[2]][""]="";
    buildHTMLTree(query);
}
//this function reads the form changes and upates the query accordingly
var updateQuery = function(){
    query.file.path=document.getElementById("model_input").value;
    query.file.type=document.getElementById("model_type").value;
}

var center = function(){
    glviewer.center({},1000,true);
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
    glviewer.center({},1000,true);
    buildHTMLTree(query);
}

//opens up the side bar
var openSide= function(){
    var width=400;
    document.getElementById("sidenav").style.width = width+"px";
    document.getElementById("menu").style.visibility="hidden";
    glviewer.translate(200,0,400,false);
    glviewer.render();
}
//closes the side bar
var closeSide= function(){
    document.getElementById("menu").style.visibility="visible";
    document.getElementById("sidenav").style.width = "0";

    glviewer.translate(-200,0,400,false);
    glviewer.render();
}