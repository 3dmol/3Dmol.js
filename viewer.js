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

        var selection_booleans ={
            surface:false,
            style:false,
            labelres:false,
        }
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
                    for(var attribute_index in model_spec_object.attributes){
                        createStyleSpec(model_spec_object.attributes[attribute_index],attribute_index,model_spec_type).appendTo(style_specs);
                    }

                    return style; 
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
                        var object_keys = Object.keys(spec[attribute_index]);
                        createAttribute(object_keys[0],spec[attribute_index][object_keys[0]]).appendTo(attributes);
                    }
                    return other;
                }
                //check for type
                if(model_spec_type=="style"){
                   model_specification = createStyle(model_spec_object.attributes)
                    var add_style_spec = $('<button/>',{
                    "class":"add_style_spec",
                    "text":"Add Style Spec",
                    "data-index":selection_index,
                    "data-type":model_spec_type,
                    "click":function(){addStyleSpec(this)},
                }).appendTo(model_specification);
                }else if(model_spec_type=="surface"){
                    model_specification = createOtherModelSpec(model_spec_object.attributes,"Surface")
                }else if(model_spec_type=="labelres"){
                    model_specification = createOtherModelSpec(model_spec_object.attributes,"LabelRes")
                }             

                return model_specification;
            }
            //check if style exists and if so create the object
            if(selection_object.style !=null && !selection_booleans.style){
                var style = createModelSpecification("style",selection_object.style);
                //add style to model_specifications
                style.appendTo(model_specifications);
                selection_booleans.style=true;
            }else if(selection_object.surface !=null && !selection_booleans.surface){
                var surface = createModelSpecification("surface", selection_object.surface)
                surface.appendTo(model_specifications);
                selection_booleans.surface=true;
            }else if(selection_object.labelres != null && !selection_booleans.labelres){
                var labelres= createModelSpecification("labelres", selection_object.labelres)
                labelres.appendTo(model_specifications);
                selection_booleans.labelres=true;
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
        var selection_count = 0;

        if(selection_object.surface !=undefined)
            selection_count++;
        if(selection_object.style != undefined)
            selection_count++;
        if(selection_object.labelres !=undefined)
            selection_count++;
        for(var i=0;i<selection_count;i++){
            var selection=createSelection()
            selection.appendTo(parent);
        }
        //creates a style tree
    }
}

Object.size = function(obj) {
    var size = 0, key;
    for (key in obj) {
        if (obj.hasOwnProperty(key)) size++;
    }
    return size;
};

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
            console.log(sub_style_object)
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
        console.log(query.selections[selection])
        objects.push(unpackSelection(query.selections[selection]))
    }
    console.log(objects);

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
    console.log(unpackQuery(query))
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
console.log($3Dmol.specStringToObject(window.location.search.substring(1)));
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
    //document.getElementById("url").value=window.location.href.substring(window.location.href.indexOf("?")+1);
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