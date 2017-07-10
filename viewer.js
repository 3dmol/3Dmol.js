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

//removes eveyrhting but style,surface, and labelres
var removeAllButValid = function(object){
    var copy = jQuery.extend(true,{},object);
    for(var i in object){
        if(i!= "surface" || i!= "labelres" || i!= "style")
            delete copy[i]
    }
    return copy;
}

Object.size = function(obj) {
    var size = 0, key;
    for (key in obj) {
        if (obj.hasOwnProperty(key)) size++;
    }
    return size;
};

var createAttribute = function(name,value,parent){
    var attribute = $('<li/>',{
        class:'attribute'
    });
    
    var other=false;
    var validNames;
    if(parent.type == "line"  || parent.type == "stick" || parent.type== "cross" || parent.type == "sphere" || parent.type == "cartoon"){
        validNames=glviewer.getModel().validAtomStyleSpecs[parent.type].validItems;
        other=false;
    }else if(parent.type.toLowerCase() == "surface"){
        validNames = glviewer.getModel().validSurfaceSpecs;
        other=true;
    }else if(parent.type.toLowerCase() == "labelres"){
        validNames = glviewer.getModel().validLabelResSpecs;
        other=true;
    }

    var attribute_name = $('<select>',{
        class:'attribute_name',
    }).appendTo(attribute);

    $(attribute_name).change(function(){
        render();
    })

    $.each(validNames,function(key,value) {
        if(value.gui){
            attribute_name.append($("<option>").attr('value',key).text(key));
        }
    });

    attribute_name.val(name.toString())
    //delete button
    var delete_selection = $("<span/>",{
        html:"&#x2715;",
        class:"delete_attribute",
        "data-index":parent.index,
        "data-attr":name,
        "data-type":parent.type.toLowerCase(),
        "click":function(){
            if(other)
                deleteOtherAttribute(this);
            else
                deleteStyleAttribute(this);
        }
    }).appendTo(attribute); 

    var itemIsDescrete = function(key){
        if(key == "")
            return false;
        var type = validNames[key].type;
        return type == "boolean" || type == "color" || type == "colorscheme" || validNames[key].validItems!=undefined
    }
 
    var attribute_value;
    if(itemIsDescrete(name)){
        var validItemsValue;
        if( validNames[name].type != undefined)
            var type = validNames[name].type.toLowerCase();
        else 
            var type = undefined
        if(type=="boolean"){
            validItemsValue = ["true","false"];
        }else if(type == "colorscheme"){
            validItemsValue =  glviewer.getModel().validColorschemeSpecs;
        }else if(type == "color"){
            validItemsValue =  glviewer.getModel().validColorSpecs;
        }else if(type == undefined){
            validItemsValue = validNames[name].validItems;
        }

        var attribute_value = $('<select/>',{
            class:'attribute_value',
        }).appendTo(attribute);

        $.each(validItemsValue,function(key,value) {
            attribute_value.append($("<option>").attr('value',value).text(value));
        });

        attribute_value.val(value.toString())
            
    }else{
        attribute_value = $('<input/>',{
            class:'attribute_value',
            value:value,
        }).appendTo(attribute);
    }

    attribute_value.change(function(){
        render();
    });

    if(name.toString()!="" &&  attribute_value.prop("tagName") == "INPUT" && validNames[name.toString()].type =="number"){
        validNames[name.toString()].type =="number"
        attribute_value.attr("type","number")
        attribute_value.attr("step",validNames[name.toString()].step)
    }
    return attribute;
}

var createOtherModelSpec = function(spec,type,selection_index){
    var attributes = $('<ul/>',{
        "class":type.toLowerCase()+'_attributes',
    });
    
    for(var attribute_index in spec){
        console.log(spec)
        createAttribute(attribute_index,spec[attribute_index],{type:type,index:selection_index}).appendTo(attributes);
    }

    var add_attribute = $('<button/>',{
        "class":"add_attribute",
        "text":"Add Attribute",
        "data-index":selection_index,
        "data-type":type,
        "click":function(){addOtherAttribute(this)},
    }).appendTo(attributes);

    return attributes;
}

var createStyleSpec = function(style_spec_object,style_spec_type,model_spec_type,selection_index){
    var style_spec=$('<li/>',{
        "class":"style_spec",
    });

    var validNames=glviewer.getModel().validAtomStyleSpecs;

    var style_spec_name = $('<select>',{
        class:'style_spec_name',
    }).appendTo(style_spec);

    style_spec_name.change(function(){
        render();
    })

    $.each(validNames,function(key,value) {
        if(value.gui){
            style_spec_name.append($("<option>").attr('value',key).text(key));
        }
    });
    
    style_spec_name.val(style_spec_type)

    var delete_selection = $("<span/>",{
        html:"&#x2715;",
        class:"delete_style_spec",
        "data-index":selection_index,
        "data-type":model_spec_type,
        "data-attr":style_spec_type,
        "click":function(){deleteStyleSpec(this)},
    }).appendTo(style_spec); 

    var style_spec_attributes = $('<ul/>',{
        class:'style_spec_attributes',
    }).appendTo(style_spec);

    for(var attribute_index in style_spec_object){
        createAttribute(attribute_index,style_spec_object[attribute_index],{type:style_spec_type,index:selection_index}).appendTo(style_spec_attributes);
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
        "class":"style",
    }); 

    var style_specs = $('<ul/>',{
        "class":'style_specs',
    }).appendTo(style);
                        
    for(var attribute_index in model_spec_object){
        console.log(model_spec_object)
        createStyleSpec(model_spec_object[attribute_index],attribute_index,model_spec_type,selection_index).appendTo(style_specs);
    }

    return style; 
}

var validNames = {
    "style":"Style",
    "surface":"Surface",
    "labelres":"LabelRes",
}

var createModelSpecification = function(model_spec_type,model_spec_object,selection_index){
    var model_specification = null;
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
        model_specification = createOtherModelSpec(model_spec_object,"Surface",selection_index)
    }else if(model_spec_type=="labelres"){
        model_specification = createOtherModelSpec(model_spec_object,"LabelRes",selection_index)
    }             

    return model_specification;
}
//this function creates the selection object
var createSelection = function(spec,object,index,type){
    //creates container
    var selection = $("<li/>",{
        class:"selection"
    });

    var createHeader = function(){
        var selection_type = $('<p>',{
            class:'selection_type',
            text:validNames[type],
        }).appendTo(selection);

        //add together sub selections
        var attribute_pairs =[];
        for(var subselection in spec){
            var obj=spec[subselection];
            attribute_pairs.push(subselection+":"+obj);
        }

        var modifier=attribute_pairs.join(";");
        if(modifier == "")
            modifier = "all"

        var selection_spec=$('<input/>', {
            class:'selection_spec',
            value:modifier,
        }).appendTo(selection); 

        selection_spec.change(function(){
            render();
        })
    }

     //delete button
    var delete_selection = $("<div/>",{
        html:"&#x2715;",
        class:"delete_selection",
        "data-index":index,
        "data-type":"",
        "click":function(){deleteSelection(this);}
    }).appendTo(selection); 

    createHeader()    
    //check if style exists and if so create the object
    var ret = createModelSpecification(type,object, index);

    delete_selection.attr("data-type",type);
    ret.appendTo(selection);

    return selection;
}
/*
builds an html tree that goes inside of the selection portion of the viewer page
*/
var buildHTMLTree = function(query){
    //get parent object for the html tree
    var parent = $('#selection_list');
    //list file type and path
    $("#model_type").attr("value",query.file.type); 
    $("#model_type").change(function(){
        render();
    })
    $("#model_input").attr("value",query.file.path);
     $("#model_input").change(function(){
        render();
    })

    //loops through selections and creates a selection tree
    for(var selection_index in query.selections){
        var selection_object = query.selections[selection_index];

        if(selection_object.style != undefined){
            parent.append(createSelection(augmentSelection(selection_object),selection_object.style,selection_index,"style"))
        }
        if(selection_object.labelres != undefined){
            parent.append(createSelection(augmentSelection(selection_object),selection_object.labelres,selection_index,"labelres"))
        }
        if(selection_object.surface != undefined){
            parent.append(createSelection(augmentSelection(selection_object),selection_object.surface,selection_index,"surface"))
        }
    }

    var spacer = $('<li><br><br><br><br></li>').appendTo(parent)
}


//takes the queyr object and creates a url for it
var queryToURL = function(query){
    var url = "";
    //unpacks everything except for style which has multiple layers 
    var unpackOther = function (object){
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

    var unpackSelection = function(object){
        var copiedObject = jQuery.extend(true,{}, object)
        var objs=[];
        var string="";

        for(var obj in object){
            if(obj == "style"){
                objs.push("style="+unpackStyle(object.style))
            }else if(obj == "labelres" || obj == "surface"){
                objs.push(obj+"="+unpackOther(object[obj]))
            }
        }
        /*
        if(copiedObject.style!=undefined){
            objs.push("style="+unpackStyle(copiedObject.style));
            delete copiedObject.style;
        }if(copiedObject.labelres!=undefined){
            objs.push("labelres="+unpackOther(copiedObject.labelres));
            delete copiedObject.labelres;
        }
        if(copiedObject.surface!=undefined){
            objs.push("surface="+unpackOther(copiedObject.surface));
            delete copiedObject.surface;
        }*/
        var select="select="+unpackOther(augmentSelection(object));
        objs.unshift(select);//prepend
        return objs.join("&");
    }

    var objects = [];

    objects.push(query.file.type+"="+query.file.path);

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
    this.selections = [];
    this.file = new File();
}

function setURL(urlPath){
    window.history.pushState({"html":"test","pageTitle":"test"},"", "viewer.html?"+urlPath);
}
//this function will look through the dictionaries defined in glmodel and validate if the types are correct and return a dictionary with flags for the types that are incorecct
var validateQuery = function(query){

}
var count = 0;
//takes the search url string and makes a query object for it 
var urlToQuery = function(url){
    var query = new Query();
    var tokens = url.split("&");
    //still using indexOf because otherwise i would need to check to see if the first substring in the string is "select" and check to see if the string isnt to small
    function stringType(string){
        if(string ==  "select")
            return "select"
        else if(string =="pdb" || string == "cid" || string == "url")
            return "file"
        else if(string == "style" || string == "surface" || string == "labelres"){
            count++;
            return string;
        }        
        throw "Illegal url string : "+string;
        return;
    }

    var currentSelection = {};
    for(var token in tokens){
        var strings = tokens[token].split("=");
        var type = stringType(strings[0]);//left side of equals
        var string = strings[1];//right side of equals
        var object = $3Dmol.specStringToObject(string);

        if(type == "file"){
            query.file = new File(string,strings[0]);
        }else if(type == "select"){
            var selection = object;
            query.selections.push(selection);
            currentSelection = selection;
        }else if(type == "style" || type=="surface" || type == "labelres"){
            currentSelection[type] = object;
        }
    }
    return query;
}

var updateQueryFromHTML = function(){
    //File/PDB/URL updating
    query.file.path=document.getElementById("model_input").value;
    query.file.type=document.getElementById("model_type").value;

    $("#model_type").change(function(){
        render();
    });

    var updateOther = function(other){
        var object={};
     
        var otherList = $(other).children(".attribute");
        otherList.each(function(li){
            object[$(otherList[li]).children(".attribute_name")[0].value]=$(otherList[li]).children(".attribute_value")[0].value
        });
        return object;
    }

    var updateStyle = function(styl){
        var object={};
        var list = $(styl).children(".style_specs");
        list = $(list).children(".style_spec")
        list.each(function(li){
            var subtype=$(list[li]).children(".style_spec_name")[0].value;
            object[subtype]={};
            var otherList =$(list[li]).children(".style_spec_attributes")[0];
            otherList=$(otherList).children(".attribute")
            otherList.each(function(li){
                var tag=object[subtype][$(otherList[li]).children(".attribute_name")[0].value]=$(otherList[li]).children(".attribute_value")[0].tagName
                object[subtype][$(otherList[li]).children(".attribute_name")[0].value]=$(otherList[li]).children(".attribute_value")[0].value;
            });
            
        });
        return object;
    }

    var updateSelectionElements = function(selection_string){
        return $3Dmol.specStringToObject(selection_string);
    }

    function arraysEqual(a, b) {
        if (a === b) return true;
        if (a == null || b == null) return false;
        if (a.length != b.length) return false;

        for (var i = 0; i < a.length; ++i) {
            if (a[i] !== b[i]) return false;
        }
        return true;
    }

    var isSame = function(obj1,obj2){
        for(var key in obj1){
            if(Array.isArray(obj1[key])){
                if(Array.isArray(obj2[key]))
                    return arraysEqual(obj1[key],obj2[key])
                return false;
            }
            if(obj2[key]==undefined || obj2[key] != obj1[key])
                return false;
        }
        return true;
    }

    function combine(obj1, src1) {
        for (var key in src1) {
            if (src1.hasOwnProperty(key)) obj1[key] = src1[key];
        }   
        return obj1;
    }      

    var selects = [];
    var listItems = $(".selection")
    listItems.each(function(index,value){
        if(listItems.hasOwnProperty(index)){
            var getSubObject = function(index){
                var attr = $(value);
                var attribute=attr[0]
                var type=$(attribute).children()[1].innerHTML.toLowerCase()

                if(type=="style"){
                    var style =updateStyle($(attribute).children(".style")[0])
                    return {"style":style}
                }else if(type=="surface"){
                    var surface = updateOther($(attribute).children(".surface_attributes")[0])
                    return {"surface":surface} 
                }else if(type == "labelres"){
                    var labelres = updateOther($(attribute).children(".labelres_attributes")[0])
                    return {"labelres":labelres} 
                }
            }

            var val = getSubObject(index);
            var selection_spec = $(listItems[index]).children(".selection_spec")[0].value;
            var selection = updateSelectionElements(selection_spec);
            var extended = combine(selection,val)
            selects.push(extended)
              

        }
    });

    var final_selections = [];
    var used=[]
    for(var sele in selects){
        var augmented = augmentSelection(selects[sele]);
        if(used.includes(sele))
            continue;
        used.push(sele)
        for(var sele1 in selects){
            if(sele== sele1 || used.includes(sele1))
                continue;
            if(isSame(augmentSelection(selects[sele1]),augmented)){
               augmented = combine(selects[sele1],augmented);
               used.push(sele1);
            }
        }
        final_selections.push(combine(selects[sele],augmented))
    }

    query.selections=final_selections;
}

var query = urlToQuery(window.location.search.substring(1));
//this function compresses the html object back into a url
var render = function(){
    //calls update query
    updateQueryFromHTML();
    setURL(queryToURL(query));
    buildHTMLTree(query);
    run();
}
//these functions all edit the query object 
var addSelection = function(type){
    count++;
    if(type == "style")      
        query.selections.push({"style":{}})
    else if(type == "surface")
        query.selections.push({"surface":{}})
    else if(type == "labelres")
        query.selections.push({"labelres":{}})

    buildHTMLTree(query);
}

var deleteSelection = function(spec){
    delete query.selections[spec.dataset.index][spec.dataset.type];
    if(query.selections[spec.dataset.index].surface == undefined && query.selections[spec.dataset.index].style == undefined && query.selections[spec.dataset.index].labelres == undefined)
        delete query.selections[spec.dataset.index]
    
    buildHTMLTree(query);
    render();
}

var addModelSpec = function(type,selection){
    var current_selection;
    current_selection = query.selections[selection.dataset.index]
    
    if(type == "style" || type == "surface" || type == "labelres"){
        if(current_selection[type]==null)
            current_selection[type]={};
        else
            console.err(type+" already defined for selection");//TODO error handling
    }
    
    buildHTMLTree(query);
}

var addStyleSpec = function(model_spec){
    var defaultKey = "";
    var defaultValue = {};
    query.selections[model_spec.dataset.index][model_spec.dataset.type][defaultKey]=defaultValue;
    
    buildHTMLTree(query);
}

var deleteStyleSpec = function(spec){
    delete query.selections[spec.dataset.index][spec.dataset.type][spec.dataset.attr]
    
    buildHTMLTree(query);
    render();
}

var addOtherAttribute= function(spec){
    var defaultKey = "";
    var defaultValue = "";
    query.selections[spec.dataset.index][spec.dataset.type.toLowerCase()][defaultKey]=defaultValue;
    
    buildHTMLTree(query);
}

var deleteOtherAttribute = function(spec){
    delete query.selections[spec.dataset.index][spec.dataset.type][spec.dataset.attr]
    
    buildHTMLTree(query);
    render();
}

var addAttribute = function(style_spec){
    var defaultKey = "";
    var defaultValue = "";
    query.selections[style_spec.dataset.index][style_spec.dataset.type][style_spec.dataset.styletype][defaultKey]=defaultValue;
    buildHTMLTree(query);
}

var deleteStyleAttribute = function(spec){
    delete query.selections[spec.dataset.index]["style"][spec.dataset.type][spec.dataset.attr]
    buildHTMLTree(query);
    render();
}
//this function reads the form changes and upates the query accordingly
var center = function(){
    glviewer.center({},1000,true);
}
//initializes the sidebar based on the given url
var initSide = function(url){
    var list = document.createElement('ul')
    document.getElementById('container').appendChild(list);
    //updating on back button
    $(window).on('popstate', function() {
        query = urlToQuery(window.location.search.substring(1));
        buildHTMLTree(query);
        render();
    });
}
//opens up the side bar
var openSide= function(){
    var width=420;
    document.getElementById("sidenav").style.width = width+"px";
    document.getElementById("menu").style.visibility="hidden";
    document.getElementById("header").style.visibility="visible";
    buildHTMLTree(query);
    glviewer.translate(width/2,0,400,false);
    glviewer.render();
}
//closes the side bar
var closeSide= function(){
    document.getElementById("menu").style.visibility="visible";
    document.getElementById("sidenav").style.width = "0";
    document.getElementById("header").style.visibility="hidden";
    //todo make resize dynamic
    glviewer.translate(-200,0,400,false);
    glviewer.render();
}