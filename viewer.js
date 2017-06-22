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

var createOtherModelSpec = function(spec,type,selection_index){
    var attributes = $('<ul/>',{
        "class":type.toLowerCase()+'_attributes',
    });

    var createAttribute = function(name,value){
        var attribute = $('<li/>',{
            class:'attribute'
        });

        var validNames;
        if(type.toLowerCase() == "surface"){
            validNames = glviewer.getModel().validSurfaceSpecs;
        }else if(type.toLowerCase() == "labelres"){
            validNames = glviewer.getModel().validLabelResSpecs;
        }

        var attribute_name = $('<select>',{
            class:'attribute_name',
        }).appendTo(attribute);

        $.each(validNames,function(key,value) {
            if(value.valid){
                attribute_name.append($("<option>").attr('value',key).text(key));
            }
        });

        attribute_name.val(name)
        //delete button
        var delete_selection = $("<span/>",{
            html:"&#x2715;",
            class:"delete_attribute",
            "data-index":selection_index,
            "data-attr":name,
            "data-type":type.toLowerCase(),
            "click":function(){deleteOtherAttribute(this)},
        }).appendTo(attribute); 

        var attribute_value = $('<span/>',{
            class:'attribute_value',
            text:value,
            contenteditable:'true',
        }).appendTo(attribute);

        return attribute;
    }

    for(var attribute_index in spec){
        createAttribute(attribute_index,spec[attribute_index]).appendTo(attributes);
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

    $.each(validNames,function(key,value) {
        if(value.valid){
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

    var createAttribute = function(name,value){
        var attribute = $('<li/>',{
            class:'attribute'
        });

        var validItems = validNames[style_spec_type].validItems;
        console.log()
        var attribute_name = $('<select>',{
            class:'attribute_name',
        }).appendTo(attribute);

        $.each(validItems,function(key,value) {
            if(value.valid){
                attribute_name.append($("<option>").attr('value',key).text(key));
            }
        });
        attribute_name.val(name)
        //delete button
        console.log(name)
        var delete_selection = $("<div/>",{
            html:"&#x2715;",
            class:"delete_attribute",
            "data-index":selection_index,
            "data-type":style_spec_type,
            "data-attr":name,
            "click":function(){deleteStyleAttribute(this)},
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
        model_specification = createOtherModelSpec(model_spec_object,"Surface",selection_index)
    }else if(model_spec_type=="labelres"){
        model_specification = createOtherModelSpec(model_spec_object,"LabelRes",selection_index)
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

    var selection_type = $('<div/>',{
        class:"selection_type",
        text:""
    }).appendTo(selection);
    //delete button
    var delete_selection = $("<div/>",{
        html:"&#x2715;",
        class:"delete_selection",
        "data-index":selection_index,
        "data-type":"",
        "click":function(){deleteSelection(this);}
    }).appendTo(selection); 
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

        selection_type.text("style");
        delete_selection.attr("data-type","style");
        style.appendTo(selection);
        selection_booleans.style=true;
    }else if(selection_object.surface !=null && !selection_booleans.surface){
        var surface = createModelSpecification("surface", selection_object.surface, selection_index);

        selection_type.text("surface");
        delete_selection.attr("data-type","surface");
        surface.appendTo(selection);
        selection_booleans.surface=true;
    }else if(selection_object.labelres != null && !selection_booleans.labelres){
        var labelres = createModelSpecification("labelres", selection_object.labelres, selection_index);

        selection_type.text("labelres");     
        delete_selection.attr("data-type","labelres");
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

Object.size = function(obj) {
    var size = 0, key;
    for (key in obj) {
        if (obj.hasOwnProperty(key)) size++;
    }
    return size;
};

//takes the queyr object and creates a url for it
var queryToURL = function(query){
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
var urlToQuery = function(url){
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

var updateQueryFromHTML = function(){
    //File/PDB/URL updating
    query.file.path=document.getElementById("model_input").value;
    query.file.type=document.getElementById("model_type").value;

    var updateOther = function(other){
        var object={};
     
        var otherList = $(other).children(".attribute");
        otherList.each(function(li){
            object[$(otherList[li]).children(".attribute_name")[0].innerHTML]=$(otherList[li]).children(".attribute_value")[0].innerHTML
        });
        return object;
    }

    var updateStyle = function(styl){
        var object={};
        var list = $(styl).children(".style_specs");
        list = $(list).children(".style_spec")
        list.each(function(li){
            var subtype=$(list[li]).children(".style_spec_name")[0].innerHTML;
            object[subtype]={};
            var otherList =$(list[li]).children(".style_spec_attributes")[0];
            otherList=$(otherList).children(".attribute")
            otherList.each(function(li){
                object[subtype][$(otherList[li]).children(".attribute_name")[0].innerHTML]=$(otherList[li]).children(".attribute_value")[0].innerHTML;
            });
            
        });
        return object;
    }

    var updateSelectionElements = function(selection_string){
        return $3Dmol.specStringToObject(selection_string)
    }

    function arraysEqual(a, b) {
        if (a === b) return true;
        if (a == null || b == null) return false;
        if (a.length != b.length) return false;

        // If you don't care about the order of the elements inside
        // the array, you should sort both arrays here.

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

    var selects = [];
    var listItems = $(".selection")
    listItems.each(function(index,value){
        if(listItems.hasOwnProperty(index)){
            var getSubObject = function(index){
                var attr = $(value);
                var attribute=attr[0]
                var type=$(attribute).children()[0].innerHTML

                if(type=="style"){
                    return {style:updateStyle($(attribute).children(".style")[0])}
                }else if(type=="surface"){
                    return {"surface":updateOther($(attribute).children(".surface_attributes")[0])} 
                }else if(type == "labelres"){
                    return {"labelres":updateOther($(attribute).children(".labelres_attributes")[0])} 
                }
            }

            function extend(obj1, src1) {
                for (var key in src1) {
                    if (src1.hasOwnProperty(key)) obj1[key] = src1[key];
                }   
                return obj1;
            }            
            var val=getSubObject(index);
            var selection = updateSelectionElements($(listItems[index]).children(".selection_spec")[0].innerHTML);
            var extended=extend(selection,val)
            console.log(extended)
            selects.push(extended)
        }
    });

    function extend(obj1, src1) {
        for (var key in src1) {
            if (src1.hasOwnProperty(key)) obj1[key] = src1[key];
        }   
        return obj1;
    }

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
               augmented = extend(selects[sele1],augmented);
               used.push(sele1);
            }
        }
        final_selections.push(extend(selects[sele],augmented))
    }

    query.selections=final_selections;
}

var query = urlToQuery(window.location.search.substring(1));
//this function compresses the html object back into a url
var render = function(){
    //calls update query
    updateQueryFromHTML();
    setURL(queryToURL(query));
    buildHTMLTree(query)

}
//these functions all edit the query object 
var addSelection = function(){
    query.selections.push({})
    buildHTMLTree(query);
}

var deleteSelection = function(spec){
    delete query.selections[spec.dataset.index][spec.dataset.type];
    if(query.selections[spec.dataset.index].surface == undefined && query.selections[spec.dataset.index].style == undefined && query.selections[spec.dataset.index].labelres == undefined)
        delete query.selections[spec.dataset.index]
    
    buildHTMLTree(query);
    render();
    run();
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
    query.selections[model_spec.dataset.index][model_spec.dataset.type][""]={};
    buildHTMLTree(query);
}

var deleteStyleSpec = function(spec){
    delete query.selections[spec.dataset.index][spec.dataset.type][spec.dataset.attr]
    
    buildHTMLTree(query);
    render();
    run();
}

var addOtherAttribute= function(spec){
    query.selections[spec.dataset.index][spec.dataset.type.toLowerCase()][""]="";
    buildHTMLTree(query)
}

var deleteOtherAttribute = function(spec){
    delete query.selections[spec.dataset.index][spec.dataset.type][spec.dataset.attr]
    
    buildHTMLTree(query);
    render();
    run();
}

var addAttribute = function(style_spec){
    query.selections[style_spec.dataset.index][style_spec.dataset.type][style_spec.dataset.styletype][""]="";
    buildHTMLTree(query);
}

var deleteStyleAttribute = function(spec){
    delete query.selections[spec.dataset.index]["style"][spec.dataset.type][spec.dataset.attr]
    
    buildHTMLTree(query);
    render();
    run();
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
        run();
    });

}
//opens up the side bar
var openSide= function(){
    var width=420;
    document.getElementById("sidenav").style.width = width+"px";
    document.getElementById("menu").style.visibility="hidden";

    buildHTMLTree(query);
    glviewer.translate(width/2,0,400,false);
    glviewer.render();
}
//closes the side bar
var closeSide= function(){
    document.getElementById("menu").style.visibility="visible";
    document.getElementById("sidenav").style.width = "0";
    //todo make resize dynamic
    glviewer.translate(-200,0,400,false);
    glviewer.render();
}