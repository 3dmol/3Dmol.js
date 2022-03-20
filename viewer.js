
// removes style,labelres, and surface from a copy of the selection object and returns it
const augmentSelection = function(selection){
    const copiedObject = jQuery.extend(true,{}, selection);// deep copy
        
    if(copiedObject.style!==undefined){
        delete copiedObject.style;
    }if(copiedObject.labelres!==undefined){
        delete copiedObject.labelres;
    }if(copiedObject.surface!==undefined){
        delete copiedObject.surface;
    }

    return copiedObject;
}

// removes eveyrhting but style,surface, and labelres
const removeAllButValid = function(object){
    const copy = jQuery.extend(true,{},object);
    for(const i in object){
        if(i!== "surface" || i!== "labelres" || i!== "style")
            delete copy[i]
    }
    return copy;
}

Object.size = function(obj) {
    let size = 0; let key;
    for (key in obj) {
        if (obj.hasOwnProperty(key)) size+=1;
    }
    return size;
};

const createAttribute = function(name,value,parent){
    const attribute = $('<li/>',{
        class:'attribute'
    });
    
    let other=false;
    let validNames;
    let type;
    if(parent.type === "line"  || parent.type === "stick" || parent.type=== "cross" || parent.type === "sphere" || parent.type === "cartoon"){
        type = "style";
        validNames=$3Dmol.GLModel.validAtomStyleSpecs[parent.type].validItems;
        other=false;
    }else if(parent.type.toLowerCase() === "surface"){
        type = "surface";
        validNames = $3Dmol.GLModel.validSurfaceSpecs;
        other=true;
    }else if(parent.type.toLowerCase() === "labelres"){
        type = "labelres";
        validNames =$3Dmol.GLModel.validLabelResSpecs;
        other=true;
    }else if(name !== ""){
        // undefined name
        return undefined;
    }

    if(validNames[name] === undefined && name!=="")
        return undefined

    const attribute_name = $('<select>',{
        class:'attribute_name',
    }).appendTo(attribute);
    const obj_type = type;

    $.each(validNames,(key,value) => {
        if(value.gui){
            attribute_name.append($("<option>").attr('value',key).text(key));
        }
    });

    attribute_name.val(name.toString())
    if(name.toString() === ""){
        let list;
        if(type === "style")
            list = query.selections[parent.index][type][parent.type];
        else
            list = query.selections[parent.index][type];

        let index;
        for(const i in validNames){
            if(validNames[i].gui && list[i] === undefined){
                index = i;
                break;
            }
        }
        name = index;

        if(name === undefined)
            return;// all of the attribute names are being used
        attribute_name.val(index)
    }

    // delete button
    const delete_selection = $("<span/>",{
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

    const itemIsDescrete = function(key){
        if(key === "")
            return false;
        const {type} = validNames[key];
        return type === "boolean" || type === "color" || type === "colorscheme" || validNames[key].validItems!==undefined
    }
    let attribute_value;
    if(itemIsDescrete(name) ){
        let validItemsValue;
        if(validNames[name].type !== undefined)
            type = validNames[name].type.toLowerCase();
        else
            type = undefined
        if(type==="boolean"){
            validItemsValue = ["false","true"];
        }else if(type === "colorscheme"){
            validItemsValue =  Object.keys($3Dmol.builtinColorSchemes).concat(['greenCarbon','cyanCarbon','yellowCarbon','whiteCarbon','magentaCarbon']);
        }else if(type === "color"){
            validItemsValue =  Object.keys($3Dmol.htmlColors);
            if(parent.type === 'cartoon') validItemsValue.unshift('spectrum');            
        }else if(type === undefined){
            validItemsValue = validNames[name].validItems;
        }

        const attribute_value = $('<select/>',{
            class:'attribute_value',
        }).appendTo(attribute);

        $.each(validItemsValue,(key,value) => {
            attribute_value.append($("<option>").attr('value',value).text(value));
        });

        attribute_value.val(value.toString());
        if(value === ""){
            attribute_value.val(validItemsValue[0])
        } 
    }else{
        if(value === "")
            value = validNames[name].default
        attribute_value = $('<input/>',{
            class:'attribute_value',
            value,
        }).appendTo(attribute);
    }
    attribute_name.change(()=> {
        let validItemsValue;
        const {type} = validNames[attribute_name.val()]
        if(type==="boolean"){
            validItemsValue = ["false","true"];
        }else if(type === "colorscheme"){
            validItemsValue =  $3Dmol.GLModel.validColorschemeSpecs;
        }else if(type === "color"){
            validItemsValue =  $3Dmol.GLModel.validColorSpecs;
        }else if(type === undefined){
            validItemsValue = validNames[name].validItems;
        }
        const defa = validNames[attribute_name.val()].default;
        let val;
        if(validItemsValue !== undefined){
            val = validItemsValue[0];
        }else{
            val = defa
        }
        if(attribute_value.children()[0]!== undefined)
            attribute_value.children()[0].value = val;
        else
            attribute_value.val(val);
        render(obj_type === "surface");
    });
    attribute_value.change(()=> {
        render(obj_type === "surface");
    });

    if(name!=="" &&attribute_value.prop("tagName") === "INPUT" && validNames[name].type ==="number"){
        validNames[name].type =="number"
        attribute_value.attr("type","number")
        attribute_value.attr("step",validNames[name].step)
        attribute_value.addClass("spinner")
        const {max} = validNames[name];
        const {min} = validNames[name];
        if(max !== undefined)
            attribute_value.attr("max",max);
        if(min !== undefined)
            attribute_value.attr("min",min);   
    }
    return attribute;
}

const createOtherModelSpec = function(spec,type,selection_index){
    const attributes = $('<ul/>',{
        "class":`${type.toLowerCase()}_attributes`,
    });
    for(const attribute_index in spec){
        const attribute=createAttribute(attribute_index,spec[attribute_index],{type,index:selection_index})
        if(attribute !== undefined)
            attribute.appendTo(attributes);
    }

    const add_attribute = $('<button/>',{
        "class":"add_attribute",
        "text":"Add Attribute",
        "data-index":selection_index,
        "data-type":type,
        "click":function(){addOtherAttribute(this)},
    }).appendTo(attributes);

    return attributes;
}

const createStyleSpec = function(style_spec_object,style_spec_type,model_spec_type,selection_index){
    const style_spec=$('<li/>',{
        "class":"style_spec",
    });

    const validNames=$3Dmol.GLModel.validAtomStyleSpecs;

    const style_spec_name = $('<select>',{
        class:'style_spec_name',
    }).appendTo(style_spec);

    style_spec_name.change(()=> {
        const obj = query.selections[selection_index].style[style_spec_type];
        for(const i in obj){
            if(!validNames[style_spec_name.val()].validItems.hasOwnProperty(i)){
                delete query.selections[selection_index].style[style_spec_type][i];
            }
        }
        query.selections[selection_index].style[style_spec_name.val()]=query.selections[selection_index].style[style_spec_type];
        delete query.selections[selection_index].style[style_spec_type];
        buildHTMLTree(query)
        render();
    });

    $.each(validNames,(key,value) => {
        if(value.gui){
            style_spec_name.append($("<option>").attr('value',key).text(key));
        }
    });
    
    style_spec_name.val(style_spec_type.toString())
    if(style_spec_type === ""){
        const list = query.selections[selection_index].style;
        let index=0
        for(const i in validNames){
            if(validNames[i].gui && list[i] === undefined){
                index = i;
                break;
            }
        }
        if(index === 0)
            return;
        style_spec_name.val(index)
    }

    const delete_selection = $("<span/>",{
        html:"&#x2715;",
        class:"delete_style_spec",
        "data-index":selection_index,
        "data-type":model_spec_type,
        "data-attr":style_spec_type,
        "click":function(){deleteStyleSpec(this)},
    }).appendTo(style_spec); 

    const style_spec_attributes = $('<ul/>',{
        class:'style_spec_attributes',
    }).appendTo(style_spec);

    for(const attribute_index in style_spec_object){
        const attribute = createAttribute(attribute_index,style_spec_object[attribute_index],{type:style_spec_type,index:selection_index})
        if(attribute !== undefined)
            attribute.appendTo(style_spec_attributes);
    }

    const add_attribute = $('<button/>',{
        "class":"add_attribute",
        "text":"Add Attribute",
        "data-index":selection_index,
        "data-type":model_spec_type,
        "data-styletype":style_spec_type,
        "click":function(){addAttribute(this)},
    }).appendTo(style_spec);
                        
    return style_spec;
}

const createStyle = function(model_spec_object,model_spec_type,selection_index){
    const style=$('<span/>',{
        "class":"style",
    }); 

    const style_specs = $('<ul/>',{
        "class":'style_specs',
    }).appendTo(style);
                        
    for(const attribute_index in model_spec_object){
        const spec = createStyleSpec(model_spec_object[attribute_index],attribute_index,model_spec_type,selection_index)
        if(spec!== undefined)
            spec.appendTo(style_specs);
    }

    const add_style_spec = $('<button/>',{
        "class":"add_style_spec",
        "text":"Add Style Spec",
        "data-index":selection_index,
        "data-type":model_spec_type,
        "click":function(){addStyleSpec(this)},
    }).appendTo(style);

    return style; 
}

const validNames = {
    "style":"Style",
    "surface":"Surface",
    "labelres":"LabelRes",
}

const createModelSpecification = function(model_spec_type,model_spec_object,selection_index){
    let model_specification = null;
    if(model_spec_type==="style"){
        model_specification = createStyle(model_spec_object,model_spec_type,selection_index)
    }else if(model_spec_type==="surface"){
        model_specification = createOtherModelSpec(model_spec_object,"Surface",selection_index)
    }else if(model_spec_type==="labelres"){
        model_specification = createOtherModelSpec(model_spec_object,"LabelRes",selection_index)
    }             

    return model_specification;
}
// this function creates the selection object
const createSelection = function(spec,object,index,type){
    // creates container
    const selection = $("<li/>",{
        class:"selection"
    });

    const createHeader = function(){
        const selection_type = $('<p>',{
            class:'selection_type',
            text:validNames[type],
        }).appendTo(selection);

        // add together sub selections
        const attribute_pairs =[];
        for(const subselection in spec){
            let obj=spec[subselection];
            if(typeof(obj) === 'object' && Object.keys(obj).length === 0)
                obj = ""; // empty object
            attribute_pairs.push(`${subselection}:${obj}`);
        }

        let modifier=attribute_pairs.join(";");
        if(modifier === "")
            modifier = "all"
        const selection_spec=$('<input/>', {
            class:'selection_spec',
            value:modifier,
        }).appendTo(selection); 

        selection_spec.change(()=> {
            render(type === "surface");
        })
    }

     // delete button
    const delete_selection = $("<div/>",{
        html:"&#x2715;",
        class:"delete_selection",
        "data-index":index,
        "data-type":"",
        "click":function(){deleteSelection(this);}
    }).appendTo(selection); 

    createHeader()    
    // check if style exists and if so create the object
    const ret = createModelSpecification(type,object, index);

    delete_selection.attr("data-type",type);
    ret.appendTo(selection);

    return selection;
}
/*
 * builds an html tree that goes inside of the selection portion of the viewer
 * page
 */
let buildHTMLTree = function(query){
    // get parent object for the html tree
    const parent = $('#selection_list');
    parent.text("");
    // list file type and path
    // $("#model_type").attr("value",query.file.type);
    document.getElementById("model_type").value = query.file.type
    $("#model_type").change(()=> {
        const val =  $("#model_type").val().toUpperCase();
        if(prev_type !== val){
            
            render(true);
            run();
        }
        prev_type = val
    })

    $("#model_input").attr("value",query.file.path);
    $("#model_input").change(()=> {
        const val =  $("#model_input").val().toUpperCase();
        if(prev_in !== val){
            if(val.match(/^[1-9][A-Za-z0-9]{3}$/) || $("#model_type").val().toLowerCase()!= "pdb"){
                render(true);
                run();
                const width = $("#sidenav").width();
            }else{
                if(prev_in!== val)
                    alert("Invalid PDB")
            }
        }
        prev_in = val;
    })
    const arr=[]
    // loops through selections and creates a selection tree
    for(const selection_index in query.selections){
        const selection_object = query.selections[selection_index];
        const aug = augmentSelection(selection_object);

        if(selection_object.style !== undefined){
            arr.push(createSelection(aug,selection_object.style,selection_index,"style"));
        }
        if(selection_object.surface !== undefined){
            arr.push(createSelection(aug,selection_object.surface,selection_index,"surface"))
        }
        if(selection_object.labelres !== undefined){
            arr.push(createSelection(aug,selection_object.labelres,selection_index,"labelres"))
        }
    }
    for(const i in arr){
        if(arr[i]!== undefined)
            parent.append(arr[i])
    }
}
// takes the queyr object and creates a url for it
const queryToURL = function(query){

    const isSame = function(obj1,obj2){
        for(const key in obj1){
            if(Array.isArray(obj1[key])){
                if(Array.isArray(obj2[key]))
                    return arraysEqual(obj1[key],obj2[key])
                return false;
            }
            if(obj2[key]===undefined || obj2[key] !== obj1[key])
                return false;
        }
        return typeof(obj1) == typeof(obj2); // {} != 0
    }
    const url = "";
    // unpacks everything except for style which has multiple layers
    const unpackOther = function (object){
        

        const objs =[]
        $.each(object, (key,value)=> {
            if(isSame(value,{}))
                value = ""
            // array values
            if(Array.isArray(value)){
                // sperate by commas
                objs.push(`${key}:${value.join(",")}`);
            }else{
                objs.push(`${key}:${value}`);
            }
        });
        return objs.join(";");
    }

    const unpackStyle = function(object){
        const subStyles=[]
        $.each(object, (sub_style,sub_style_object)=> {
            let string="";
            string+=sub_style;
            if(Object.size(sub_style_object)!==0)
                string+=":";
            const assignments =[]
            $.each(sub_style_object, (key,value)=> {
                assignments.push(`${key}~${value}`);
            });
            string+=assignments.join(",");
            subStyles.push(string)
        });

        return subStyles.join(";");
    }

    const unpackSelection = function(object){
        const copiedObject = jQuery.extend(true,{}, object)
        const objs=[];
        const string="";

        for(const obj in object){
            if(obj === "style"){
                objs.push(`style=${unpackStyle(object.style)}`)
            }else if(obj === "labelres" || obj === "surface"){
                objs.push(`${obj}=${unpackOther(object[obj])}`)
            }
        }
        const unpacked =unpackOther(augmentSelection(object));
        let select=`select=${ unpacked}`
        if(select === "select=")
            select = "select=all"
        objs.unshift(select);// prepend
        return objs.join("&");
    }

    const objects = [];
    let str = `${query.file.type}=${query.file.path}`;
    if(query.file.helper !== "")
        str+=`&type=${query.file.helper}`
    objects.push(str);

    for(const selection in query.selections){
        objects.push(unpackSelection(query.selections[selection]))
    }

    return objects.join("&");
}

function File(path,type){
    this.path=path;
    this.type=type;
    this.helper="";
}

const Query = function(){
    this.selections = [];
    this.file = new File();
}

function setURL(urlPath){
    window.history.pushState('page2',"Title", `viewer.html?${urlPath}`);
}
// this function will look through the dictionaries defined in glmodel and
// validate if the types are correct and return a dictionary with flags for the
// types that are incorecct

let count = 0;
// takes the search url string and makes a query object for it
const urlToQuery = function(url){
    // url= decodeURIComponent(url)
    if(url === "" || url.startsWith('session=') || url.startsWith('SESSION='))
        return new Query();
    
    const query = new Query();
    const tokens = url.split("&");
    // still using indexOf because otherwise i would need to check to see if the
    // first substring in the string is "select" and check to see if the string
    // isnt to small
    function stringType(string){
        if(string ===  "select")
            return "select"
        if(string ==="pdb" || string === "cid" || string === "url")
            return "file"
        if(string === "style" || string === "surface" || string === "labelres"){
            count+=1;
            return string;
        }if(string === "type"){
            return string
        }    
        throw new Error(`Illegal url string : ${string}`);
        return;
    }

    let currentSelection = null;
    for(const token in tokens){
	const uri = decodeURIComponent(tokens[token]);
	const i = uri.indexOf('=');
	const left = uri.slice(0,i);
        const type = stringType(left);// left side of first equals
        const string = uri.slice(i+1);// right side of equals
        const object = $3Dmol.specStringToObject(string);
        if(type === "file"){
            query.file = new File(string,left);
        }else if(type === "select"){
            currentSelection = object
            query.selections.push(currentSelection);
        }else if(type == "style" || type=="surface" || type == "labelres"){
            if(currentSelection == null){
                currentSelection = {}
                query.selections.push(currentSelection)
            }
            currentSelection[type] = object;
        }else if(type == type){
            query.file.helper = string;
        }
    }
    if(query.selections[0] === {})
        delete query.selections[0]
    return query;
}

const updateQueryFromHTML = function(){
    // File/PDB/URL updating
    query.file.path= $("#model_input").val(); 
    query.file.type=$("#model_type").val();

   
    const updateOther = function(other){
        const object={};
     
        const otherList = $(other).children(".attribute");
        otherList.each((li)=> {
            object[$(otherList[li]).children(".attribute_name")[0].value]=$(otherList[li]).children(".attribute_value")[0].value
        });
        return object;
    }

    const updateStyle = function(styl){
        const object={};
        let list = $(styl).children(".style_specs");
        list = $(list).children(".style_spec")
        list.each((li)=> {
            const subtype=$(list[li]).children(".style_spec_name")[0].value;
            object[subtype]={};
            let otherList =$(list[li]).children(".style_spec_attributes")[0];
            otherList=$(otherList).children(".attribute")
            otherList.each((li)=> {
                const tag=object[subtype][$(otherList[li]).children(".attribute_name")[0].value]=$(otherList[li]).children(".attribute_value")[0].tagName
                object[subtype][$(otherList[li]).children(".attribute_name")[0].value]=$(otherList[li]).children(".attribute_value")[0].value;
            });
        });
        return object;
    }

    const updateSelectionElements = function(selection_string){
        return $3Dmol.specStringToObject(selection_string);
    }

    function arraysEqual(a, b) {
        if (a === b) return true;
        if (a == null || b == null) return false;
        if (a.length != b.length) return false;

        for (let i = 0; i < a.length; ++i) {
            if (a[i] !== b[i]) return false;
        }
        return true;
    }

    const isSame = function(obj1,obj2){
        for(const key in obj1){
            if(Array.isArray(obj1[key])){
                if(Array.isArray(obj2[key]))
                    return arraysEqual(obj1[key],obj2[key])
                return false;
            }
            if(obj2[key]==undefined || obj2[key] != obj1[key])
                return false;
        }
        return typeof(obj1) == typeof(obj2); // 0 != {}
    }

    function combine(obj1, src1) {
        for (const key in src1) {
            if (src1.hasOwnProperty(key)) obj1[key] = src1[key];
        }   
        return obj1;
    }      

    const selects = [];
    const listItems = $(".selection")
    listItems.each((index,value)=> {
        if(listItems.hasOwnProperty(index) && listItems[index].id!="spacer"){
            const getSubObject = function(){
                const attr = $(value);
                const attribute=attr[0]
                const type=$(attribute).children()[1].innerHTML.toLowerCase()

                if(type=="style"){
                    const style =updateStyle($(attribute).children(".style")[0])
                    return {"style":style}
                }if(type=="surface"){
                    const surface = updateOther($(attribute).children(".surface_attributes")[0])
                    return {"surface":surface} 
                }if(type == "labelres"){
                    const labelres = updateOther($(attribute).children(".labelres_attributes")[0])
                    return {"labelres":labelres} 
                }
            }

            const val = getSubObject();
            const selection_spec = $(listItems[index]).children(".selection_spec")[0].value;
            const selection = updateSelectionElements(selection_spec);
            const extended = combine(selection,val)
            selects.push(extended)
        }
    });

    query.selections=selects;
}

var query = urlToQuery(window.location.search.substring(1));
// this function compresses the html object back into a url
var render = function(surfaceEdited){
    surfaceEdited = surfaceEdited == undefined ? false : surfaceEdited;
    // calls update query
    updateQueryFromHTML();
    const url = queryToURL(query);
    setURL(url);
    buildHTMLTree(query);
    glviewer.setStyle({},{line:{}});
    runcmds(url.split("&"),glviewer,surfaceEdited);
    glviewer.render();
}
// these functions all edit the query object
const addSelection = function(type){
    const surface  = type == "surface"
    if(type == "style")      
        query.selections.push({"style":{line:{}}})
    else if(type == "surface")
        query.selections.push({"surface":{}})
    else if(type == "labelres")
        query.selections.push({"labelres":{}})
    buildHTMLTree(query);
    render(surface);
}

var deleteSelection = function(spec){
    delete query.selections[spec.dataset.index][spec.dataset.type];
    if(query.selections[spec.dataset.index].surface == undefined && query.selections[spec.dataset.index].style == undefined && query.selections[spec.dataset.index].labelres == undefined)
        delete query.selections[spec.dataset.index]
    
    buildHTMLTree(query);
    render(spec.dataset.type == "surface");
}

const addModelSpec = function(type,selection){
    let current_selection;
    current_selection = query.selections[selection.dataset.index]
    
    if(type == "style" || type == "surface" || type == "labelres"){
        if(current_selection[type]==null)
            current_selection[type]={};
        else
            console.err(`${type} already defined for selection`);// TODO error
                                                                // handling
    }
    
    buildHTMLTree(query);
    render();
}

var addStyleSpec = function(model_spec){
    const defaultKey = "";
    const defaultValue = {};
    query.selections[model_spec.dataset.index][model_spec.dataset.type][defaultKey]=defaultValue;
    
    buildHTMLTree(query);
    render();
}

var deleteStyleSpec = function(spec){
    delete query.selections[spec.dataset.index][spec.dataset.type][spec.dataset.attr]
    
    buildHTMLTree(query);
    render();
}

var addOtherAttribute= function(spec){
    const defaultKey = "";
    const defaultValue = "";
    query.selections[spec.dataset.index][spec.dataset.type.toLowerCase()][defaultKey]=defaultValue;
    
    buildHTMLTree(query);
    render();
}

var deleteOtherAttribute = function(spec){
    delete query.selections[spec.dataset.index][spec.dataset.type][spec.dataset.attr]
    
    buildHTMLTree(query);
    render(spec.dataset.type == "surface");
}

var addAttribute = function(style_spec){
    const defaultKey = "";
    const defaultValue = "";
    query.selections[style_spec.dataset.index][style_spec.dataset.type][style_spec.dataset.styletype][defaultKey]=defaultValue;

    buildHTMLTree(query);
    render();
}

var deleteStyleAttribute = function(spec){
    delete query.selections[spec.dataset.index].style[spec.dataset.type][spec.dataset.attr]
    buildHTMLTree(query);
    render();
}
// this function reads the form changes and upates the query accordingly
const center = function(){
    glviewer.center({},1000,true);
}

const vrml = function() {
    const filename = "3dmol.wrl";
    const text = glviewer.exportVRML();
    const blob = new Blob([text], {type: "text/plain;charset=utf-8"});
    saveAs(blob, filename);
}
const savePng = function() {
    const filename = "3dmol.png";
    const text = glviewer.pngURI();
    const ImgData = text;
    const link = document.createElement('a');
    link.href = ImgData;
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
}
// initializes the sidebar based on the given url
const initSide = function(url){
    const list = document.createElement('ul')
    document.getElementById('container').appendChild(list);
    // updating on back button
    $(window).on('popstate', () => {
        query = urlToQuery(window.location.search.substring(1));
        buildHTMLTree(query);
        render(true);
    });

    buildHTMLTree(query);
}
const showCreateSession = function(){
    $('#session_list2').hide();
    $('#session_list1').toggle();
}
const showJoinSession = function(){
    $('#session_list1').hide();
  $('#session_list2').toggle();
}

let toggle = true;
let width=420;
var prev_in = $("#model_input").val();
var prev_type = $("#model_type").val();
const toggleHide =  function(){
    if(toggle){        
        $("#menu").css("display","none");
        $("#sidenav").css("width",`${width}px`);
        $('#createSession,#joinSession,#addStyle,#addSurface,#addLabelRes,#centerModel,#savePng,#vrmlExport').css("display","inline")
        glviewer.translate(width/2,0,400,false);
    }else{
        $("#sidenav").css("width","0");
        $('#createSession,#joinSession,#addStyle,#addSurface,#addLabelRes,#centerModel,#savePng,#header,#vrmlExport').css("display","none")
        $("#menu").css("display","inline");
        width = $("#sidenav").width();
        glviewer.translate(-width/2,0,400,false);
    }
    toggle = !toggle;
}

var glviewer = null;
// http://localhost/$3Dmol/viewer.html?pdb=1ycr&style=cartoon&addstyle=line&select=chain~A&colorbyelement=whiteCarbon&style=surface,opacity~.8&select=chain~B&addstyle=stick&select=chain~B,resn~TRP&style=sphere
// Process commands, in order, and run on viewer (array of strings split on '&')
var runcmds = function(cmds, viewer,renderSurface) {
            console.log("rendering")
    renderSurface = renderSurface == undefined ? true : renderSurface;
    if(renderSurface)
        viewer.removeAllSurfaces();
    viewer.removeAllLabels();
    let currentsel = {};

    for (let i = 0; i < cmds.length; i++) {
        const kv = cmds[i].split('=');
        const cmdname = kv[0];
        const cmdobj = $3Dmol.specStringToObject(kv[1]);

        if (cmdname == 'select')
            currentsel = cmdobj;
        else if (cmdname == 'surface' && renderSurface){
            viewer.addSurface($3Dmol.SurfaceType.VDW, cmdobj, currentsel,
                    currentsel);
        } else if (cmdname == 'style'){
            viewer.setStyle(currentsel, cmdobj);
        } else if (cmdname == 'addstyle'){
            viewer.addStyle(currentsel, cmdobj);
        } else if (cmdname == 'labelres'){
            viewer.addResLabels(currentsel, cmdobj);
        } else if (cmdname == 'colorbyelement'){
            if (typeof ($3Dmol.elementColors[cmdobj.colorscheme]) != "undefined")
                viewer.setColorByElement(currentsel,
                        $3Dmol.elementColors[cmdobj.colorscheme]);
        } else if (cmdname == 'colorbyproperty'){
            if (typeof (cmdobj.prop) != "undefined"
                    && typeof ($3Dmol.Gradient[cmdobj.scheme]) != "undefined"){
                viewer.setColorByProperty(currentsel, cmdobj.prop,
                        new $3Dmol.Gradient[cmdobj.scheme]());
            }
        }

    }

};
function run() {
        try {
            let url = window.location.search.substr(1);
            url= decodeURIComponent(url)
            const cmds = url.split("&");
            const first = cmds.splice(0, 1)[0];
            const pos = first.indexOf('=');
            const src = first.substring(0, pos); let data = first
                    .substring(pos + 1);
            let type = "pdb";

            if(glviewer === null) {
                glviewer = $3Dmol.createViewer("gldiv", {
                    defaultcolors : $3Dmol.rasmolElementColors
                });
                glviewer.setBackgroundColor(0xffffff);
            } else {
                glviewer.clear();
            }

            if (src == 'session' || src == 'SESSION') {
                // join a session
                joinSession(data);
            }
            if (src == 'pdb') {
                console.log(data)
                data = data.toUpperCase();
                if (!data.match(/^[1-9][A-Za-z0-9]{3}$/)) {
                    return;
                }
                data = `http://files.rcsb.org/view/${  data
                         }.pdb`;
                type = "pdb";
            } if (src == 'cif') {
                data = data.toUpperCase();
                if (!data.match(/^[1-9][A-Za-z0-9]{3}$/)) {
                    return;
                }
                data = `http://files.rcsb.org/view/${  data
                         }.cif`;
                type = "cif";
            } else if (src == 'cid') {
                type = "sdf";
                data = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${
                         data  }/SDF?record_type=3d`;
            } else if (src == 'mmtf') {
                data = data.toUpperCase();
                data = `http://mmtf.rcsb.org/full/${  data  }.mmtf`;
                type = 'mmtf';
            } else { // url
                // try to extract extension
                type = data.substr(data.lastIndexOf('.') + 1);
                if(type == 'gz') {
                    const base = data.substr(0,data.lastIndexOf('.'));
                    type = `${base.substr(base.lastIndexOf('.'))  }.gz`;
                }
            }
            if (cmds[0] && cmds[0].indexOf('type=') == 0) {
                type = cmds[0].split('=')[1];
            }

            const start = new Date();

            if (/\.gz$/.test(data) || type == 'mmtf') { // binary
                                                        // data
                $.ajax({url:data, 
                    type: "GET",
                    dataType: "binary",
                    responseType: "arraybuffer",
                    processData: false,
                    success(ret, txt, response) {
                        console.log(`mtf fetch ${  +new Date() - start  }ms`);
                        const time = new Date();
                        glviewer.addModel(ret, type);
                        runcmds(cmds, glviewer);
                        glviewer.zoomTo();
                        glviewer.render();
                        console.log(`mtf load ${  +new Date() - time  }ms`);

                }}).fail(() => {
                    // if couldn't get url natively, go through echo
                    // server
                    $.ajax({ url:"echo.cgi", 
                        data: { 'url' : data },
                        processData: false,
                        responseType: "arraybuffer",
                        dataType: "binary",
                        success(ret, txt, response) {

                        glviewer.addModel(ret, type);
                        runcmds(cmds, glviewer);
                        glviewer.zoomTo();
                        glviewer.render();
                    }})
                });
            } else {
                $.get(data, (ret, txt, response) => {
                    console.log(`alt fetch ${  +new Date() - start  }ms`);
                    const time = new Date();
                    glviewer.addModel(ret, type);
                    runcmds(cmds, glviewer);
                    glviewer.zoomTo();
                    glviewer.render();
                    console.log(`alt load ${  +new Date() - time  }ms`);

                }).fail(() => {
                    // if couldn't get url natively, go through echo
                    // server
                    $.post("echo.cgi", {
                        'url' : data
                    }, (ret, txt, response) => {
                        if(src == 'pdb' && (ret.search("We're sorry, but the requested") >= 0 || ret == "")) {
                            // really large files aren't available
                            // in pdb format
                            type = 'cif';
                            data = data.replace(/pdb$/,'cif');
                            $.post("echo.cgi",{
                                'url' : data
                            }, (ret, txt, response) => {

                                glviewer.addModel(ret, type);
                                runcmds(cmds, glviewer);
                                glviewer.zoomTo();
                                glviewer.render();
                            })
                        } else {
                            glviewer.addModel(ret, type);
                            runcmds(cmds, glviewer);
                            glviewer.zoomTo();
                            glviewer.render();
                        }
                    });
                });
            }
        }

        catch (e) {
            console
                    .error(`Could not instantiate viewer from supplied url: '${
                             e  }'`);
            window.location = "http://get.webgl.org";

        }
}

$(document).ready(()=> {
    const url=window.location.href.substring(window.location.href.indexOf("?")+1);
    initSessions(); 
    
    run();
    let start_width;
    $("#sidenav").resizable({
        handles: 'e',
        minWidth: 300,
        maxWidth: 1000,
        start(event,ui){
            start_width=$("#sidenav").width();
        },
        resize(event,ui){
            glviewer.center();
            glviewer.translate(($("#sidenav").width()-start_width)/2,0,0,false);
            start_width=$("#sidenav").width();
        }
    });
    $( "#selection_list" ).sortable({
      items: ".selection:not(#spacer)",
      update (event, ui) {
            render(true);
        },
    });// $("#selection_list").accordion();
    $("#selection_list").disableSelection();
    
    initSide(url);
});
