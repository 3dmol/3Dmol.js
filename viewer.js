
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

    const attributeName = $('<select>',{
        class:'attributeName',
    }).appendTo(attribute);
    const objType = type;

    $.each(validNames,(key,value) => {
        if(value.gui){
            attributeName.append($("<option>").attr('value',key).text(key));
        }
    });

    attributeName.val(name.toString())
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
        attributeName.val(index)
    }

    // delete button
    const deleteSelection = $("<span/>",{
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
    let attributeValue;
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

        const attributeValue = $('<select/>',{
            class:'attributeValue',
        }).appendTo(attribute);

        $.each(validItemsValue,(key,value) => {
            attributeValue.append($("<option>").attr('value',value).text(value));
        });

        attributeValue.val(value.toString());
        if(value === ""){
            attributeValue.val(validItemsValue[0])
        } 
    }else{
        if(value === "")
            value = validNames[name].default
        attributeValue = $('<input/>',{
            class:'attributeValue',
            value,
        }).appendTo(attribute);
    }
    attributeName.change(()=> {
        let validItemsValue;
        const {type} = validNames[attributeName.val()]
        if(type==="boolean"){
            validItemsValue = ["false","true"];
        }else if(type === "colorscheme"){
            validItemsValue =  $3Dmol.GLModel.validColorschemeSpecs;
        }else if(type === "color"){
            validItemsValue =  $3Dmol.GLModel.validColorSpecs;
        }else if(type === undefined){
            validItemsValue = validNames[name].validItems;
        }
        const defa = validNames[attributeName.val()].default;
        let val;
        if(validItemsValue !== undefined){
            val = validItemsValue[0];
        }else{
            val = defa
        }
        if(attributeValue.children()[0]!== undefined)
            attributeValue.children()[0].value = val;
        else
            attributeValue.val(val);
        render(objType === "surface");
    });
    attributeValue.change(()=> {
        render(objType === "surface");
    });

    if(name!=="" &&attributeValue.prop("tagName") === "INPUT" && validNames[name].type ==="number"){
        // eslint-disable-next-line no-unused-expressions
        validNames[name].type ==="number" // can this be commented out?
        attributeValue.attr("type","number")
        attributeValue.attr("step",validNames[name].step)
        attributeValue.addClass("spinner")
        const {max} = validNames[name];
        const {min} = validNames[name];
        if(max !== undefined)
            attributeValue.attr("max",max);
        if(min !== undefined)
            attributeValue.attr("min",min);   
    }
    return attribute;
}

const createOtherModelSpec = function(spec,type,selectionIndex){
    const attributes = $('<ul/>',{
        "class":`${type.toLowerCase()}_attributes`,
    });
    for(const attributeIndex in spec){
        const attribute=createAttribute(attributeIndex,spec[attributeIndex],{type,index:selectionIndex})
        if(attribute !== undefined)
            attribute.appendTo(attributes);
    }

    const addAttribute = $('<button/>',{
        "class":"addAttribute",
        "text":"Add Attribute",
        "data-index":selectionIndex,
        "data-type":type,
        "click":function(){addOtherAttribute(this)},
    }).appendTo(attributes);

    return attributes;
}

const createStyleSpec = function(styleSpecObject,styleSpecType,modelSpecType,selectionIndex){
    const styleSpec=$('<li/>',{
        "class":"styleSpec",
    });

    const validNames=$3Dmol.GLModel.validAtomStyleSpecs;

    const styleSpecName = $('<select>',{
        class:'styleSpecName',
    }).appendTo(styleSpec);

    styleSpecName.change(()=> {
        const obj = query.selections[selectionIndex].style[styleSpecType];
        for(const i in obj){
            if(!validNames[styleSpecName.val()].validItems.hasOwnProperty(i)){
                delete query.selections[selectionIndex].style[styleSpecType][i];
            }
        }
        query.selections[selectionIndex].style[styleSpecName.val()]=query.selections[selectionIndex].style[styleSpecType];
        delete query.selections[selectionIndex].style[styleSpecType];
        buildHTMLTree(query)
        render();
    });

    $.each(validNames,(key,value) => {
        if(value.gui){
            styleSpecName.append($("<option>").attr('value',key).text(key));
        }
    });
    
    styleSpecName.val(styleSpecType.toString())
    if(styleSpecType === ""){
        const list = query.selections[selectionIndex].style;
        let index=0
        for(const i in validNames){
            if(validNames[i].gui && list[i] === undefined){
                index = i;
                break;
            }
        }
        if(index === 0)
            return;
        styleSpecName.val(index)
    }

    const deleteSelection = $("<span/>",{
        html:"&#x2715;",
        class:"delete_style_spec",
        "data-index":selectionIndex,
        "data-type":modelSpecType,
        "data-attr":styleSpecType,
        "click":function(){deleteStyleSpec(this)},
    }).appendTo(styleSpec); 

    const styleSpecAttributes = $('<ul/>',{
        class:'styleSpecAttributes',
    }).appendTo(styleSpec);

    for(const attributeIndex in styleSpecObject){
        const attribute = createAttribute(attributeIndex,styleSpecObject[attributeIndex],{type:styleSpecType,index:selectionIndex})
        if(attribute !== undefined)
            attribute.appendTo(styleSpecAttributes);
    }

    const addAttribute = $('<button/>',{
        "class":"addAttribute",
        "text":"Add Attribute",
        "data-index":selectionIndex,
        "data-type":modelSpecType,
        "data-styletype":styleSpecType,
        "click":function(){addAttribute(this)},
    }).appendTo(styleSpec);
                        
    return styleSpec;
}

const createStyle = function(modelSpecObject,modelSpecType,selectionIndex){
    const style=$('<span/>',{
        "class":"style",
    }); 

    const styleSpecs = $('<ul/>',{
        "class":'styleSpecs',
    }).appendTo(style);
                        
    for(const attributeIndex in modelSpecObject){
        const spec = createStyleSpec(modelSpecObject[attributeIndex],attributeIndex,modelSpecType,selectionIndex)
        if(spec!== undefined)
            spec.appendTo(styleSpecs);
    }

    const addStyleSpec = $('<button/>',{
        "class":"addStyleSpec",
        "text":"Add Style Spec",
        "data-index":selectionIndex,
        "data-type":modelSpecType,
        "click":function(){addStyleSpec(this)},
    }).appendTo(style);

    return style; 
}

const validNames = {
    "style":"Style",
    "surface":"Surface",
    "labelres":"LabelRes",
}

const createModelSpecification = function(modelSpecType,modelSpecObject,selectionIndex){
    let modelSpecification = null;
    if(modelSpecType==="style"){
        modelSpecification = createStyle(modelSpecObject,modelSpecType,selectionIndex)
    }else if(modelSpecType==="surface"){
        modelSpecification = createOtherModelSpec(modelSpecObject,"Surface",selectionIndex)
    }else if(modelSpecType==="labelres"){
        modelSpecification = createOtherModelSpec(modelSpecObject,"LabelRes",selectionIndex)
    }             

    return modelSpecification;
}
// this function creates the selection object
const createSelection = function(spec,object,index,type){
    // creates container
    const selection = $("<li/>",{
        class:"selection"
    });

    const createHeader = function(){
        const selectionType = $('<p>',{
            class:'selectionType',
            text:validNames[type],
        }).appendTo(selection);

        // add together sub selections
        const attributePairs =[];
        for(const subselection in spec){
            let obj=spec[subselection];
            if(typeof(obj) === 'object' && Object.keys(obj).length === 0)
                obj = ""; // empty object
            attributePairs.push(`${subselection}:${obj}`);
        }

        let modifier=attributePairs.join(";");
        if(modifier === "")
            modifier = "all"
        const selectionSpec=$('<input/>', {
            class:'selectionSpec',
            value:modifier,
        }).appendTo(selection); 

        selectionSpec.change(()=> {
            render(type === "surface");
        })
    }

     // delete button
    const deleteSelection = $("<div/>",{
        html:"&#x2715;",
        class:"deleteSelection",
        "data-index":index,
        "data-type":"",
        "click":function(){deleteSelection(this);}
    }).appendTo(selection); 

    createHeader()    
    // check if style exists and if so create the object
    const ret = createModelSpecification(type,object, index);

    deleteSelection.attr("data-type",type);
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
        if(prevType !== val){
            
            render(true);
            run();
        }
        prevType = val
    })

    $("#model_input").attr("value",query.file.path);
    $("#model_input").change(()=> {
        const val =  $("#model_input").val().toUpperCase();
        if(prevIn !== val){
            if(val.match(/^[1-9][A-Za-z0-9]{3}$/) || $("#model_type").val().toLowerCase()!== "pdb"){
                render(true);
                run();
                const width = $("#sidenav").width();
            }else{
                if(prevIn!== val)
                    alert("Invalid PDB")
            }
        }
        prevIn = val;
    })
    const arr=[]
    // loops through selections and creates a selection tree
    for(const selectionIndex in query.selections){
        const selectionObject = query.selections[selectionIndex];
        const aug = augmentSelection(selectionObject);

        if(selectionObject.style !== undefined){
            arr.push(createSelection(aug,selectionObject.style,selectionIndex,"style"));
        }
        if(selectionObject.surface !== undefined){
            arr.push(createSelection(aug,selectionObject.surface,selectionIndex,"surface"))
        }
        if(selectionObject.labelres !== undefined){
            arr.push(createSelection(aug,selectionObject.labelres,selectionIndex,"labelres"))
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
        $.each(object, (subStyle,subStyleObject)=> {
            let string="";
            string+=subStyle;
            if(Object.size(subStyleObject)!==0)
                string+=":";
            const assignments =[]
            $.each(subStyleObject, (key,value)=> {
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
        }else if(type === "style" || type==="surface" || type === "labelres"){
            if(currentSelection == null){
                currentSelection = {}
                query.selections.push(currentSelection)
            }
            currentSelection[type] = object;
        // eslint-disable-next-line no-self-compare
        }else if(type === type){
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
            object[$(otherList[li]).children(".attributeName")[0].value]=$(otherList[li]).children(".attributeValue")[0].value
        });
        return object;
    }

    const updateStyle = function(styl){
        const object={};
        let list = $(styl).children(".styleSpecs");
        list = $(list).children(".styleSpec")
        list.each((li)=> {
            const subtype=$(list[li]).children(".styleSpecName")[0].value;
            object[subtype]={};
            let otherList =$(list[li]).children(".styleSpecAttributes")[0];
            otherList=$(otherList).children(".attribute")
            otherList.each((li)=> {
                const tag=object[subtype][$(otherList[li]).children(".attributeName")[0].value]=$(otherList[li]).children(".attributeValue")[0].tagName
                object[subtype][$(otherList[li]).children(".attributeName")[0].value]=$(otherList[li]).children(".attributeValue")[0].value;
            });
        });
        return object;
    }

    const updateSelectionElements = function(selectionString){
        return $3Dmol.specStringToObject(selectionString);
    }

    function arraysEqual(a, b) {
        if (a === b) return true;
        if (a == null || b == null) return false;
        if (a.length !== b.length) return false;

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
            if(obj2[key]===undefined || obj2[key] !== obj1[key])
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
        if(listItems.hasOwnProperty(index) && listItems[index].id!=="spacer"){
            const getSubObject = function(){
                const attr = $(value);
                const attribute=attr[0]
                const type=$(attribute).children()[1].innerHTML.toLowerCase()

                if(type==="style"){
                    const style =updateStyle($(attribute).children(".style")[0])
                    return {"style":style}
                }if(type==="surface"){
                    const surface = updateOther($(attribute).children(".surface_attributes")[0])
                    return {"surface":surface} 
                }if(type === "labelres"){
                    const labelres = updateOther($(attribute).children(".labelres_attributes")[0])
                    return {"labelres":labelres} 
                }
            }

            const val = getSubObject();
            const selectionSpec = $(listItems[index]).children(".selectionSpec")[0].value;
            const selection = updateSelectionElements(selectionSpec);
            const extended = combine(selection,val)
            selects.push(extended)
        }
    });

    query.selections=selects;
}

let query = urlToQuery(window.location.search.substring(1));
// this function compresses the html object back into a url
let render = function(surfaceEdited){
    surfaceEdited = surfaceEdited === undefined ? false : surfaceEdited;
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
    const surface  = type === "surface"
    if(type === "style")      
        query.selections.push({"style":{line:{}}})
    else if(type === "surface")
        query.selections.push({"surface":{}})
    else if(type === "labelres")
        query.selections.push({"labelres":{}})
    buildHTMLTree(query);
    render(surface);
}

const deleteSelection = function(spec){
    delete query.selections[spec.dataset.index][spec.dataset.type];
    if(query.selections[spec.dataset.index].surface === undefined && query.selections[spec.dataset.index].style === undefined && query.selections[spec.dataset.index].labelres === undefined)
        delete query.selections[spec.dataset.index]
    
    buildHTMLTree(query);
    render(spec.dataset.type === "surface");
}

const addModelSpec = function(type,selection){
    let currentSelection;
    currentSelection = query.selections[selection.dataset.index]
    
    if(type === "style" || type === "surface" || type === "labelres"){
        if(currentSelection[type]==null)
            currentSelection[type]={};
        else
            console.err(`${type} already defined for selection`);// TODO error
                                                                // handling
    }
    
    buildHTMLTree(query);
    render();
}

const addStyleSpec = function(modelSpec){
    const defaultKey = "";
    const defaultValue = {};
    query.selections[modelSpec.dataset.index][modelSpec.dataset.type][defaultKey]=defaultValue;
    
    buildHTMLTree(query);
    render();
}

let deleteStyleSpec = function(spec){
    delete query.selections[spec.dataset.index][spec.dataset.type][spec.dataset.attr]
    
    buildHTMLTree(query);
    render();
}

let addOtherAttribute= function(spec){
    const defaultKey = "";
    const defaultValue = "";
    query.selections[spec.dataset.index][spec.dataset.type.toLowerCase()][defaultKey]=defaultValue;
    
    buildHTMLTree(query);
    render();
}

let deleteOtherAttribute = function(spec){
    delete query.selections[spec.dataset.index][spec.dataset.type][spec.dataset.attr]
    
    buildHTMLTree(query);
    render(spec.dataset.type === "surface");
}

const addAttribute = function(styleSpec){
    const defaultKey = "";
    const defaultValue = "";
    query.selections[styleSpec.dataset.index][styleSpec.dataset.type][styleSpec.dataset.styletype][defaultKey]=defaultValue;

    buildHTMLTree(query);
    render();
}

let deleteStyleAttribute = function(spec){
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
let prevIn = $("#model_input").val();
let prevType = $("#model_type").val();
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

let glviewer = null;
// http://localhost/$3Dmol/viewer.html?pdb=1ycr&style=cartoon&addstyle=line&select=chain~A&colorbyelement=whiteCarbon&style=surface,opacity~.8&select=chain~B&addstyle=stick&select=chain~B,resn~TRP&style=sphere
// Process commands, in order, and run on viewer (array of strings split on '&')
let runcmds = function(cmds, viewer,renderSurface) {
            console.log("rendering")
    renderSurface = renderSurface === undefined ? true : renderSurface;
    if(renderSurface)
        viewer.removeAllSurfaces();
    viewer.removeAllLabels();
    let currentsel = {};

    for (let i = 0; i < cmds.length; i++) {
        const kv = cmds[i].split('=');
        const cmdname = kv[0];
        const cmdobj = $3Dmol.specStringToObject(kv[1]);

        if (cmdname === 'select')
            currentsel = cmdobj;
        else if (cmdname === 'surface' && renderSurface){
            viewer.addSurface($3Dmol.SurfaceType.VDW, cmdobj, currentsel,
                    currentsel);
        } else if (cmdname === 'style'){
            viewer.setStyle(currentsel, cmdobj);
        } else if (cmdname === 'addstyle'){
            viewer.addStyle(currentsel, cmdobj);
        } else if (cmdname === 'labelres'){
            viewer.addResLabels(currentsel, cmdobj);
        } else if (cmdname === 'colorbyelement'){
            if (typeof ($3Dmol.elementColors[cmdobj.colorscheme]) != "undefined")
                viewer.setColorByElement(currentsel,
                        $3Dmol.elementColors[cmdobj.colorscheme]);
        } else if (cmdname === 'colorbyproperty'){
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

            if (src === 'session' || src === 'SESSION') {
                // join a session
                joinSession(data);
            }
            if (src === 'pdb') {
                console.log(data)
                data = data.toUpperCase();
                if (!data.match(/^[1-9][A-Za-z0-9]{3}$/)) {
                    return;
                }
                data = `http://files.rcsb.org/view/${  data
                         }.pdb`;
                type = "pdb";
            } if (src === 'cif') {
                data = data.toUpperCase();
                if (!data.match(/^[1-9][A-Za-z0-9]{3}$/)) {
                    return;
                }
                data = `http://files.rcsb.org/view/${  data
                         }.cif`;
                type = "cif";
            } else if (src === 'cid') {
                type = "sdf";
                data = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${
                         data  }/SDF?record_type=3d`;
            } else if (src === 'mmtf') {
                data = data.toUpperCase();
                data = `http://mmtf.rcsb.org/full/${  data  }.mmtf`;
                type = 'mmtf';
            } else { // url
                // try to extract extension
                type = data.substr(data.lastIndexOf('.') + 1);
                if(type === 'gz') {
                    const base = data.substr(0,data.lastIndexOf('.'));
                    type = `${base.substr(base.lastIndexOf('.'))  }.gz`;
                }
            }
            if (cmds[0] && cmds[0].indexOf('type=') === 0) {
                type = cmds[0].split('=')[1];
            }

            const start = new Date();

            if (/\.gz$/.test(data) || type === 'mmtf') { // binary
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
                        if(src === 'pdb' && (ret.search("We're sorry, but the requested") >= 0 || ret === "")) {
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
    let startWidth;
    $("#sidenav").resizable({
        handles: 'e',
        minWidth: 300,
        maxWidth: 1000,
        start(event,ui){
            startWidth=$("#sidenav").width();
        },
        resize(event,ui){
            glviewer.center();
            glviewer.translate(($("#sidenav").width()-startWidth)/2,0,0,false);
            startWidth=$("#sidenav").width();
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
