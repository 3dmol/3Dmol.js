        $3Dmol.UI.Form = (function(){
            Form.Color = function(outerControl){
                // var label = $('<div></div>');
                // label.text(outerControl.key);

                var redDot = $('<div></div>');
                redDot.height(10);
                redDot.width(10);
                redDot.css('border-radius', '50%');
                redDot.css('background', 'red');
                redDot.css('margin-right', '3px');

                var blueDot = redDot.clone();
                blueDot.css('background', 'blue');

                var greenDot = redDot.clone();
                greenDot.css('background', 'green');

                var control = this.control = {
                    R : { value : 0, min : 0, max : 255, label : redDot },
                    G : { value : 0, min : 0, max : 255, label : greenDot },
                    B : { value : 0, min : 0, max : 255, label : blueDot },
                };

                var surroundingBox = this.ui = $('<div></div>')
                var boundingBox = $('<div></div>');
                // surroundingBox.append(label);
                surroundingBox.append(boundingBox);
                
                var RValue = new Form.Slider(control.R);
                var GValue = new Form.Slider(control.G);
                var BValue = new Form.Slider(control.B);
                
                var sliders = $('<div></div>');
                sliders.append(RValue.ui, GValue.ui, BValue.ui);
                
                var color = $('<div></div>');
                // color.width(50);
                // color.height(color.width());
                
                boundingBox.append(sliders);
                boundingBox.append(color);
                
                
                // CSS
                // boundingBox.width(400);
                // label.css('width', '100%');
                // boundingBox.css('display', 'flex');
                // boundingBox.css('align-items', 'center');
                
                // RValue.ui.css('display', 'block');
                RValue.slide.css('color', 'red');
                
                // GValue.ui.css('display', 'block');
                GValue.slide.css('color', 'green');
                
                // BValue.ui.css('display', 'block');
                BValue.slide.css('color', 'blue');
                
                // color.css('display', 'inline-block');
                // sliders.css('display', 'inline-block');

                color.height(15);
                // color.width(50);
                color.css('margin-top', '6px');
                color.css('margin-bottom', '6px');
                color.css('border','1px solid grey');
                color.css('border-radius', '500px');

                this.update = function(control){};
                var self = this;
                // Functionality
                function updatePreview(c){
                    var c = `rgb(${control.R.value}, ${control.G.value}, ${control.B.value})`;
                    color.css('background', c);
                    outerControl.value = c;
                    // console.log('rn');
                    self.update(control);
                }

                RValue.update = GValue.update = BValue.update = updatePreview;
                updatePreview();

                this.getValue = function(){
                    return outerControl;
                }

                this.validate = function(){
                    return true;                
                }
            }

            Form.ListInput = function(control, listElements){
                // var label = $('<div></div>');
                // label.text(control.key);

                var surroundingBox = this.ui = $('<div></div>');
                var boundingBox = $('<div></div>');
                var itemList = listElements;
                // surroundingBox.append(label);
                surroundingBox.append(boundingBox);

                var select = $('<select></select>');

                boundingBox.append(select);

                var showAlertBox = this.showAlertBox = true;
                var failMessage = $('<div></div>');
                failMessage.text('Please select some value');
                failMessage.css({
                    'color' : 'crimson',
                    'font-family': 'Arial',
                    'font-weight' : 'bold',
                    'font-size' : '10px'
                });
                failMessage.hide();
                boundingBox.append(failMessage);

                this.update = function(control){
                    console.log("From::Input:update", "Default Update", control);
                }
                
                select.on('click', {parent: this}, (event)=>{
                    control.value = select.children('option:selected').val();
                    event.data.parent.update(control);
                });

                this.getValue = ()=>{
                    return control;
                }

                // this.preventAlertBox = function(){
                //     show
                // }

                this.validate = function(){
                    if(control.value == 'select' || control.value == null){
                        (this.showAlertBox)? failMessage.show() : null ;
                        boundingBox.css({
                            'box-shadow' : '0px 0px 2px red'
                        });
                        return false;
                    }
                    else{
                        failMessage.hide();
                        boundingBox.css({
                            'box-shadow' : 'none'
                        });
                        return true;
                    }
                }

                this.setValue = function(val){
                  if( listElements.indexOf(val) != -1){
                    
                  }
                }

                this.updateList = function(newList){
                  select.empty();

                  var defaultOption = $('<option></option>');
                  defaultOption.text('select');
                  defaultOption.attr('value', 'select');

                  select.append(defaultOption);

                  itemList = newList;
                  itemList.forEach((item)=>{
                    var option = $('<option></option>');
                    option.text(item);
                    option.attr('value', item );
                    select.append(option);
                  });
                }

                this.updateList(itemList);
            }
            
            Form.Input = function(control){
                // var label = $('<div></div>');
                // label.text(control.key);    
                
                var surroundingBox = this.ui = $('<div></div>');
                var boundingBox = $('<div></div>');
                // surroundingBox.append(label);
                surroundingBox.append(boundingBox);

                var validationType = this.validationType = 'text';

                surroundingBox.css({
                    'width' : '100%',
                    'box-sizing':'border-box'
                })

                var input = this.domElement = $('<input type="text">');
                boundingBox.append(input);

                var alertBox = $('<div></div>');
                alertBox.css({
                    'border' : '1px solid darkred',
                    'border-radius' : '3px',
                    'font-family' : 'Arial',
                    'font-size' : '10px',
                    'font-weight' : 'bold',
                    'margin' : '2px',
                    'margin-left' : '4px',
                    'padding' : '2px',
                    'color' : 'darkred',
                    'background' : 'lightcoral'
                });

                var alertMessage = {
                    'invalid-input' : 'Invalid input please check the value entered',
                }
                
                boundingBox.append(alertBox);
                alertBox.hide();

                this.setWidth = function(width){
                    input.width(width);
                }

                this.setWidth(75);

                input.css({
                    'margin-left' : '4px'
                });

                this.update = function(control){
                    console.log("From::Input:update", "Default Update", control);
                }

                // $(document).on('ready', ()=>{
                input.on('change', { parent:this, control: control }, (event)=>{
                    inputString = input.val();
                    
                    if(inputString[inputString.length - 1] == ',') {
                        inputString = inputString.slice(0,-1);
                    }

                    if(validationType == 'range'){
                        control.value = inputString.split(',');
                    }
                    else {
                        control.value = inputString;
                    }
                    
                    // calling update function 
                    event.data.parent.update(control);
                });

                var selectedText = null;
                
                input.on('select', (e)=>{
                    console.log('Selection in input', e, input, e.target.selectionStart, e.target.selectionEnd, input.val().substring(e.target.selectionStart, e.target.selectionEnd) );
                    
                    selectedText = input.val().substring(e.target.selectionStart, e.target.selectionEnd);
                });
                
                
                this.getValue = ()=>{
                    // console.log("Form Inut", control);
                    return control;
                }
                
                var error = this.error = function(msg){
                    alertBox.show();
                    alertBox.text(msg)
                }

                this.setValue = function(val){
                    input.val(val)
                }

                

                function checkInputFloat(){
                    var inputString = input.val();
                    
                    var dots = inputString.match(/\./g) || [];
                    var checkString = inputString.replaceAll(/\./g, '').replaceAll(/[0-9]/g, '');

                    console.log('Check Float', dots, checkString)
                    if(dots.length > 1){ return false };

                    if(checkString != '') return false;

                    if(isNaN(parseFloat(inputString))){
                        return false;
                    }
                    else {
                        return true;
                    }
                }

                function checkInputNumber(){
                    var inputString = input.val();

                    var checkString = inputString.replaceAll(/[0-9]/g, '');

                    if(checkString != '') return false;

                    if(isNaN(parseInt(inputString))){
                        return false;
                    }
                    else {
                        return true;
                    }
                }

                // Parse Input Range Functions

                // Checks only number, comma and hyphen present
                function checkRangeTokens(inputString){
                    var finalString = inputString.replaceAll(',','').replaceAll('-','').replaceAll(/[0-9]/g, '').replaceAll(' ', '');

                    console.log('Check string in input range', finalString );

                    if(finalString == '')
                        return true;
                    else 
                        return false;
                }

                function checkList(inputString, submit=false){
                    inputString = inputString.replaceAll(' ', '');
                    
                    if(inputString[inputString.length - 1] == ',') {
                        inputString = inputString.slice(0,-1);
                    }


                    var rangeList = inputString.split(',');

                    // If dublicate comma return false;
                    if(/,,/g.exec(inputString)) return false;

                    // If first element not a number return false;
                    if(isNaN(parseInt(rangeList[0]))) return false;

                    var validRangeList = rangeList.map((rangeInput)=>{
                        return checkRangeInput(rangeInput);
                    });

                    console.log("Checking Range List", rangeList, validRangeList);
                    return validRangeList.find( (e)=>{ return e == false}) == undefined ? true : false; 
                }

                function checkRangeInput(inputString, submit=false){  
                    var rangeInputs = inputString.split('-');
                    if(rangeInputs.length > 2){
                        return false;
                    }
                    else {
                        if(rangeInputs.length == 0){
                            return true;
                        }
                        else if (rangeInputs.length == 1){
                            if(isNaN(parseInt(rangeInputs[0])))
                                return false;
                            else 
                                return true;
                        }
                        else if (rangeInputs.length == 2){
                            if(isNaN(parseInt(rangeInputs[0])) || isNaN(parseInt(rangeInputs[1])))
                                return false;
                            else 
                                return true;
                        }
                        else 
                            return false;
                    }
                }

                var checkInput = this.checkInput = function(){
                    var inputString = input.val();
                    console.log("CheckInput", this.validationType, validationType)
                    
                    if(validationType == 'number'){
                        if(checkInputNumber()){
                            alertBox.hide();
                            return true;
                        }else{
                            error(alertMessage['invalid-input']);
                            return false;
                        }
                    }
                    else if(validationType == 'float'){
                        if(checkInputFloat()){
                            alertBox.hide();
                            return true;
                        }else{
                            error(alertMessage['invalid-input']);
                            return false;
                        }
                    }
                    else if(validationType == 'range') {
                        if(checkRangeTokens(inputString)){
                            if(checkList(inputString)){
                                alertBox.hide();
                                return true;
                            }
                            else {
                                error(alertMessage['invalid-input']);
                                return false;   
                            }
                        }
                        else{
                            error(alertMessage['invalid-input']);
                            return false;
                        }

                    }
                    else {
                        return true;
                    }
                }

                this.validateOnlyNumber = function(floatType=false){
                    if(floatType){
                        validationType = 'float';
                    } else {
                        validationType = 'number';
                    }

                    input.on('keydown keyup paste cut', function(event){
                        console.log('Validation Type', validationType);
                        checkInput();
                    });
                }


                this.validateInputRange = function(){
                    validationType = 'range';

                    input.on('keydown keyup paste cut', ()=>{
                        console.log('Validation Type', validationType);
                        checkInput();
                    });

                }

                this.isEmpty = function(){
                    if(control.value == ""){
                        return true;
                    }
                }

                this.validate = function(){
                    if( (control.active == true && control.value != null && control.value != "" && checkInput()) || (control.active == false)){
                        input.css('box-shadow', 'none');
                        console.log("Form::Input:validate", 'success', control);
                        return true
                    }
                    else {
                        input.css('box-shadow', '0px 0px 2px red');
                        console.log("Form::Input:validate", 'failure', control);
                        return false;
                    }
                }

            }
            
            Form.Checkbox = function(control){
                var label = $('<div></div>');
                label.text(control.key); 

                var surroundingBox = this.ui = $('<div></div>');
                var boundingBox = $('<div></div>');
                surroundingBox.append(boundingBox);
                surroundingBox.append(label);

                var checkbox = $('<input type="checkbox" />');
                boundingBox.append(checkbox);

                this.click = ()=>{};

                this.update = function(control){
                    console.log("From::Input:Checkbox", "Default Update", control);
                }

                this.getValue = ()=>{
                    return control;
                }

                checkbox.on('click',{ parent: this}, (event)=>{
                    control.value = checkbox.prop('checked');
                    event.data.parent.update(control);
                });

                // CSS
                label.css('display', 'inline-block');
                boundingBox.css('display', 'inline-block')

                this.validate = function(){
                    return true;
                }

                this.setValue = function(val){
                    checkbox.prop('checked', val);
                }
            }

            Form.Slider = function(control, style = {}){
                // var label = $('<div></div>');
                // label.text(control.key); 

                var surroundingBox = this.ui = $('<div></div>');
                // surroundingBox.append(label);
                
                var boundingBox = $('<div></div>');
                surroundingBox.append(boundingBox);
                
                boundingBox.css('display', 'flex');
                var slide = this.slide = $('<input type="range">');
                slide.css('width' , '100%');
                
                var min = control.min || 0;
                var max = control.max || 100;
                var step = control.step || 1;
                var defaultValue = control.default || min;
                var labelContent = control.label || '';

                var label = $('<div></div>');
                label.append(labelContent);
                boundingBox.append(label);

                slide.attr('min', min);
                slide.attr('max', max);
                slide.attr('step', step);
                console.log('step::', step);
                slide.attr('value', defaultValue);
                control.value = defaultValue;
                boundingBox.append(slide);

                var setValue = false;
                
                this.update = function(control){
                    console.log("Form::Slider:update", "Default Update", control);
                };

                this.getValue = ()=>{
                    return control;
                }

                slide.on('mousedown', ()=>{
                    setValue = true;
                });

                slide.on('mousemove', { parent: this } ,(event)=>{
                    if(setValue){
                        control.value = slide.val();
                        event.data.parent.update(control);
                        // console.log('working');
                    }
                });

                slide.on('mouseup', ()=>{
                    setValue = false;
                });

                // CSS
                boundingBox.css('align-items', 'center');
                boundingBox.height('21px');
                // boundingBox.css('border-radius', '2px');
                // label.css('line-height', '21px');
                slide.css('padding', '0px');
                slide.css('margin', '0px');
                
                this.validate = function(){
                    return true;
                }

                this.setValue = function(val){
                    slider.val(val);
                    control.value = slider.val();
                }


            }

            Form.EmptyElement = function(control){
                this.ui = $('<div></div>');

                var update = ()=>{};
                
                this.onUpdate = (callback)=>{
                    update = callback;
                }

                this.getValue = ()=>{
                    return control;
                }

                this.validate = function(){
                    return true;
                }
            }

            // mainControl param will be used to take in specName
            // in the form of key 
            // type will be 'form'
            // active will be used to activate deactivate form if more than one form
            function Form(specs, mainControl ){
                specs = specs || {};
                var boundingBox = this.ui = $('<div></div>');

                var heading = $('<div></div>');
                heading.text(mainControl.key);

                // Styling heading 
                heading.css({
                    'border-bottom':'1px solid black',
                    'font-family':'Arial',
                    'font-size':'14px',
                    'font-weight':'bold',
                    'padding-top':'2px',
                    'padding-bottom':'4px'
                });

                boundingBox.append(heading);
                boundingBox.addClass('form');
                
                var inputs = this.inputs = [];
                // body.append(boundingBox);
                
                var keys = Object.keys(specs);
                keys.forEach((key)=>{
                    if(specs[key].gui){
                        var prop = new Property(key, specs[key].type);
                        inputs.push(prop);
                        boundingBox.append(prop.ui);
                    }

                    // console.log('Checking Specs', prop.placeholder, key, specs[key]);
                });
                var self = this;

                this.update = function(){
                }

                var update = (control)=>{
                    console.log("Updated Form Values", this.getValue());
                    this.update();
                };


                inputs.forEach((input)=>{
                    input.update = update;
                })
                
                this.getValue = function(){
                    mainControl.value = {};
                    // console.log("Form Inputs", inputs);
                    inputs.forEach((input)=>{
                        var inputValue = input.getValue();
                        // console.log("Current Input Value", inputValue)
                        if(inputValue.active){
                            mainControl.value[inputValue.key] = inputValue.value;
                        }
                    });

                    return mainControl;
                }
                
                var updateValues = function(inputControl){
                    mainControl.value[inputControl.key] = control.value;
                    update(mainControl);     
                }

                this.validate = function(){
                    var validations = inputs.map((i)=>{
                        // console.log("Checking inputs", i)
                        if(i.active.getValue().value){
                            return i.placeholder.validate();
                        }
                        else {
                            return true;
                        }
                    });

                    console.log('Checking Validation', validations);

                    if(validations.find( e => e == false) == undefined )
                        return true;
                    else {
                        return false;
                    }

                }

                this.setValue = function(val){
                    var keys = Object.keys(val);
                    for(var i = 0; i < keys.length; i++){
                        var input = inputs.find((e)=>{
                            if(e.control.key == keys[i])
                                return e;
                        });

                        e.placeholder.setValue(val[key[i]]);
                    }

                    console.log('Setting Key Value', val, inputs);
                }
                
                function Property(key, type) {
                    var control = this.control = { value : null, type : type, key : key, active: false };
                    var boundingBox = this.ui = $('<div></div>');
                    this.placeholder = { ui : $('<div></div>') }; // default value for ui element 
                    this.active = new Form.Checkbox({value: false, key:key});


                    if(specs[key].type == 'string' || specs[key].type == 'element'){
                        this.placeholder = new Form.Input(control);
                        this.placeholder.ui.attr('type', 'text');
                    }
                    else if(specs[key].type == 'number'){
                        console.log('Form::Property', specs[key]);
                        var slider = false;
                        
                        if(specs[key].min != undefined && specs[key].max != undefined && specs[key].default != undefined){
                            slider = true;
                        }

                        if(slider){
                        // if( specs[key].min && spec[key].max){
                            control.min = specs[key].min;
                            control.max = specs[key].max;
                            control.default = specs[key].default;
                            control.step = specs[key].step || ((control.max - control.max )/ 1000);
                            this.placeholder = new Form.Slider(control);
                        }
                        else{
                            this.placeholder = new Form.Input(control);
                            this.placeholder.ui.attr('type', 'text');
                            this.placeholder.validateOnlyNumber(specs[key].floatType);
                        }
                    }
                    else if(specs[key].type == 'array_range'){
                        this.placeholder = new Form.Input(control);
                        this.placeholder.ui.attr('type', 'text');
                        this.placeholder.validateInputRange();
                    }
                    else if(specs[key].type == 'color'){
                        this.placeholder = new Form.Color(control);

                    }else if(specs[key].type == 'boolean'){
                        this.placeholder = new Form.Checkbox(control);

                    }else if(specs[key].type == 'properties'){
                        this.placeholder = new Form.Input(control);
                        this.placeholder.ui.attr('type', 'text');

                    }else if(specs[key].type == 'colorscheme'){
                        this.placeholder = new Form.ListInput(control, Object.keys($3Dmol.builtinColorSchemes));
                                this.placeholder.ui.attr('type', 'text');
                            
                    }else if(specs[key].type == undefined){
                        if( specs[key].validItems) {
                                    this.placeholder = new Form.ListInput(control, specs[key].validItems);
                        }
                    
                    }else if(specs[key].type == 'form'){
                                this.placeholder = new Form(specs[key].validItems, control);
                                this.placeholder.ui.append($('<hr />'));
                    }
                    else {
                        this.placeholder = new Form.EmptyElement(control);
                                // return new Form.EmptyElement(control);
                    }

                    // console.log("Checking placeholder", this.placeholder);
                    
                    this.getValue = function(){
                        // console.log("Property Input", key ,this.placeholder);
                        if(this.placeholder.getValue)
                            return this.placeholder.getValue();
                        else
                            return null;
                    }


                    // Adding active control for the property
                    var placeholder = this.placeholder;
                    // console.log("Property", placeholder);
                    placeholder.ui.hide();
                    // this.active.ui.width(200);


                    if(type !='boolean') {
                        boundingBox.append(this.active.ui);
                        this.active.update = function(c){
                            (c.value)? placeholder.ui.show() : placeholder.ui.hide();
                            control.active = c.value;
                        }
                    }
                    
                    boundingBox.append(this.placeholder.ui);

                    if(this.placeholder.onUpdate)
                        this.placeholder.onUpdate(updateValues);
                }
                

            }

            
            
            return Form;
        })();