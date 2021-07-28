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
                    G : { value : 0, min : 0, max : 255, label : blueDot },
                    B : { value : 0, min : 0, max : 255, label : greenDot },
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
                boundingBox.css('display', 'flex');
                boundingBox.css('align-items', 'center');
                
                // RValue.ui.css('display', 'block');
                RValue.slide.css('color', 'red');
                
                // GValue.ui.css('display', 'block');
                GValue.slide.css('color', 'green');
                
                // BValue.ui.css('display', 'block');
                BValue.slide.css('color', 'blue');
                
                color.css('display', 'inline-block');
                sliders.css('display', 'inline-block');

                color.height(50);
                color.width(50);
                color.css('margin-left', '6px');
                color.css('border','1px solid black');
                color.css('border-radius', '2px');

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

                this.validate = function(){
                    if(control.value == 'default')
                        return false;
                    else 
                        return true;
                }

                this.setValue = function(val){
                  if( listElements.indexOf(val) != -1){
                    
                  }
                }

                this.updateList = function(newList){
                  select.empty();

                  var defaultOption = $('<option></option>');
                  defaultOption.text('select');
                  defaultOption.attr('value', 'default');

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

                var input = this.domElement = $('<input type="text">');
                boundingBox.append(input);

                this.update = function(control){
                    console.log("From::Input:update", "Default Update", control);
                }

                // $(document).on('ready', ()=>{
                input.on('change', { parent:this, control: control }, (event)=>{
                    control.value = input.val();
                    
                    // calling update function 
                    event.data.parent.update(control);
                });
                // });
                
                
                this.getValue = ()=>{
                    // console.log("Form Inut", control);
                    return control;
                }
                
                this.error = function(){
                    $(document).trigger('Error', ["Incorrect Input", "Please enter the input as per the direction"]);  
                }

                this.setValue = function(val){
                    input.val(val)
                }

                this.validateOnlyNumber = function(){
                    var decimalEntered = false;

                    input.on('keydown', function(event){
                        event.preventDefault();

                        if((event.key.charCodeAt(0) >= "0".charCodeAt(0) && event.key.charCodeAt(0) <= "9".charCodeAt(0)) || event.key == ".") {

                            if(decimalEntered && event.key == ".") {
                                console.log('Again Entering decimal');
                                return;
                            }
                            $(this).val( $(this).val() + event.key);
                            // event.preventDefault();
                            console.log("Input Number Only", event.key.charCodeAt(), event.key);
                            if(event.key == ".")
                                decimalEntered = true;
                        }
                        else if(event.key == "Backspace"){
                            var val = $(this).val();
                            var toRemove = val.slice(-1);
                            if(toRemove == '.')
                                decimalEntered = false;
                            $(this).val(val.slice(0,-1));
                        }
                        else{
                            console.log('Not a number in the input');
                        }

                        //  var text = $(this).val().replace(/[^\d]+/, "");
                        // $(this).val(text);
                        input.trigger('change')
                        console.log("Text", control.value);
                    });
                }

                this.validateInputRange = function(){
                    var dashEntered = false;
                    var commaEntered = false;
                    input.on('keydown', function(event){
                        event.preventDefault();
                        
                        if( $(this).val().length < 1 && (event.key == '-' || event.key == ',')){
                            console.log('Cannnot Enter "-" or "," at the beginning');
                            return;
                        }

                        if((event.key.charCodeAt(0) >= "0".charCodeAt(0) && event.key.charCodeAt(0) <= "9".charCodeAt(0)) || event.key == '-' || event.key == ',') {

                            if(dashEntered && event.key == "-") {
                                console.log('Again dash entered');
                                return;
                            }
                            else if(dashEntered && event.key != "-" && event.key != ","){
                                dashEntered = false;
                            }
                            else if(dashEntered && event.key == ",") {
                                return;
                            }

                            if(commaEntered && event.key == ",") {
                                console.log('Again comma entered');
                                return;
                            }
                            else if(commaEntered && event.key != "," && event.key != "-"){
                                commaEntered = false;
                            }
                            else if(commaEntered && event.key == "-") {
                                return;
                            }

                            $(this).val( $(this).val() + event.key);
                            // event.preventDefault();
                            console.log("Input Number Only", event.key.charCodeAt(), event.key);
                            if(event.key == "-")
                                dashEntered = true;

                            if(event.key == ",")
                                commaEntered = true;
                        }
                        else if(event.key == "Backspace"){
                            var val = $(this).val();
                            var toRemove = val.slice(-1);
                            var beforeToRemove = val.slice(-2);

                            if(beforeToRemove == '-')
                                dashEntered = true;
                            else if (toRemove == '-')
                                dashEntered = false;

                            if(beforeToRemove == ',')
                                commaEntered = true;
                            else if (toRemove == ',')
                                commaEntered = false;


                            $(this).val(val.slice(0,-1));
                        }
                        else{
                            console.log('Not a number in the input');
                        }

                        //  var text = $(this).val().replace(/[^\d]+/, "");
                        // $(this).val(text);
                        input.trigger('change')
                        console.log("Text", control.value);

                    });
                }

                this.isEmpty = function(){
                    if(control.value == ""){
                        return true;
                    }
                }

                this.validate = function(){
                    if((control.active == true && (control.value != null && control.value != "")) || (control.active == false)){
                        input.css('box-shadonw', 'none');
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

            Form.Slider = function(control){
                // var label = $('<div></div>');
                // label.text(control.key); 

                var surroundingBox = this.ui = $('<div></div>');
                // surroundingBox.append(label);
                
                var boundingBox = $('<div></div>');
                surroundingBox.append(boundingBox);
                
                boundingBox.css('display', 'flex');
                var slide = this.slide = $('<input type="range">');
                

                var min = control.min || 0;
                var max = control.max || 100;
                var step = control.step || 1;
                var defaultValue = control.value || 0;
                var labelContent = control.label || '';

                var label = $('<div></div>');
                label.append(labelContent);
                boundingBox.append(label);

                slide.attr('min', min);
                slide.attr('max', max);
                slide.attr('step', step);
                slide.attr('value', defaultValue);
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

                boundingBox.append($(`<p><b>${mainControl.key}</b></p><div style="border-top:1px solid black"></div>`));
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
                        return i.placeholder.validate();
                    });

                    if(validations.find( e => e == false) == undefined )
                        return true;
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
                        if(false){
                        // if( specs[key].min && spec[key].max){
                            control.min = specs[key].min;
                            control.max = specs[key].max;
                            control.default = specs[key].default;
                            this.placeholder = new Form.Slider(control);
                        }
                        else{
                            this.placeholder = new Form.Input(control);
                            this.placeholder.ui.attr('type', 'text');
                            this.placeholder.validateOnlyNumber();
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
                    this.active.ui.width(200);


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