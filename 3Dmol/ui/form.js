    $3Dmol.UI.form = (function(){
        function form(specs){
            specs = specs || {};

            // Test area
            body = $('body');
            // var changingColor = new Color();
            // body.append(changingColor.ui);
            // console.log('added changing color', changingColor);

            // Parse Specs and create array of inputs 
            var boundingBox = this.ui = $('<div></div>');
            boundingBox.addClass('form');
            var inputs = this.inputs = [];
            // body.append(boundingBox);
            
            var keys = Object.keys(specs);
            keys.forEach((key)=>{
                var prop = new Property(key, specs[key].type);
                inputs.push(prop);
                boundingBox.append(prop.placeholder.ui);

                // console.log('Checking Specs', prop.placeholder, key, specs[key]);
            });

            function Property(key, type){
                var control = { value : null, type : type, key : key };
                this.placeholder = $('<div></div>'); // default value for ui element 

                if(specs[key].type == 'string'){
					this.placeholder = new Input(control);
                    this.placeholder.ui.attr('type', 'text');
                    console.log('String type')
                    console.log('string input', this.placeholder);
				}
				else if(specs[key].type == 'number'){
                    this.placeholder = new Input(control);
                    this.placeholder.ui.attr('type', 'text');
				}
				else if(specs[key].type == 'color'){
                    this.placeholder = new Color(control);

				}else if(specs[key].type == 'boolean'){
                    this.placeholder = new Checkbox(control);

				}else if(specs[key].type == 'array'){
					console.log('Form Generation ERROR: Arrays will be implemented later');

				}else if(specs[key].type == 'properties'){
					this.placeholder = new Input(control);
                    this.placeholder.ui.attr('type', 'text');

				}else if(specs[key].type == 'callback'){
                    console.log('Form Generation ERROR: Callback will be implemented later');
					
				}else if(specs[key].type == 'function'){
					console.log('Form Generation ERROR: Function will be implemented later');
					
				}else if(specs[key].type == 'invalid'){
                    console.log('Form Generation ERROR: Invalid will be implemented later');
					
				}else if(specs[key].type == 'colorscheme'){
					this.placeholder = new Input(control);
                    this.placeholder.ui.attr('type', 'text');
                }else if(specs[key].type == undefined){
					if( specs[key].validItems) {
                        this.placeholder = new ListInput(control, specs[key].validItems);
					}
				}else {
					this.placeholder.text('Unknown property type');
				}

                console.log("Checking placeholder", this.placeholder);
                
                this.getValue = function(){
                    return {   
                        key : key,
                        control : control,
                    }
                }
            }

            function Color(outerControl){
                var label = $('<div></div>');
                label.text(outerControl.key);

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
                surroundingBox.append(label);
                surroundingBox.append(boundingBox);
                
                var RValue = new Slider(control.R);
                var GValue = new Slider(control.G);
                var BValue = new Slider(control.B);
                
                var sliders = $('<div></div>');
                sliders.append(RValue.ui, GValue.ui, BValue.ui);
                
                var color = $('<div></div>');
                // color.width(50);
                // color.height(color.width());
                
                boundingBox.append(sliders);
                boundingBox.append(color);
                
                
                // CSS
                boundingBox.width(400);
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

                // Functionality
                function updatePreview(){
                    var c = `rgb(${control.R.value}, ${control.G.value}, ${control.B.value})`;
                    color.css('background', c);
                    outerControl.value = c;
                    // console.log('rn');
                }

                RValue.update = GValue.update = BValue.update = updatePreview;
                updatePreview();

            }

            function ListInput(control, listElements){
                var label = $('<div></div>');
                label.text(control.key);

                var surroundingBox = this.ui = $('<div></div>');
                var boundingBox = $('<div></div>');
                surroundingBox.append(label);
                surroundingBox.append(boundingBox);

                var select = $('<select></select>');

                boundingBox.append(select);
                
                // Elements of the list 
                var defaultOption = $('<option></option>');
                defaultOption.text('select');

                select.append(defaultOption);

                listElements.forEach((listElement)=>{
                    var option = $('<option></option>');
                    option.text(listElement);
                    option.attr('value', 'value');
                    select.append(option);
                });

                select.on('click', ()=>{
                    control.value = boundingBox.find(":selected").val();
                });
            }
            
            function Input(control){
                var label = $('<div></div>');
                label.text(control.key);    
                
                var surroundingBox = this.ui = $('<div></div>');
                var boundingBox = $('<div></div>');
                surroundingBox.append(label);
                surroundingBox.append(boundingBox);

                var input =  $('<input type="text">');
                boundingBox.append(input);

                input.on('change', ()=>{
                    control.value = input.val();
                });
            }
            
            function Checkbox(control){
                var label = $('<div></div>');
                label.text(control.key); 

                var surroundingBox = this.ui = $('<div></div>');
                var boundingBox = $('<div></div>');
                surroundingBox.append(label);
                surroundingBox.append(boundingBox);

                var checkbox = $('<input type="checkbox" />');
                boundingBox.append(checkbox);

                checkbox.on('click', ()=>{
                    control.value = checkbox.prop('checked');
                });

                // CSS
                label.css('display', 'inline-block');
                boundingBox.css('display', 'inline-block')
            }

            function Slider(control){
                var label = $('<div></div>');
                label.text(control.key); 

                var surroundingBox = this.ui = $('<div></div>');
                surroundingBox.append(label);
                
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
                
                this.update = ()=>{};
                
                slide.on('mousedown', ()=>{
                    setValue = true;
                });

                slide.on('mousemove', ()=>{
                    if(setValue){
                        control.value = slide.val();
                        this.update();
                        console.log('working');
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
                

            }

            function Submit(){
                var boundingBox = this.ui = $('<div></div>');
                var button = $('<button></button>');
                button.text('Submit');

                this.update = ()=>{};

                button.on('click', ()=>{
                    this.update();
                });
            }

        }
        
        return form;
    })();