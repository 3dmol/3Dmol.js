    $3Dmol.UI.form = (function(){
        function form(specs){
            specs = specs || {};

            // Test area
            body = $('body');
            var changingColor = new Color();
            body.append(changingColor.ui);
            console.log('added changing color', changingColor);

            // Parse Specs and create array of inputs 
            var inputs = [];
            
            var keys = Object.keys(specs);
            keys.forEach((key)=>{
                
                // add input element object to the list

            });

            function Property(key, type){
                config = config || {};
                
                
            }

            function Color(){
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


                var control = {
                    R : { value : 0, min : 0, max : 255, label : redDot },
                    G : { value : 0, min : 0, max : 255, label : blueDot },
                    B : { value : 0, min : 0, max : 255, label : greenDot },
                }

                var boundingBox = this.ui = $('<div></div>');
                var RValue = new slider(control.R);
                var GValue = new slider(control.G);
                var BValue = new slider(control.B);

                var sliders = $('<div></div>');
                sliders.append(RValue.ui, GValue.ui, BValue.ui);

                var color = $('<div></div>');
                // color.width(50);
                // color.height(color.width());

                boundingBox.append(sliders);
                boundingBox.append(color);


                // CSS
                boundingBox.width(400);
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
                    console.log('rn');
                }

                RValue.update = GValue.update = BValue.update = updatePreview;
                updatePreview();
            }

            function list(control, listElements){
                var boundingBox = this.ui = $('<select></select>');
                
                // Elements of the list 
                var defaultOption = $('<option></option>');
                defaultOption.text('select');

                boundingBox.append(defaultOption);

                listElements.forEach((listElement)=>{
                    var option = $('<option></option>');
                    option.text(listElement);
                    option.attr('value', 'value');
                    boundingbox.append(option);
                });

                boundingBox.on('click', ()=>{
                    control.value = boundingBox.find(":selected").val();
                });
            }
            
            function input(control){
                var boundingBox = this.ui = $('<input type="text">');

                boundingBox.on('change', ()=>{
                    control.value = boundingBox.val();
                });
            }
            
            function checkbox(control){
                var boundingBox = this.ui = $('<input type="checkbox" />');

                boundingBox.on('click', ()=>{
                    control.value = checkbox.prop('checked');

                });
            }

            function slider(control){
                var boundingBox = this.ui = $('<div></div>');
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

            function submit(control){

            }

        }
        
        return form;
    })();