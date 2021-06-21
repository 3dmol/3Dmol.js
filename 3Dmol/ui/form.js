    $3Dmol.UI.form = (function(){
        function form(specs){
            specs = specs || {};

            // Test area
            body = $('body');
            var changingColor = new Color();



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
                var control = {
                    R : { value : 0, min : 0, max : 255 },
                    G : { value : 0, min : 0, max : 255 },
                    B : { value : 0, min : 0, max : 255 },
                }

                var boundingBox = $('<div></div>');
                var RValue = new slider(control.R);
                var GValue = new slider(control.G);
                var BValue = new slider(control.B);

                var sliders = $('<div></div>');
                slider.append(RValue, GValue, BValue);

                var color = $('<div></div>');
                color.width(10);
                color.height(color.width());
                boundingBox.append(slider);
                boundingBox.append(color);
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
                var boundingBox = this.ui = $('<input type="range">');
                var setValue = false;
                
                this.update = ()=>{};
                
                boundingBox.on('mousedown', ()=>{
                    setValue = true;
                });

                boundingBox.on('mousemove', ()=>{
                    if(setValue){
                        control.value = boundingBox.val();
                        this.update();
                    }
                });

                boundingBox.on('mouseup', ()=>{
                    setValue = false;
                });

            }

            function submit(control){

            }

        }
        
        return form;
    })();