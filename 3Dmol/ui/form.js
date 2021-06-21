    $3Dmol.UI.form = (function(){
        function form(specs){
            specs = specs || {};

            // Parse Specs and create array of inputs 
            var inputs = [];
            
            var keys = Object.keys(specs);
            keys.forEach((key)=>{
                
                // add input element object to the list

            });

            function property(key, type){
                config = config || {};
                
                
            }

            function color(){
                var R = 0;
                var G = 0;
                var B = 0;

                var boundingBox = $('<div></div>');
                var RValue = new 
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
                var checkbox = this.ui = $('<input type="checkbox" />');

                checkbox.on('click', ()=>{
                    control.value = checkbox.prop('checked');

                });
            }

            function slider(config){
                config = config || {}
                var step = config.step || 1;
                var min = config.min || 0;
                var max = config.max || 100;
                var current = config.default || min;

                var boundingBox = this.ui = $('<div></div>');
                var sliderBox = $('<div></div>');
                var slide = $('<div></div>');
                var marker = $('<div></div>');
        
                var height = config.height = config.height || 12;
                var width = config.width = config.width || 100;
                sliderBox.css('box-sizing', 'border-box');
        
                marker.css('height', height + 'px');
                marker.css('width', height + 'px');
                marker.css('background', 'rgb(68, 90, 235)');
                marker.css('border-radius', height/2 + 'px');
        
                slide.css('height', '3px');
                slide.css('width', width + 'px');
                slide.css('background', 'grey');
                slide.css('border-radius', '1.5px');
        
                sliderBox.append(slide);
                sliderBox.append(marker);
        
                sliderBox.css('height', height + 'px');
                sliderBox.css('position', 'relative');
        
                slide.css('position', 'absolute');
                slide.css('top', height/2 - 1.5 + 'px');
                marker.css('position', 'absolute');
        
                boundingBox.append(sliderBox);
                boundingBox.css('width', width);
                // boundingBox.css('position','absolute');
        
                //vertically centering the slide
        
                this.demo = sliderBox;

                this.updateSlider = function(control){
                    control.value = current;
                }

                // Adding Mouse move 
            }

            function submit(control){

            }

        }
        
        return form;
    })();