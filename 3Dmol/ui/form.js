/**
 * This is a colection of contructor to make different input element
 * @function $3Dmol.UI#Form
 */
$3Dmol.UI.Form = (function () {
    /**
     * Create Color input
     * @function $3Dmol.UI#Form.Color
     * @param {Object} outerControl Reference object to store the value
     */
    Form.Color = function (outerControl) {
        const redDot = $('<div></div>');
        redDot.height(10);
        redDot.width(10);
        redDot.css('border-radius', '50%');
        redDot.css('background', 'red');
        redDot.css('margin-right', '3px');

        const blueDot = redDot.clone();
        blueDot.css('background', 'blue');

        const greenDot = redDot.clone();
        greenDot.css('background', 'green');

        const control = this.control = {
            R: {
                value: 0,
                min: 0,
                max: 255,
                label: redDot
            },
            G: {
                value: 0,
                min: 0,
                max: 255,
                label: greenDot
            },
            B: {
                value: 0,
                min: 0,
                max: 255,
                label: blueDot
            },
        };

        const surroundingBox = this.ui = $('<div></div>')
        const boundingBox = $('<div></div>');

        surroundingBox.append(boundingBox);

        const spectrumControl = {
            key: 'Spectrum',
            value: null
        }

        const spectrum = new Form.Checkbox(spectrumControl);

        boundingBox.append(spectrum.ui);

        spectrum.ui.css({
            'margin-left': '2px'
        })


        const RValue = new Form.Slider(control.R);
        const GValue = new Form.Slider(control.G);
        const BValue = new Form.Slider(control.B);

        const sliders = $('<div></div>');
        sliders.append(RValue.ui, GValue.ui, BValue.ui);

        const color = $('<div></div>');

        boundingBox.append(sliders);
        boundingBox.append(color);


        // CSS

        RValue.slide.css('color', 'red');

        // GValue.ui.css('display', 'block');
        GValue.slide.css('color', 'green');

        // BValue.ui.css('display', 'block');
        BValue.slide.css('color', 'blue');

        color.height(15);
        // color.width(50);
        color.css('margin-top', '6px');
        color.css('margin-bottom', '6px');
        color.css('border', '1px solid grey');
        color.css('border-radius', '500px');

        this.update = function (control) {};
        const self = this;
        // Functionality
        function updatePreview(c) {
            c = `rgb(${control.R.value}, ${control.G.value}, ${control.B.value})`;
            color.css('background', c);
            outerControl.value = c;
            self.update(control);
        }

        RValue.update = GValue.update = BValue.update = updatePreview;
        updatePreview();

        spectrum.update = function (v) {
            sliders.toggle();

            if (v.value) {
                color.css({
                    'background': 'linear-gradient(to right, red, orange, yellow, green, blue, indigo, violet)'
                });
                outerControl.value = 'spectrum';
            } else {
                updatePreview();
            }
        }

        this.getValue = function () {
            return outerControl;
        }

        this.validate = function () {
            return true;
        }

        this.setValue = function (colorValue) {

            if (colorValue === 'spectrum') {
                spectrum.setValue(true);
                spectrum.update(spectrumControl);
                sliders.hide();

                outerControl.value = 'spectrum';
            }
        }

        spectrum.ui.hide();

        this.enableSpectrum = function () {
            spectrum.ui.show();
        }
    }

    /**
     * Create ListInput input
     * @function $3Dmol.UI#Form.ListInput
     * @param {Object} control Reference object to store the value
     * @param {Array} listElements list of the elements through which options are generated
     */
    Form.ListInput = function (control, listElements) {
        // var label = $('<div></div>');
        // label.text(control.key);

        const surroundingBox = this.ui = $('<div></div>');
        const boundingBox = $('<div></div>');
        let itemList = listElements;
        // surroundingBox.append(label);
        surroundingBox.append(boundingBox);

        const select = $('<select></select>');
        select.css($3Dmol.defaultCSS.ListInput.select);

        boundingBox.append(select);

        const showAlertBox = this.showAlertBox = true;
        const failMessage = $('<div></div>');
        failMessage.text('Please select some value');
        failMessage.css({
            'color': 'crimson',
            'font-family': 'Arial',
            'font-weight': 'bold',
            'font-size': '10px'
        });
        failMessage.hide();
        boundingBox.append(failMessage);

        this.update = function (control) {}

        select.on('click', {
            parent: this
        }, (event) => {
            control.value = select.children('option:selected').val();
            event.data.parent.update(control);
        });

        this.getValue = () => control

        // this.preventAlertBox = function(){
        //     show
        // }

        this.validate = function () {
            if (control.value === 'select' || control.value == null) {
                // eslint-disable-next-line no-unused-expressions
                (this.showAlertBox) ? failMessage.show(): null;
                select.css({
                    'box-shadow': '0px 0px 2px red'
                });
                return false;
            } 
                failMessage.hide();
                boundingBox.css({
                    'box-shadow': 'none'
                });
                return true;
            
        }

        this.setValue = function (val) {
            if (listElements.indexOf(val) !== -1) {
                select.empty();
                const defaultOption = $('<option></option>');
                defaultOption.text('select');

                itemList.forEach((item) => {
                    const option = $('<option></option>');
                    option.text(item);
                    option.attr('value', item);
                    select.append(option);

                    if (val === item) {
                        option.prop('selected', true);
                    }
                });

                control.value = select.children('option:selected').val();
            } else {
                console.error('UI::Form::ListInput:incorrect value', val);
            }
        }

        this.updateList = function (newList) {
            select.empty();

            const defaultOption = $('<option></option>');
            defaultOption.text('select');
            defaultOption.attr('value', 'select');

            select.append(defaultOption);

            itemList = newList;
            itemList.forEach((item) => {
                const option = $('<option></option>');
                option.text(item);
                option.attr('value', item);
                select.append(option);
            });
        }

        this.updateList(itemList);
    }

    /**
     * Create text, numeric or range Input
     * @function $3Dmol.UI#Form.Input
     * @param {Object} control Reference object to store the value
     */
    Form.Input = function (control) {
        const surroundingBox = this.ui = $('<div></div>');
        const boundingBox = $('<div></div>');
        // surroundingBox.append(label);
        surroundingBox.append(boundingBox);

        let validationType = this.validationType = 'text';

        surroundingBox.css({
            'width': '100%',
            'box-sizing': 'border-box'
        })

        const input = this.domElement = $('<input type="text">');
        boundingBox.append(input);

        const alertBox = $('<div></div>');
        alertBox.css({
            'border': '1px solid darkred',
            'border-radius': '3px',
            'font-family': 'Arial',
            'font-size': '10px',
            'font-weight': 'bold',
            'margin': '2px',
            'margin-left': '4px',
            'padding': '2px',
            'color': 'darkred',
            'background': 'lightcoral'
        });

        const alertMessage = {
            'invalid-input': 'Invalid input please check the value entered',
        }

        boundingBox.append(alertBox);
        alertBox.hide();

        this.setWidth = function (width) {
            input.width(width - 6);
        }

        this.setWidth(75);

        input.css({
            // 'margin-left': '4px'
        });

        this.update = function (control) {

        }

        input.on('change', {
            parent: this,
            control
        }, (event) => {
            inputString = input.val();

            if (inputString[inputString.length - 1] === ',') {
                inputString = inputString.slice(0, -1);
            }

            if (validationType === 'range') {
                control.value = inputString.split(',');
            } else {
                control.value = inputString;
            }

            // calling update function 
            event.data.parent.update(control);
        });

        let selectedText = null;

        input.on('select', (e) => {
            selectedText = input.val().substring(e.target.selectionStart, e.target.selectionEnd);
        });


        this.getValue = () => control

        const error = this.error = function (msg) {
            alertBox.show();
            alertBox.text(msg)
        }

        this.setValue = function (val) {

            if (validationType === 'range') {
                const text = val.join(',');
                input.val(text);
            } else {
                input.val(val);
            }

            control.value = val;
        }



        function checkInputFloat() {
            const inputString = input.val();

            const dots = inputString.match(/\./g) || [];
            const checkString = inputString.replaceAll(/\./g, '').replaceAll(/[0-9]/g, '');

            if (dots.length > 1) {
                return false
            };

            if (checkString !== '') return false;

            if (isNaN(parseFloat(inputString))) {
                return false;
            } 
                return true;
            
        }

        function checkInputNumber() {
            const inputString = input.val();

            const checkString = inputString.replaceAll(/[0-9]/g, '');

            if (checkString !== '') return false;

            if (isNaN(parseInt(inputString))) {
                return false;
            } 
                return true;
            
        }

        // Parse Input Range Functions

        // Checks only number, comma and hyphen present
        function checkRangeTokens(inputString) {
            const finalString = inputString.replaceAll(',', '').replaceAll('-', '').replaceAll(/[0-9]/g, '').replaceAll(' ', '');;

            if (finalString === '')
                return true;
            return false;
        }

        function checkList(inputString, submit = false) {
            inputString = inputString.replaceAll(' ', '');

            if (inputString[inputString.length - 1] === ',') {
                inputString = inputString.slice(0, -1);
            }


            const rangeList = inputString.split(',');

            // If dublicate comma return false;
            if (/,,/g.exec(inputString)) return false;

            // If first element not a number return false;
            if (isNaN(parseInt(rangeList[0]))) return false;

            const validRangeList = rangeList.map((rangeInput) => checkRangeInput(rangeInput));

            return validRangeList.find((e) => e === false) === undefined;
        }

        function checkRangeInput(inputString, submit = false) {
            const rangeInputs = inputString.split('-');
            if (rangeInputs.length > 2) {
                return false;
            } 
                if (rangeInputs.length === 0) {
                    return true;
                } if (rangeInputs.length === 1) {
                    if (isNaN(parseInt(rangeInputs[0])))
                        return false;
                    return true;
                } if (rangeInputs.length === 2) {
                    if (isNaN(parseInt(rangeInputs[0])) || isNaN(parseInt(rangeInputs[1])))
                        return false;
                    return true;
                } return false;
            
        }

        const checkInput = this.checkInput = function () {
            const inputString = input.val();

            if (validationType === 'number') {
                if (checkInputNumber()) {
                    alertBox.hide();
                    return true;
                } 
                    error(alertMessage['invalid-input']);
                    return false;
                
            } if (validationType === 'float') {
                if (checkInputFloat()) {
                    alertBox.hide();
                    return true;
                } 
                    error(alertMessage['invalid-input']);
                    return false;
                
            } if (validationType === 'range') {
                if (checkRangeTokens(inputString)) {
                    if (checkList(inputString)) {
                        alertBox.hide();
                        return true;
                    } 
                        error(alertMessage['invalid-input']);
                        return false;
                    
                } 
                    error(alertMessage['invalid-input']);
                    return false;
                

            } 
                return true;
            
        }

        this.validateOnlyNumber = function (floatType = false) {
            if (floatType) {
                validationType = 'float';
            } else {
                validationType = 'number';
            }

            input.on('keydown keyup paste cut', (event) => {
                checkInput();
            });
        }


        this.validateInputRange = function () {
            validationType = 'range';

            input.on('keydown keyup paste cut', () => {
                checkInput();
            });

        }

        this.isEmpty = function () {
            if (control.value === "") {
                return true;
            }
        }

        this.validate = function () {
            if ((control.active === true && control.value != null && control.value !== "" && checkInput()) || (control.active === false)) {
                input.css('box-shadow', 'none');
                return true
            } 
                input.css('box-shadow', '0px 0px 2px red');
                return false;
            
        }

        // CSS 

        input.css($3Dmol.defaultCSS.Input.input);
        boundingBox.css($3Dmol.defaultCSS.Input.boundingBox);

    }

    /**
     * Create Checkbox input for boolean values
     * @function $3Dmol.UI#Form.Checkbox
     * @param {Object} control Reference object to store the value
     */
    Form.Checkbox = function (control) {
        const label = $('<div></div>');
        label.text(control.key);
        label.css($3Dmol.defaultCSS.TextDefault);

        const surroundingBox = this.ui = $('<div></div>');
        const boundingBox = $('<div></div>');
        surroundingBox.append(boundingBox);
        surroundingBox.append(label);

        const checkbox = $('<input type="checkbox" />');
        boundingBox.append(checkbox);

        this.click = () => {};

        this.update = function (control) {

        }

        this.getValue = () => control

        checkbox.on('click', {
            parent: this
        }, (event) => {
            control.value = checkbox.prop('checked');
            event.data.parent.update(control);
        });

        // CSS
        label.css('display', 'inline-block');
        boundingBox.css('display', 'inline-block')

        this.validate = function () {
            return true;
        }

        this.setValue = function (val) {
            checkbox.prop('checked', val);
            this.update(control);
            control.value = val;
        }
    }

    /**
     * Create input for values between two numbers
     * @function $3Dmol.UI#Form.Slider
     * @param {Object} control Reference object to store the value
     */
    Form.Slider = function (control, style = {}) {
        const surroundingBox = this.ui = $('<div></div>');

        const boundingBox = $('<div></div>');
        surroundingBox.append(boundingBox);

        boundingBox.css('display', 'flex');
        const slide = this.slide = $('<input type="range">');
        slide.css('width', '100%');

        const min = control.min || 0;
        const max = control.max || 100;
        const step = control.step || 1;
        const defaultValue = control.default || min;
        const labelContent = control.label || '';

        const label = $('<div></div>');
        label.append(labelContent);
        boundingBox.append(label);

        slide.attr('min', min);
        slide.attr('max', max);
        slide.attr('step', step);
        slide.attr('value', defaultValue);
        control.value = defaultValue;
        boundingBox.append(slide);

        let setValue = false;

        this.update = function (control) {

        };

        this.getValue = () => control

        slide.on('mousedown', () => {
            setValue = true;
        });

        slide.on('mousemove', {
            parent: this
        }, (event) => {
            if (setValue) {
                control.value = slide.val();
                event.data.parent.update(control);
            }
        });

        slide.on('mouseup', () => {
            setValue = false;
        });

        // CSS
        boundingBox.css('align-items', 'center');
        boundingBox.height('21px');
        // boundingBox.css('border-radius', '2px');
        // label.css('line-height', '21px');
        slide.css('padding', '0px');
        slide.css('margin', '0px');

        this.validate = function () {
            return true;
        }

        this.setValue = function (val) {
            slider.val(val);
            control.value = slider.val();
        }


    }

    /**
     * Create empty element used for property that whose input cannot be taken
     * @function $3Dmol.UI#Form.EmptyElement
     * @param {Object} control Reference object to store the value
     */
    Form.EmptyElement = function (control) {
        this.ui = $('<div></div>');

        let update = () => {};

        this.onUpdate = (callback) => {
            update = callback;
        }

        this.getValue = () => control

        this.validate = function () {
            return true;
        }
    }

    // mainControl param will be used to take in specName
    // in the form of key 
    // type will be 'form'
    // active will be used to activate deactivate form if more than one form
    /**
     * Creates Form input that takes input from different input element 
     * 
     * @function $3Dmol.UI#Form
     * @param {validSelectionSpec|validStyleSpec|validAtomSpec} specs the defination of spec is used as an input to generate the form
     * @param {Object} mainControl Reference of variable to store the value from the form
     */
    function Form(specs, mainControl) {
        specs = specs || {};
        const boundingBox = this.ui = $('<div></div>');

        const heading = $('<div></div>');
        heading.text(mainControl.key);

        // Styling heading 
        heading.css({
            'border-bottom': '1px solid black',
            'font-family': 'Arial',
            'font-size': '14px',
            'font-weight': 'bold',
            'padding-top': '2px',
            'padding-bottom': '4px'
        });

        boundingBox.append(heading);
        boundingBox.addClass('form');

        const inputs = this.inputs = [];
        // body.append(boundingBox);

        const keys = Object.keys(specs);
        keys.forEach((key) => {
            if (specs[key].gui) {
                const prop = new Property(key, specs[key].type);
                inputs.push(prop);
                boundingBox.append(prop.ui);
            }

        });
        const self = this;

        this.update = function () {}

        const update = (control) => {

        };


        inputs.forEach((input) => {
            input.update = update;
        })

        this.getValue = function () {
            mainControl.value = {};

            inputs.forEach((input) => {
                const inputValue = input.getValue();

                if (inputValue.active) {
                    mainControl.value[inputValue.key] = inputValue.value;
                }
            });

            return mainControl;
        }

        const updateValues = function (inputControl) {
            mainControl.value[inputControl.key] = control.value;
            update(mainControl);
        }

        this.validate = function () {
            const validations = inputs.map((i) => {

                if (i.active.getValue().value) {
                    return i.placeholder.validate();
                } 
                    return true;
                
            });


            if (validations.find(e => e === false) === undefined)
                return true;
            
                return false;
            

        }

        this.setValue = function (val) {
            const keys = Object.keys(val);
            for (let i = 0; i < keys.length; i++) {
                const input = inputs.find((e) => {
                    if (e.control.key === keys[i])
                        return e;
                });



                input.placeholder.setValue(val[keys[i]]);
                input.active.setValue(true);
                input.placeholder.ui.show();
                input.control.active = true;
            }

            // mainControl.value = val;
            this.update(mainControl);
            const v = this.getValue();

        }

        this.getInputs = function () {
            return inputs;
        }

        function Property(key, type) {
            const control = this.control = {
                value: null,
                type,
                key,
                active: false
            };
            const boundingBox = this.ui = $('<div></div>');
            this.placeholder = {
                ui: $('<div></div>')
            }; // default value for ui element 
            this.active = new Form.Checkbox({
                value: false,
                key
            });


            if (specs[key].type === 'string' || specs[key].type === 'element') {
                this.placeholder = new Form.Input(control);
                this.placeholder.ui.attr('type', 'text');
            } else if (specs[key].type === 'number') {

                let slider = false;

                if (specs[key].min !== undefined && specs[key].max !== undefined && specs[key].default !== undefined) {
                    slider = true;
                }

                if (slider) {
                    // if( specs[key].min && spec[key].max){
                    control.min = specs[key].min;
                    control.max = specs[key].max;
                    control.default = specs[key].default;
                    control.step = specs[key].step || ((control.max - control.max) / 1000);
                    this.placeholder = new Form.Slider(control);
                } else {
                    this.placeholder = new Form.Input(control);
                    this.placeholder.ui.attr('type', 'text');
                    this.placeholder.validateOnlyNumber(specs[key].floatType);
                }
            } else if (specs[key].type === 'array_range') {
                this.placeholder = new Form.Input(control);
                this.placeholder.ui.attr('type', 'text');
                this.placeholder.validateInputRange();
            } else if (specs[key].type === 'color') {
                this.placeholder = new Form.Color(control);
                if (specs[key].spectrum) {
                    this.placeholder.enableSpectrum();
                }

            } else if (specs[key].type === 'boolean') {
                this.placeholder = new Form.Checkbox(control);

            } else if (specs[key].type === 'properties') {
                this.placeholder = new Form.Input(control);
                this.placeholder.ui.attr('type', 'text');

            } else if (specs[key].type === 'colorscheme') {
                this.placeholder = new Form.ListInput(control, Object.keys($3Dmol.builtinColorSchemes));
                this.placeholder.ui.attr('type', 'text');

            } else if (specs[key].type === undefined) {
                if (specs[key].validItems) {
                    this.placeholder = new Form.ListInput(control, specs[key].validItems);
                }

            } else if (specs[key].type === 'form') {
                this.placeholder = new Form(specs[key].validItems, control);
                this.placeholder.ui.append($('<div></div>').css($3Dmol.defaultCSS.LinkBreak));
            } else {
                this.placeholder = new Form.EmptyElement(control);
                // return new Form.EmptyElement(control);
            };

            this.getValue = function () {

                if (this.placeholder.getValue)
                    return this.placeholder.getValue();
                return null;
            }


            // Adding active control for the property
            const {placeholder} = this;

            if (type !== 'boolean') {
                placeholder.ui.hide();
                boundingBox.append(this.active.ui);
                this.active.update = function (c) {
                    (c.value) ? placeholder.ui.show(): placeholder.ui.hide();
                    control.active = c.value;
                }
            } else {
                this.placeholder.update = function (c) {
                    control.active = c.value;
                }
            }

            boundingBox.append(this.placeholder.ui);

            if (this.placeholder.onUpdate)
                this.placeholder.onUpdate(updateValues);
        }


    }



    return Form;
})();