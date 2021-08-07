  /*
  * This is ui that we are trying to draw at the time of generationo
  * of viewer.
  * Refer to line number 43 of the glViewer to see the viewer code
  */

  /**
   * $3Dmol.UI - UI creates panels of viewer to assist control of the viewport in
   * 3dmol.js
   *
   * @param  {$3Dmol.StateManager} Manages states of the UI
   * @return {object} configuration for the UI
   */
  $3Dmol.UI = (function(){
    function UI(stateManager, config, parentElement){
      config = config || {}

      // Extract the viewer and then render it
      var icons =new $3Dmol.UI.Icons();
      var body = $('body');
      
      var mainParent = $(parentElement[0]);
      console.log("Main Parent", mainParent.css('position'), config);
      // Generates the necessary UI elements
      var HEIGHT = config.height;
      var WIDTH = config.width;
      var uiElements = this.tools = generateUI(config);
      
      /**
       * @function generateUI creates all the jquery object of different UI features
       * @param  {object} config
       */
      function generateUI(config){    
        var movieControl = new MovieBar();
        mainParent.append(movieControl.ui);
        // movieControl.ui.css('top', 200);
        // movieControl.ui.css('z-index', 99);
        setLocation(mainParent, movieControl.ui, 'center', 'bottom');
        movieControl.ui.hide();
        
        var topbar = new Toolbar();
        mainParent.append(topbar.ui);
        setLocation(mainParent, topbar.ui, 'left', 'top');
        topbar.ui.hide();

        var selectionBox = new SelectionBox(icons.select);
        mainParent.append(selectionBox.ui);
        setLocation(mainParent, selectionBox.ui, 'left', 'top');

        var dialog = new DialogBox({ height: 300, width: 300 });
        mainParent.append(dialog.ui);
        setLocation(mainParent, dialog.ui, 'center', 'center');
        dialog.ui.hide();
        
        var alertBox = new AlertBox({ width : 100 });
        // console.log(alertBox);
        mainParent.append(alertBox.ui);
        setLocation(mainParent, alertBox.ui, 'right', 'top');

        var contextMenu = new ContextMenu();
        mainParent.append(contextMenu.ui);
        // setLocation(mainParent, contextMenu.ui, 'center', 'center');
        setPosition(contextMenu.ui, 100, 100)
        
        var surfaceMenu = new SurfaceMenu();
        mainParent.append(surfaceMenu.ui);
        setLocation(mainParent, surfaceMenu.ui, 'right', 'top', );

        return {
          topbar : topbar,
          selectionBox : selectionBox,
          dialog : dialog,
          alertBox : alertBox,
          movieControl : movieControl,
          contextMenu : contextMenu,
          surfaceMenu : surfaceMenu
        } 
      }
      /**
       * @function UI_Overlay adds overlay on the glviewer to assists placement of different UI elements on the viewer.
       * @param  {object} config
       */

      /**
      * @function Toobar creates horizontal toolbar for the UI
      * @param {Object} config : Stores the top, left, width, padding property for the toolbar
      * @return {Object} Jquery div object
      **/
      function Toolbar(config){
        config = config || {};

        var top = config.top || 10;
        var left = config.left || 10;
        var padding = config.padding || '3px';

        var toolbar = this.ui = $('<div></div>');

        // Toolbar content
        var optionButton = new button(icons.option, 20)
        var moveButton = new button(icons.move, 20);
        var rotateButton = new button(icons.rotate, 20);
        toolbar.append(optionButton.ui);
        toolbar.append(moveButton.ui);
        toolbar.append(rotateButton.ui);

        // Toolbar design
        toolbar.css('border', 'none');
        toolbar.css('padding', padding);
        toolbar.css('background', 'none');
        toolbar.css('color', 'white');
        toolbar.css('box-sizing', 'border-box');
        toolbar.css('position', 'absolute');

        rotateButton.ui.css('margin-right', '0px');

        // Button Action
        optionButton.ui.click(()=>{
          console.log('option button clicked');
        });

        moveButton.ui.click(()=>{
          console.log('move button clicked');
        });

        rotateButton.ui.click(()=>{
          console.log('click button clicked');
        });
      }

      /**
       * @function SelectionBox - Draws the box where all the selections on atoms are listed
       * This will be used to modify style for specific set of selection
       * @param {SVG} $3Dmol.UI.Icons stores all the SVG elements needed of the UI 
       * @return {Object}  Jquery element of div
       */
      function SelectionBox(icon, side='left') {
        var selectionBox = this.ui = $('<div></div>');
        _editingForm = false;

        var selections = $('<div></div>');
        var scrollBox = $('<div></div>');

        selections.css('opacity', '0.9');
        
        var showArea = $('<div></div>');
        var addArea = $('<div></div>');
        var plusButton = new button(icons.plus, 20);
        plusButton.ui.css('margin','0px');
        
        var hideButton = new button(icon, 20);
        this.selectionObjects = [];

        // Content
        selectionBox.append(hideButton.ui);
        selectionBox.append(showArea);
        selectionBox.css('position', 'absolute');
        
        scrollBox.append(selections);
        showArea.append(scrollBox);
        addArea.append(plusButton.ui);
        showArea.append(addArea);

        // CSS
        if(side == 'left'){
          selectionBox.css('text-align', 'left');
        }
        else if(side == 'right') {
          selectionBox.css('text-align', 'right');
        }
        else {
          // Add alert box code
          selectionBox.css('text-align', 'right');
        }

        showArea.css('box-sizing', 'border-box');
        showArea.css('padding', '3px');
        // showArea.css('width', '162px');

        // scrollBox.css('max-height', HEIGHT*0.8);
        // scrollBox.css('overflow', 'hidden');
        
        // selections.css('max-height', HEIGHT*0.8);
        // selections.css('overflow', 'auto');
        // selections.css('overflow-x', 'visible');
        selections.css('box-sizing', 'content-box');

        // Action
        var hidden = true;
        showArea.hide();

        hideButton.ui.click(toggleHide);

        function toggleHide(){
          if(hidden){
            showArea.show(100);
          }
          else {
            showArea.hide(100);
          }
          hidden = !hidden;
        }

        function Selection(){
          var boundingBox = this.ui = $('<div></div>');
          var sid = this.id = null;
          boundingBox.css({
            'background' : '#e8e8e8',
            'padding' : '4px 4px 2px 4px',
            'border-radius' : '6px',
            'margin-bottom' : '3px',
            'position':'relative'
          });

          var header = $('<div></div>');
          boundingBox.append(header);
          var heading = $('<div></div>');
          var controls = $('<div></div>');

          header.append(heading, controls);
          heading.css({
            'font-family' : 'Arial',
            'font-weight': 'bold',
            'font-size':'12px',
            'display':'inline-block',
            'width' : '60px'
          });

          controls.css({
            'display' : 'inline-block'
          });

          header.hide();
          controls.editMode = false;

          var removeButton = new button(icons.minus, 16, { bfr:0.5, backgroundColor:'#f06f6f'});
          var editButton = new button(icons.pencil, 16);
          var visibleButton = new button(icons.visible, 16);

          controls.append(removeButton.ui)
          controls.append(editButton.ui);
          controls.append(visibleButton.ui);

          var parameters = $('<div></div>');
          boundingBox.append(parameters);

          removeButton.ui.on('click', function(){
            stateManager.removeSelection(sid);
            boundingBox.detach();
            delete this;
          });

          editButton.ui.on('click', function(){
            parameters.toggle();
          });

          var hidden = false;
          visibleButton.ui.on('click', ()=>{
            stateManager.toggleHide(sid);
            if(hidden){
              hidden = false;
              visibleButton.setSVG(icons.visible);
            }
            else {
              hidden = true;
              visibleButton.setSVG(icons.invisible);
            }
          });

          var styleBox = new StyleBox();

          var showStyle = false;
          boundingBox.append(styleBox.ui);
          styleBox.ui.css({
            'position' : 'absolute',
            'left' : '100%',
            'top' : 0,
            'border-radius' : '4px'
          });

          styleBox.ui.hide();

          var allControl = this.allSelector = {
            key : 'Select All Atom',
            value : null,
            active : true
          }

          var allCheckBox = new $3Dmol.UI.Form.Checkbox(allControl);
          parameters.append(allCheckBox.ui);


          var selectionFormControl = this.selectionValue = {
            key : 'Selection Spec',
            value : null,
            active : true
          }
          
          var selectionSpecForm = new $3Dmol.UI.Form($3Dmol.GLModel.validAtomSelectionSpecs, selectionFormControl);
          parameters.append(selectionSpecForm.ui);

          var submitControls = $('<div></div>');
          var submit = new button(icons.tick, 16, { backgroundColor : 'lightgreen'});
          var cancel = new button(icons.cross, 16, { backgroundColor : 'lightcoral'});
          submitControls.append(submit.ui, cancel.ui);

          
          var alertBox = $('<div></div>');
          parameters.append(alertBox);
          alertBox.css({
            'color': 'darkred',
            'border':'1px solid darkred',
            'border-radius': '3px',
            'background-color':'lightcoral',
            'padding':'3px',
            'text-align':'center',
            'font-family':'Arial',
            'font-size':'12px',
            'font-weight':'bold'
          });
          alertBox.hide();
          
          parameters.append(submitControls);

          
          allCheckBox.update = function(){
            selectionSpecForm.ui.toggle();
          }

          function finalizeSelection(id){
            header.show();
            controls.editMode = true;
            sid = id;
            heading.text('Sel#' + id);
            parameters.hide();
            showStyle = true;
            styleBox.setSid(id);
          }

          function checkAndAddSelection(sid = null){
            var validate = selectionSpecForm.validate();
            if(validate){
              selectionSpecForm.getValue();
              console.log('Selection Form Value', selectionFormControl, sid);
              var checkAtoms = stateManager.checkAtoms(selectionFormControl.value);
              
              if(checkAtoms){
                var id = stateManager.addSelection(selectionFormControl.value, sid);
                finalizeSelection(id);
                if(sid == null ) _editingForm = false;
              }
              else {
                alertBox.text('No atom selected');
                alertBox.show();
                setTimeout(()=>{
                  alertBox.hide();
                }, 5000);
              }
            }
            else {
              alertBox.text('Invalid Input');
              alertBox.show();
              setTimeout(()=>{
                alertBox.hide();
              }, 5000);
            }
          }

          submit.ui.on('click', ()=>{
            if(controls.editMode == false){
              if(allControl.value){
                var id = stateManager.addSelection({});
                finalizeSelection(id);
                _editingForm = false;
              }
              else{
                checkAndAddSelection(); 
              }

            }
            else {
              if(allControl.value){
                var id = sid
                stateManager.addSelection({}, id);
                finalizeSelection(id);
              }
              else{
                var id = sid;
                checkAndAddSelection(id);
              }
            }
          });

          cancel.ui.on('click', ()=>{
            if(controls.editMode){
              parameters.hide();
            }
            else {
              boundingBox.detach();
              delete this;
              _editingForm = false;
            }
          });

          boundingBox.on('mouseenter', ()=>{
            if(showStyle){
              styleBox.ui.show();
              console.log("entering", showStyle)
            }
          });

          boundingBox.on('mouseleave', ()=>{
            if(showStyle){
              styleBox.ui.hide();
              console.log("exiting", showStyle)
            }
          });

          boundingBox.on('keyup', (e)=>{
            if(e.key == 'Enter'){
              submit.ui.trigger('click');
            }
          });


        }

        plusButton.ui.on('click', ()=>{
          if(!_editingForm){
            var newSelection = new Selection();
            selections.append(newSelection.ui);
            _editingForm = true;
          }
          
        });   
      }

      
      /**
       * @function StyleBox creates styleBox for listing out different styles for a particular selection
       * @return {Object} Jquery dom object
       */
       function StyleBox(selId, side='left') {
        var styleBox = this.ui = $('<div></div>');
        _editingForm = false;
        var sid = this.sid = selId; // selection id

        this.setSid = function(id){
          sid = this.sid = id;
        }

        var styles = $('<div></div>');
        var scrollBox = $('<div></div>');

        styles.css('opacity', '0.9');
        
        var showArea = $('<div></div>');
        var addArea = $('<div></div>');
        addArea.css('text-align' , 'center');

        var plusButton = new button(icons.plus, 20);
        plusButton.ui.css('margin','0px');
      
        this.selectionObjects = [];

        // Content
        styleBox.append(showArea);
        styleBox.css('position', 'absolute');
        
        scrollBox.append(styles);
        showArea.append(scrollBox);
        addArea.append(plusButton.ui);
        showArea.append(addArea);

        // CSS
        if(side == 'left'){
          styleBox.css('text-align', 'left');
        }
        else if(side == 'right') {
          styleBox.css('text-align', 'right');
        }
        else {
          // Add alert box code
          styleBox.css('text-align', 'right');
        }

        showArea.css('box-sizing', 'border-box');
        showArea.css('padding', '3px');
        showArea.css('width', '162px');
        showArea.css('background-color', '#a4a4a4')
        showArea.css('border-radius', '4px');

        scrollBox.css('max-height', HEIGHT*0.8);
        scrollBox.css('overflow', 'hidden');
        
        styles.css('max-height', HEIGHT*0.8);
        styles.css('overflow', 'auto');
        styles.css('box-sizing', 'content-box');


        function Style(sid){
          var boundingBox = this.ui = $('<div></div>');
          var stid = this.id = null; // style id 
          boundingBox.css({
            'background' : '#e8e8e8',
            'padding' : '4px 4px 2px 4px',
            'border-radius' : '6px',
            'margin-bottom' : '3px',
            'position':'relative'
          });

          var header = $('<div></div>');
          boundingBox.append(header);
          var heading = $('<div></div>');
          var controls = $('<div></div>');

          header.append(heading, controls);
          heading.css({
            'font-family' : 'Arial',
            'font-weight': 'bold',
            'font-size':'12px',
            'display':'inline-block',
            'width' : '60px'
          });

          controls.css({
            'display' : 'inline-block'
          });

          header.hide();
          controls.editMode = false;

          var removeButton = new button(icons.minus, 16, { bfr:0.5, backgroundColor:'#f06f6f'});
          var editButton = new button(icons.pencil, 16);
          var visibleButton = new button(icons.visible, 16);

          controls.append(removeButton.ui)
          controls.append(editButton.ui);
          controls.append(visibleButton.ui);

          var parameters = $('<div></div>');
          boundingBox.append(parameters);

          removeButton.ui.on('click', { parent: this, stid : stid }, function(e){
            console.log(e.data.parent, e.data.stid, stid)
            stateManager.removeStyle(sid, stid);
            boundingBox.detach();
            delete this;
          });

          editButton.ui.on('click', function(){
            parameters.toggle();
          });

          var hidden = false;
          visibleButton.ui.on('click', ()=>{
            stateManager.toggleHideStyle(sid, stid);
            if(hidden){
              hidden = false;
              visibleButton.setSVG(icons.visible);
            }
            else {
              hidden = true;
              visibleButton.setSVG(icons.invisible);
            }
          });

          var styleFormControl = this.selectionValue = {
            key : 'Selection Spec',
            value : null,
            active : true
          }
          
          var styleSpecForm = new $3Dmol.UI.Form($3Dmol.GLModel.validAtomStyleSpecs, styleFormControl);
          parameters.append(styleSpecForm.ui);

          var submitControls = $('<div></div>');
          var submit = new button(icons.tick, 16, { backgroundColor : 'lightgreen'});
          var cancel = new button(icons.cross, 16, { backgroundColor : 'lightcoral'});
          submitControls.append(submit.ui, cancel.ui);

          
          var alertBox = $('<div></div>');
          parameters.append(alertBox);
          alertBox.css({
            'color': 'darkred',
            'border':'1px solid darkred',
            'border-radius': '3px',
            'background-color':'lightcoral',
            'padding':'3px',
            'text-align':'center',
            'font-family':'Arial',
            'font-size':'12px',
            'font-weight':'bold'
          });
          alertBox.hide();
          
          parameters.append(submitControls);

          function finalizeStyle(id){
            header.show();
            controls.editMode = true;
            stid = id;
            heading.text('Sty#' + id);
            parameters.hide();
          }

          function checkAndAddStyle(stid = null){
            var validate = styleSpecForm.validate();
            if(validate){
              styleSpecForm.getValue();
              console.log('Style Form Value', styleFormControl, sid, stid);
              var id = stateManager.addStyle(styleFormControl.value, sid, stid);
              finalizeStyle(id);
              if(stid == null) _editingForm = false;
            }
            else {
              alertBox.text('Invalid Input');
              alertBox.show();
              setTimeout(()=>{
                alertBox.hide();
              }, 5000);
            }
          }

          submit.ui.on('click', ()=>{
            if(controls.editMode == false){
              checkAndAddStyle(); 
            }
            else {
              var id = stid
              styleSpecForm.getValue();
              stateManager.addStyle(styleFormControl.value, sid, stid);
              finalizeStyle(id);
            }
          });

          cancel.ui.on('click', ()=>{
            if(controls.editMode){
              parameters.hide();
            }
            else {
              boundingBox.detach();
              delete this;
            }
          });

          boundingBox.on('keyup', (e)=>{
            if(e.key == 'Enter'){
              submit.ui.trigger('click');
            }
          });

        }

        plusButton.ui.on('click', ()=>{
          if( !_editingForm){
            var newStyle = new Style(sid);
            styles.append(newStyle.ui);
            _editingForm = true;
          }
        });   
      }

      function MovieBar(){
        var boundingBox = this.ui = $('<div></div>');
        var slide = new slider({width: 120 });
        var controlButtons = $('<div></div>');

        boundingBox.append(slide.ui);
        boundingBox.append(controlButtons);

        var play = new button(icons.movie.play, 20);
        var stop = new button(icons.movie.stop, 20);
        var previous = new button(icons.movie.previous, 20);
        var next = new button(icons.movie.next, 20);

        controlButtons.append(previous.ui);
        controlButtons.append(play.ui);
        controlButtons.append(stop.ui);
        controlButtons.append(next.ui);

        // Style
        boundingBox.css('position','absolute');
        boundingBox.css('padding','3px');
        boundingBox.css('text-align', 'center');
        // boundingBox.css('width','');



        // state variable
        var playing = false;

        play.ui.on('click', ()=>{
          // console.log('trying to change the background', play);
          if(playing){
            play.setSVG(icons.movie.play);
          }
          else {
            play.setSVG(icons.movie.pause);
          }
          playing = !playing;
        });

        stop.ui.on('click', ()=>{

        });

        previous.ui.click(()=>{

        });

        next.ui.click(()=>{

        });

      }

      function AlertBox(config){
        config = config || {};
        var width = config.width || 400;
        var timeout = config.timeout || 5000;
    
        var boundingBox = this.ui = $('<div></div>');
        var displayBox = $('<div></div>');
        boundingBox.append(displayBox);
    
        boundingBox.css('border-radius', '2px');
        // boundingBox.height(100);
        boundingBox.width(100);
        displayBox.css('border-radius', '2px');
        boundingBox.css('position', 'absolute');
        boundingBox.css('padding', '10px');
        displayBox.css('padding', '10px');
        boundingBox.hide();
        
        this.alert = function(type, message){
          boundingBox.show();

          if(type=='success'){
            displayBox.css('border', '1px solid green');
            displayBox.css('background', 'lightgreen');
            displayBox.css('color', 'green');
          } else if( type == 'warning'){
            displayBox.css('border', '1px solid orange');
            displayBox.css('background', 'yellow');
            displayBox.css('color', 'orange');
          } else if(type == 'error'){
            displayBox.css('border', '1px solid red');
            displayBox.css('background', 'lightred');
            displayBox.css('color', 'red');
          }

          displayBox.text(message);
    
          setTimeout(()=>{ 
            boundingBox.hide(); 
            // console.log('Hiding');
          }, timeout);
        }
    
      }
    
      function DialogBox(config){
        config = config || {};
        var height = config.height || 400;
        var width = config.width || 400;
    
        
        var boundingBox = this.ui = $('<div></div>');
        var mainForm = null;

        var semiDisplayBox = $('<div></div>');
        boundingBox.append(semiDisplayBox);
        
        var displayBox = $('<div></div>');
        boundingBox.append(displayBox);

        var controlBar = $('<div></div>');
        boundingBox.append(controlBar);
    
        var done = new button(icons.tick, 20);
        var drop = new button(icons.cross, 20);
        controlBar.width(20);
        controlBar.append(done.ui);
        controlBar.append(drop.ui);
        done.ui.css('display', 'block');
        done.ui.css('background', 'green');
        done.ui.css('margin-bottom', '5px');
        done.ui.css('margin-top', '5px');
        
        drop.ui.css('display', 'block');
        drop.ui.css('background', 'red');
        // drop.ui.css('background', 'red');
        
        controlBar.css('position', 'absolute');
        controlBar.css('top', '0px');
        // controlBar.css('display', 'flex');
        // controlBar.css('flex-direction', 'row');
    
        boundingBox.css('background', 'lightgrey');
        boundingBox.height(height);
        boundingBox.width(width);
        boundingBox.css('position', 'absolute');
        boundingBox.css('overflow-y', 'auto');
        boundingBox.css('border-radius', '4px');
        boundingBox.css('opacity', '0.4');
        boundingBox.css('padding', '6px');
    
        displayBox.css('display', 'flex');
        displayBox.css('flex-direction', 'row');
        displayBox.css('align-items', 'center');
        displayBox.css('position', 'relative');
    
        // semiDisplayBox.css('display', 'flex');
        // semiDisplayBox.css('flex-direction', 'row');
        semiDisplayBox.css('text-align', 'center');
        // semiDisplayBox.css('position', 'relative');
    
        boundingBox.on('mouseenter', ()=>{
          boundingBox.css('opacity', '0.8');
        });
    
        boundingBox.on('mouseleave', ()=>{
          boundingBox.css('opacity', '0.4');
        });
    
        done.ui.on('click', ()=>{
          if(mainForm.validate()){
            boundingBox.hide();
            console.log("Form Validated in the UI", mainForm.validate());
            // console.log("UI:DialogBox:Finalize", mainForm)
            stateManager.finalize(mainForm.getValue(), mainForm);
            // semiDisplayBox.html(null);
            // displayBox.html(null);
            mainForm.ui.detach();
            console.log("UI::DialogBox:done.ui.on->click:mainForm", mainForm, mainForm.getValue());
          }
          else{
            $(document).trigger('Error', ['Incorrect Input', 'Please check the input provided for selection, and resubmit the form']);
          }

        });
    
        drop.ui.on('click', ()=>{
          boundingBox.hide();
          stateManager.cancel();
          mainForm.ui.detach();
          mainForm = null;
        });
    
        var addForm = this.addForm = function(form, formType=null){
            form.ui.detach();  
            // displayBox.children().detach();
            displayBox.append(form.ui);
            form.ui.css('margin', 'auto');
            boundingBox.show();
            mainForm = form;
            // console.log("Form", mainForm);
          
        }
    
      }

      function slider(config){
        config = config || {};
        min = config.min || 0;
        max = config.max || 100;
        width = config.width || 400;
        
        var boundingBox = this.ui = $('<div></div>');
        var slide = $('<input type="range" min="0">');
  
        boundingBox.append(slide);
      }
    

      
      function ContextMenu(){
        var boundingBox = this.ui = $('<div></div>');

        boundingBox.css('position', 'absolute');
        // boundingBox.css('border', '1px solid black');
        boundingBox.css('border-radius', '3px');
        boundingBox.css('background', '#f1f1f1');
        boundingBox.css('z-index', 99);
        var contentBox = $('<div></div>');
        contentBox.css('position', 'relative');
        boundingBox.css('opacity', '0.85');

        boundingBox.append(contentBox);
        contentBox.css({
          'background':  '#f1f1f1',
          'border-radius': '4px',
          'padding': '4px',
          'width':'140px'
        });
        // Context Box
        // Remove Label Button 
        
        var labelMenuStyle = {
          'background' : '#d3e2ee',
          'padding' : '2px',
          'font-family':'Arial',
          'font-weight':'bold',
          'font-size':'12px',
          'border-radius':'2px',
          // 'margin-top':'3px'
        }

        var removeLabelMenu = $('<div></div>');
        removeLabelMenu.text('Remove Label');
        removeLabelMenu.css(labelMenuStyle);
        removeLabelMenu.css('margin-bottom', '3px');

        contentBox.append(removeLabelMenu);
        removeLabelMenu.hide();

        // Label Property List 
        var propertyKeys = Object.keys($3Dmol.GLModel.validAtomSpecs);
        var propertyList = [];
        var propertyObjectList = [];

        propertyKeys.forEach((prop)=>{
          var propObj = $3Dmol.GLModel.validAtomSpecs;
          if(propObj[prop].prop === true){
            propertyList.push(prop);
          }
        });

        // Property Menu 
        var propertyMenu = $('<div></div>');
        contentBox.append(propertyMenu);
        
        function Property(key, value){
          this.row = $('<tr></tr>');
          var propLabelValue = this.control = {
            key : '',
            value : null,
            active : true,
            name : key,
          }

          this.key = key;
          this.value = value;

          var checkbox = new $3Dmol.UI.Form.Checkbox(propLabelValue);
          var checkboxHolder = $('<td></td>');
          checkboxHolder.append(checkbox.ui); 
          var keyHolder = $('<td></td>');
          var separatorHolder = $('<td></td>').text(':');
          var valueHolder = $('<td></td>');

          this.row.append(checkboxHolder, keyHolder, separatorHolder, valueHolder);

          keyHolder.text(key);
          valueHolder.text(value);
        }

        function setProperties(atom){        
          propertyMenu.empty();
          propertyObjectList = [];
          
          var propertyTable = $('<table></table>');
          
          propertyList.forEach((prop)=>{
            var propObj = new Property(prop, atom[prop]);
            propertyTable.append(propObj.row);
            propertyObjectList.push(propObj);
          });

          propertyMenu.append(propertyTable);
          
          var submit = new button(icons.tick, 18, { backgroundColor: 'lightgreen'});
          var cancel = new button(icons.cross, 18, { backgroundColor: 'indianred'});

          var controlButtons = $('<div></div>');
          controlButtons.append(submit.ui, cancel.ui);
          // controlButtons.css('text-align', 'center');

          var alertBox = $('<div></div>');
          propertyMenu.append(alertBox);
          alertBox.css({
            'color': 'darkred',
            'border':'1px solid darkred',
            'border-radius': '3px',
            'background-color':'lightcoral',
            'padding':'3px',
            'text-align':'center',
            'font-family':'Arial',
            'font-size':'12px',
            'font-weight':'bold'
          });

          alertBox.hide();

          propertyMenu.append(controlButtons);


          submit.ui.on('click', ()=>{
            var props = processPropertyList();
            if(props !=null){
              stateManager.addAtomLabel(props, atom);
              stateManager.exitContextMenu(false);
            }
            else {
              alertBox.show();
              alertBox.text('No value selected for label');
            }
          });

          cancel.ui.on('click', ()=>{
            stateManager.exitContextMenu();
          });
        }

        // Previous Labels 
        var labelHolder = $('<div></div>');
        contentBox.append(labelHolder);

        // Add Menu 
        var addMenu = $('<div></div>');
        contentBox.append(addMenu);
        addMenu.css('width', '100%');

        var addLabelMenu = $('<div></div>');
        addMenu.append(addLabelMenu);


        addLabelMenu.text('Add Label');
        addLabelMenu.css(labelMenuStyle);
        addLabelMenu.css('margin-bottom', '3px');
        
        // Edit Menu
        var editMenu = $('<div></div>');
        contentBox.append(editMenu);

        contentBox.css({
          'position' : 'relative',
        });

        editMenu.css({
          'background':'#dfdfdf',
          'border-radius':'3px',
          'font-family':'Arial',
          'font-weight':'bold',
          'font-size':'12px',
          // 'position': 'absolute',
          // 'left' : '105%',
          // 'top' : '0',,
          'box-sizing' : 'border-box',
          'width' : '100%',
          
        });
        editMenu.hide();

        // Add Label Inputs 

        function generateAddLabelForm(){
          var addLabelForm = $('<div></div>');

          var addLabelValue = {
            text : {
              key : 'Label Text',
              value : null,
              active : true,
            },
  
            style : {
              key : 'Style',
              value : null,
              active : true,
            },
  
            sel : {
              key : 'Selection',
              value : null,
              active : true,
            }
          }
          var formModifierControl = $('<div></div>');
          var removeButton = new button(icons.minus, 16);
          var tick = new button(icons.tick, 16);
          var cross = new button(icons.cross, 16);
          formModifierControl.append(removeButton.ui, tick.ui, cross.ui);
          removeButton.ui.hide();
          addLabelForm.append(formModifierControl);
  
          var addLabelTextBox = $('<div></div>');
          var lt = $('<div></div>').text('Label Text');
          var addLabelTextInput = new $3Dmol.UI.Form.Input(addLabelValue.text);
          addLabelTextBox.append(lt, addLabelTextInput.ui);
          var width = 126//editMenu.innerWidth()*0.8;
          addLabelTextInput.setWidth(width);
          addLabelForm.append(addLabelTextBox);

          var addLabelStyleBox = $('<div></div>');
          var ls = $('<div></div>').text('Label Style');
          var addLabelStyleInput = new $3Dmol.UI.Form.ListInput(addLabelValue.style, Object.keys($3Dmol.labelStyles));
          addLabelStyleBox.append(ls, addLabelStyleInput.ui);
          addLabelForm.append(addLabelStyleBox);

          var selectionList = stateManager.getSelectionList();
          
          var addLabelSelectionBox = $('<div></div>');
          var lsl = $('<div></div>').text('Label Selection');
          var addLabelSelectionInput = new $3Dmol.UI.Form.ListInput(addLabelValue.sel, selectionList);
          addLabelSelectionBox.append(lsl, addLabelSelectionInput.ui);
          addLabelForm.append(addLabelSelectionBox);

          // CSS 
          addLabelForm.css({
            'padding' : '2px',
            
          });

          tick.ui.on('click', ()=>{
            var validate = true;

            if(!addLabelStyleInput.validate())
              validate = false;
            
            if(!addLabelTextInput.validate())
              validate = false;
            
            if(validate){
              stateManager.addLabel(addLabelValue);
            } else {
              console.log('Please Check the input');
            }       
          });

          cross.ui.on('click', ()=>{
            stateManager.exitContextMenu();
          });

          removeButton.ui.on('click', ()=>{
            stateManager.removeLabel()
          });

          addLabelForm.on('keyup', (e)=>{
            if(e.key == 'Enter'){
              tick.ui.trigger('click');
            }
          });
      
          return {
            boundingBox : addLabelForm,
            text : addLabelTextInput,
            style : addLabelStyleInput,
            selection: addLabelSelectionInput,
            editMode : function(){
              removeButton.ui.show();
            }
          }
        }


        function processPropertyList(){
          var propsForLabel = {};
          
          propertyObjectList.forEach((propObj)=>{
            if(propObj.control.value === true){
              propsForLabel[propObj.key] = propObj.value;
            }
          });

          if(Object.keys(propsForLabel).length != 0){
            return propsForLabel
          }
          else{
            return null;
          }
        }

        // Context Menu UI Funciton 
        boundingBox.hide();
        this.hidden = true;
        this.atom = null;

        removeLabelMenu.on('click', { atom : this.atom }, function(e){
          stateManager.removeAtomLabel(removeLabelMenu.atom);
          // console.log(removeLabelMenu.atom, "Atom to remove");
        });


        this.show = function(x, y, atom, atomExist){
          console.log('Context Menu open and atom Exist', atomExist);

          if(atomExist){
            removeLabelMenu.show();
            removeLabelMenu.atom = atom;
          }
          else {
            removeLabelMenu.hide();
            removeLabelMenu.atom = null;
          }

          unsetForm();
          setPosition(boundingBox, x, y);
          console.log('CONTEXT MENU::Atom Selected', atom);
          boundingBox.show();
          this.hidden = false;
          
          if(atom){
            setProperties(atom);
            this.atom = atom;
          }
          else{

            propertyMenu.empty();
          }
        }
        
        this.hide = function(processContextMenu){
          if(processContextMenu){
            var propsForLabel = processPropertyList();
            if(propsForLabel != null){
              stateManager.addAtomLabel(propsForLabel, this.atom);
              // console.log("These property will be used to add label", propsForLabel);
            }
          }

          boundingBox.hide();
          this.hidden = true;
          unsetForm();
        }

        addLabelMenu.on('click', function(){
          var addLabelMenuForm = generateAddLabelForm();
          setForm(addLabelMenuForm);
        });

        function setForm(form){
          editMenu.children().detach();
          editMenu.append(form.boundingBox);
          editMenu.show();
        }

        function unsetForm(){
          editMenu.children().detach();
          editMenu.hide();
        }
      }


      function SurfaceMenu(){
        var boundingBox = this.ui = $('<div></div>');
        var _editingForm = false;
        // Selection Layout

        // var contentBox = $('<div></div>');
        // contentBox.css('background', 'red');
        // boundingBox.append(contentBox);

        boundingBox.css({
          'position': 'absolute',
          'width': '140px',
          'text-align': 'right'   
        });

        // contentBox.css({
        //   'position':'relative',
        //   'right':'100%',
        //   'text-align':'right'
        // })

        
        var surfaceButton = new button(icons.surface, 20, { tooltip : 'Open Surface Menu'});

        boundingBox.append(surfaceButton.ui);


        var displayBox = $('<div></div>');
        boundingBox.append(displayBox);

        // Overflow fix 
        boundingBox.css({
          'overflow':'hidden',
        });

        console.log('Surface box display box height', HEIGHT);

        displayBox.css({
          'max-height' : HEIGHT*0.8,
          'overflow-y': 'auto',
          'overflow-x' : 'hidden'
        });

        var newSurfaceSpace = $('<div></div>');
        
        // newSurfaceSpace.append(controlButton);
        // controlButton.hide();

        displayBox.append(newSurfaceSpace);

        var addArea = $('<div></div>');
        var addButton = new button(icons.plus, 20);
        addArea.append(addButton.ui);
        displayBox.append(addArea);
        displayBox.hide();

        var surfaces = this.surfaces = [];
        var currentSurface = this.currentSurface = null;

        function Surface(){
          var control = {
            surfaceType: {
              key : 'Surface Type',
              value : null
            },
            surfaceStyle: {
              key: 'Surface Style',
              value: null
            },
            surfaceFor: {
              key: 'Selection Atoms',
              value: null
            },
            surfaceOf: {
              key: 'Surface Generating Atoms',
              value: null,
            },
          };
  
          var surfaceBox = this.ui = $('<div></div>');
          surfaceBox.css({
            'margin-top':'3px',
            'padding':'6px',
            'border-radius':'3px',
            'background-color': '#e8e8e8',
            // 'position':'relative',
            'width':'100%',
            'box-sizing':'border-box',
            // 'left': "-100%",
            'opacity' : 0.9,
            'text-align':'left'
          });
        
          var heading = this.heading = $('<div></div>');
          var header = $('<div></div>');

          header.css({
            'text-align' : 'right'
          })

          // Control Buttons
          var toolButtons = $('<div></div>');
          
          var editButton = new button(icons.pencil, 16);
          var removeButton = new button(icons.minus, 16, { bfr:0.5, backgroundColor:'#f06f6f'});
          
          toolButtons.append(removeButton.ui);
          toolButtons.append(editButton.ui);

          toolButtons.editButton = editButton;
          toolButtons.removeButton = removeButton;
          toolButtons.editMode = false;
          
          var defaultTextStyle = {
            'font-weight': 'bold',
            'font-family' : 'Arial',
            'font-size' : '12px'
          }
          
          heading.css('display', 'inline-block');
          heading.css(defaultTextStyle);
          
          toolButtons.css('display', 'inline-block');
          header.hide();
          
          header.append(heading, toolButtons);
          surfaceBox.append(header);

          // toolButtons.hide();
          var surfacePropertyBox = $('<div></div>');
          surfaceBox.append(surfacePropertyBox);

          // Surface Type
          var surfaceType = $('<div></div>');

          var labelSurfaceType = $('<div></div>');
          labelSurfaceType.text('Surface Type');
          labelSurfaceType.css(defaultTextStyle);

          var listSurfaceType =new $3Dmol.UI.Form.ListInput(control.surfaceType, Object.keys($3Dmol.SurfaceType));

          surfaceType.append(labelSurfaceType, listSurfaceType.ui);
          surfacePropertyBox.append(surfaceType);
          
          // Surface Style
          var surfaceStyle = $('<div></div>');

          var labelSurfaceStyle = $('<div></div>');
          // labelSurfaceStyle.text('Surface Style');

          var formSurfaceStyle = new $3Dmol.UI.Form($3Dmol.GLModel.validSurfaceSpecs, control.surfaceStyle);
          
          surfaceStyle.append(labelSurfaceStyle, formSurfaceStyle.ui);
          surfacePropertyBox.append(surfaceStyle);
          
          // Surface Of
          var surfaceOf = $('<div></div>');
          
          var labelSurfaceOf = $('<div></div>');
          labelSurfaceOf.text('Surface Atoms');
          labelSurfaceOf.css(defaultTextStyle);
          
          var surfaceGeneratorAtomType = ['self', 'all'];
          
          var listSurfaceOf = new $3Dmol.UI.Form.ListInput(control.surfaceOf, surfaceGeneratorAtomType);
          
          surfaceOf.append(labelSurfaceOf, listSurfaceOf.ui);
          surfacePropertyBox.append(surfaceOf);
          
          // Surface For
          var selectionListElement = stateManager.getSelectionList();
          var surfaceFor = $('<div></div>');
          
          var labelSurfaceFor = $('<div></div>');
          labelSurfaceFor.text('Show Atoms');
          labelSurfaceFor.css(defaultTextStyle);

          var listSurfaceFor = new $3Dmol.UI.Form.ListInput(control.surfaceFor, selectionListElement);
          
          surfaceFor.append(labelSurfaceFor, listSurfaceFor.ui);
          surfacePropertyBox.append(surfaceFor);
          
          // Control Button
          var controlButton = $('<div></div>');
          var submit = new button(icons.tick, 20);
          var cancel = new button(icons.cross, 20);
          controlButton.append(submit.ui);
          controlButton.append(cancel.ui);
          surfacePropertyBox.append(controlButton);

          // Functionality 
          removeButton.ui.on('click', { surfaceBox : surfaceBox }, function(e){
            var id = e.data.surfaceBox.data('surf-id');
            surfaceBox.remove();
            stateManager.removeSurface(id);
          });

          editButton.ui.on('click', function(){
            surfacePropertyBox.toggle();
            
            // After creation of the surface box all the changes will be edit to the surfaces so on first submit toolButtons.editMode == true;
          });

          


          // Form Validation 

          var validateInput = this.validateInput = function(){
            var validated = true;
            
            if( listSurfaceType.getValue().value == 'default'){
              validated = false;
              console.log('Surface::InvalidSurfaceType');
            }

            if( listSurfaceOf.getValue().value == 'default'){
              validated = false;
              console.log('Surface::InvalidSurfaceAtom Selected');
            }

            if(listSurfaceFor.getValue().value == 'default'){
              validate = false;
              console.log('Surface::InvalidSurfaceGenerationAtom Selected');
            }

            return validated;
          }

          

          // boundingBox.on('mouseenter', function(){
          //   selections = stateManager.getSelectionList();
          //   selectionListElement = selections.map( (m)=>{
          //     return m.id;
          //   });
          //   listSurfaceFor.updateList(selectionListElement);
          //   listSurfaceOf.updateList(selectionListElement);
          // });

          var updateSelection = this.updateSelection = function(){
            selections = stateManager.getSElectionList();
            selectionListElement = selections.map( (m)=>{
              return m.id;
            });
            listSurfaceFor.updateList(selectionListElement);
            listSurfaceOf.updateList(selectionListElement);
          }

          // Submit 
          submit.ui.on('click', {}, function(){
            if(validateInput()){ 
              if(toolButtons.editMode === false){
                formSurfaceStyle.getValue();
                var id = stateManager.addSurface(control);
                control.id = id;

                // element properties
                surfaceBox.data('surf-id', id);
                
                heading.text('surf#' + id);

                header.show();
                toolButtons.editMode = true;
                surfacePropertyBox.hide();

                surfaces.push(this);
                _editingForm = false;
              }
              else{
                console.log('Edit Surface called');

                formSurfaceStyle.getValue();
                control.id = surfaceBox.data('surf-id');
                stateManager.editSurface(control); // -> add updateSurface funciton to surfaceMenu
                surfacePropertyBox.hide();
              }
            }
          });

          // Cancel Edit
          cancel.ui.on('click', {}, function(){
            if(toolButtons.editMode == false){
              surfaceBox.detach();
              surfaceBox.remove();
              _editingForm = false;
              // Add statemanager handle to remove the object itself
            }
            else {
              surfacePropertyBox.hide();
              toolButtons.editMode = false;
            }
          });

          surfaceBox.on('keyup', (e)=>{
            if(e.key == 'Enter'){
              submit.ui.trigger('click');
            }
          });

        }

        // Functionality

        // Surface addition

        addButton.ui.on('click', { surfaces: this }, function(e){
          if(!_editingForm){
            var newSurface = new Surface();
            newSurfaceSpace.append(newSurface.ui);
            _editingForm = true;
          }
        });

        surfaceButton.ui.on('click', ()=>{
          displayBox.toggle();
        });

      }

      /**
       * position : Sets the css position property : absolute
       * @param {object} jquery html element
       * @param {number} left : css left property
       * @param {number} top : css top peroperty
       */
      function setPosition(ele, left, top){
        // ele.css('position', 'absolute');
        ele.css('left', left);
        ele.css('top', top);
      }


      /**
        * setLocation - Sets the location of the element relative to the parseInt
        * as per position types
        *
        * @param  {object} parent jquery object
        * @param  {object} child  jquery object
        * @param  {string}   x_type 'left|right'
        * @param  {string}   y_type 'top|bottm'
        *
        */
      function setLocation(parent, child, x_type='left', y_type='top', x_offset=0, y_offset=0){
        // console.log('Setting location', parent.offset(), child.offset(), parent.height(), parent.width(), child.height(), child.width());

        // p_ stands for parent
        child.css('z-index', 99);


        var p_position = parent.position();
        var p_width = getWidth(parent);
        var p_height = getHeight(parent);

        // c_ stand for child
        var c_width = child.outerWidth(); // includes padding and margin
        var c_height = child.outerHeight(); // includes padding and margin

        var padding = parseInt(parent.css('padding').replace('px', ''));
        padding = (padding)? padding: 0;
        console.log("Checking Padding", padding == NaN, padding)
        var p_top = getTop(parent) + parseInt(parent.css('margin-top').replace('px',''));
        var p_left = getLeft(parent) + parseInt(parent.css('margin-left').replace('px',''));


        // Setting position
        var c_position = {
          left: 0,
          top: 0
        };

        if(x_type == 'left'){
          // console.log('left is called');
          c_position.left = padding + x_offset;
        }
        else if(x_type == 'center'){
          c_position.left = p_width/2 - c_width/2 + x_offset;
        }
        else if(x_type == 'right'){
          c_position.left = p_width  - c_width - padding + x_offset;
        }
        else {
          c_position.left = x_offset + padding;
        }

        if(y_type == 'top'){
          c_position.top = y_offset + padding;
        }
        else if(y_type == 'center'){
          c_position.top = p_height/2 - c_height/2 + y_offset;
        }
        else if(y_type == 'bottom'){
          c_position.top = p_height - c_height - y_offset - padding;
        }
        else {
          c_position.top = y_offset + padding;
        }

        setPosition(child, c_position.left, c_position.top);
        console.log('Setting Location', c_position, p_height, c_height, y_offset, padding);
      }

      // Copied from glviewer.js
      function getRect(container) {
        let div = container[0];
        let rect = div.getBoundingClientRect();
        if(rect.width == 0 && rect.height == 0 && div.style.display === 'none' ) {
          let oldpos = div.style.position;
          let oldvis = div.style.visibility;
          div.style.display = 'block';
          div.style.visibility = 'hidden';
          div.style.position = 'absolute';
          rect = div.getBoundingClientRect();
          div.style.display = 'none';
          div.style.visibility = oldvis;
          div.style.position = oldpos;
        }
        return rect;
      };

      function getTop(container) {
        return container.offset().top;
      }

      function getLeft(container) {
        return container.offset().left;
      }

      function getHeight(container) {
        return getRect(container).height;
      }

      function getWidth(container) {
        return getRect(container).width;
      }

      /**
        * button - generates button with the given markup as contents
        * @param {string} SVG markup string that contains the content of the button
        * @param {number} Height of the content
        * @param {}
        * @return {object}  Returns the jquery object for button
        */
      function button(svg, height, config){
        config = config || {};
        var borderRadius = config.bfr*height || (height/4); // body radius factor
        var bgColor = config.backgroundColor || 'rgb(177, 194, 203)';
        var color = config.color || 'black'; 
        var hoverable = config.hoverable || 'true';
        var tooltipText = config.tooltip || null;

        // Button instance
        var button = this.ui = $('<div></div>');
        var innerButton = $('<div></div>');
        button.append(innerButton);

        // CSS
        button.css('box-sizing', 'border-box');
        button.css('display', 'inline-block');
        button.css('margin', '3px');
        button.css('height', height);
        button.css('width', height);
        button.css('border-radius', borderRadius + 'px');
        
        //  button.css('padding', '3px');
        button.css('color', color);
        button.css('background', bgColor);

        innerButton.css('display', 'flex');
        innerButton.css('justify-content','center');
        innerButton.css('align-items', 'center');
        innerButton.css('padding', '2px');
        
        // content
        this.setSVG = function(svg){
          innerButton.empty();
          var formatted_content = $(svg);
          innerButton.append(formatted_content);

          }

        this.setSVG(svg);

        // Hover
        
        // Setting up tool tip
        button.css({
          'position' : 'relative'
        });


        // Tool tips will only be created if tooltip description is set during creation
        var tooltipBox = null

        if(tooltipText != null ){
          tooltipBox = $('<div></div>');
          tooltipBox.text(tooltipText);
  
          button.append(tooltipBox);
          tooltipBox.css({
            'background-color': 'black',
            'color': '#fff',
            'text-align': 'center',
            'border-radius': '6px',
            'padding': '4px', 
            'width' : '120px',

            'font-size':'10px',
            'font-family':'Arial',
            
            'position' : 'absolute',
            'top': '0%',
            'left':'110%',
            'z-index':'100'
          });

          tooltipBox.hide();

          button.on('hover', ()=>{
            tooltipBox.toggle();
          })
        }

        if(hoverable == 'true'){
          button.on('mouseenter',
            ()=>{
              button.css('box-shadow', '0px 0px 3px black');
              if(tooltipText != null){
                tooltipBox.show();

                setTimeout(()=>{
                  tooltipBox.hide();
                }, 3000);
              }
            }).on('mouseleave', 
            ()=>{
              button.css('box-shadow', 'none');
              if(tooltipText != null) {
                tooltipBox.hide();
              }
            }
          );
            
          var longPressTime = 0;
          var mouseX = 0;
          var mouseY = 0;

          // click
          button.on('mousedown', (e)=>{
            button.css('box-shadow', '0px 0px 1px black');
            // var startX = e.clientX;
            // var startY = e.clientY;

            // setTimeout(()=>{
            //   if(mouseX == startX && mouseY == startY){
            //     console.log('Show tooltip');
            //   }
            //   else{
            //     console.log('mouse moved');
            //   }
            // }, $3Dmol.longPressDuration);
          });
  
          button.on('mouseup', ()=>{
            button.css('box-shadow', '0px 0px 3px black');
          });

          button.on('mousemove', (e)=>{
            // mouseX = e.clientX;
            // mouseY = e.clientY;
          });

          var timer;
          button.on('touchstart', (e)=>{
            timer = setTimeout(()=>{
              (tooltipText != null)? tooltipBox.show() : null;
            }, $3Dmol.longPressTime)
          });

          button.on('touchend', ()=>{
            (tooltipText != null)? tooltipBox.hide() : null;
            if(timer)
              clearTimeout(timer);
          })

        }



      }

    }

    return UI;
  })();
