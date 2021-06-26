/*
* This is ui that we are trying to draw at the time of generationo
* of viewer.
* Refer to line number 43 of the glViewer to see the viewer code
*/

/**
 * $3Dmol.UI - UI creates panels of viewer to assist control of the viewport in
 * 3dmol.js
 *
 * @param  {object} contains selector for viewer element
 * @return {object} contains the configuration for making the viewer with
 * different sets of controls
 */
$3Dmol.UI = (function(){
  function UI(stateManager, config){
    // Extract the viewer and then render it
    var icons =new $3Dmol.UI.Icons();
    var body = $('body');
    
    // Generates the necessary UI elements
    var uiElements = this.tools = generateUI(config);
    
    function generateUI(config){
      
      
      var ui_overlay = new UI_Overlay(config);
      body.append(ui_overlay.ui);
      // body = ui_overlay.ui;


      var topbar = new Toolbar();
      body.append(topbar.ui);
      setLocation(ui_overlay.ui, topbar.ui, 'left', 'top');

      var selectionBox = new SelectionBox(icons.select);
      setLocation(ui_overlay.ui, selectionBox.ui, 'left', 'top', 0, topbar.ui.outerHeight());
      body.append(selectionBox.ui);

      // var styleBox = new SelectionBox(icons.paintbrush);
      // setLocation(ui_overlay.ui, styleBox.ui, 'left', 'top', 0, topbar.ui.outerHeight() + selectionBox.ui.outerHeight());
      // body.append(styleBox.ui);

      var movieControl = new MovieBar();
      body.append(movieControl.ui);
      setLocation(ui_overlay.ui, movieControl.ui, 'center', 'bottom');


      var dialog = new DialogBox({ height: 300, width: 300 });
      body.append(dialog.ui);
      setLocation(ui_overlay.ui, dialog.ui, 'center', 'center');
      dialog.ui.hide();
      
      // Testing form
      // var selectionForm = new $3Dmol.UI.form($3Dmol.GLModel.validAtomSelectionSpecs, 'Atom Selection');
      // dialog.addForm(selectionForm);

      var alertBox = new AlertBox({ width : 100 });
      // console.log(alertBox);
      body.append(alertBox.ui);
      setLocation(ui_overlay.ui, alertBox.ui, 'right', 'top');

      return {
        uiOverlay : ui_overlay,
        topbar : topbar,
        selectionBox : selectionBox,
        dialog : dialog,
        alertBox : alertBox,
        movieControl : movieControl
      } 
    }

    function UI_Overlay(config){
      config = config || {};

      var top = (config.offset.top !=undefined) ? config.offset.top : 10;
      var left = (config.offset.left !=undefined) ? config.offset.left : 10;
      var width = config.width || '400px';
      var height = config.height || '400px';

      var ui_overlay = this.ui = $('<div></div>');
      ui_overlay.height(height);
      ui_overlay.width(width);
      setPosition(ui_overlay, left, top);
      ui_overlay.css('box-sizing','border-box');
      // ui_overlay.css('border', '1px solid black');
      ui_overlay.css('background', 'rbga(0,0,0,0)');
      ui_overlay.css('padding', '0px');
      ui_overlay.css('pointer-events', 'none');
      ui_overlay.css('position', 'absolute');

      this.resize = function(config){
        var top = (config.offset.top !=undefined) ? config.offset.top : 10;
        var left = (config.offset.left !=undefined) ? config.offset.left : 10;
        var width = config.width || '400px';
        var height = config.height || '400px';
        ui_overlay.height(height);
        ui_overlay.width(width);
        setPosition(ui_overlay, left, top);
      }

    }

    this.resize = function(config){
      uiElements.uiOverlay.resize(config);
      this.orient();
    }

    this.orient = function(){
      // Position of the UI elements can later be orgnazided by using a config


      var ui_overlay = uiElements.uiOverlay;
      var topbar = uiElements.topbar;
      var dialog = uiElements.dialog;
      var alertBox = uiElements.alertBox;
      var selectionBox = uiElements.selectionBox;
      var movieControl = uiElements.movieControl;


      setLocation(ui_overlay.ui, alertBox.ui, 'right', 'top');
      setLocation(ui_overlay.ui, dialog.ui, 'center', 'center');
      setLocation(ui_overlay.ui, movieControl.ui, 'center', 'bottom');
      setLocation(ui_overlay.ui, selectionBox.ui, 'left', 'top', 0, topbar.ui.outerHeight());
      setLocation(ui_overlay.ui, topbar.ui, 'left', 'top');
    }



    /**
    * function to make the html toolbar div
    * @param {object} config : Stores the top, left, width, padding property1
    * for the toolbar
    **/
    function Toolbar(config){
      config = config || {};

      var top = config.top || 10;
      var left = config.left || 10;
      // var width = config.width || '300px';
      var padding = config.padding || '3px';
      // var height = config.height || '20';

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
      // toolbar.css('width', width);
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
     * SelectionBox - Draws the box where all the selections on atoms are listed
     * This will be used to modify style for specific set of selection
     *
     * @return {Object}  Jquery element of div
     */
    function SelectionBox(icon, side='left') {
      var selectionBox = this.ui = $('<div></div>')
      var selections = $('<div></div>');
      var showArea = $('<div></div>');
      var addArea = $('<div></div>');
      var plusButton = new button(icons.plus, 20);
      var hideButton = new button(icon, 20);
      var selectionObjects = [];

      // Content
      selectionBox.append(hideButton.ui);
      selectionBox.append(showArea);
      selectionBox.css('position', 'absolute');

      showArea.append(selections);
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

      selectionBox.css('padding', '3px');

      showArea.css('box-sizing', 'border-box');
      showArea.css('padding', '3px');
      showArea.css('width', '120px');
      showArea.css('height', '300px');
      showArea.css('border-radius', '6px');
      showArea.css('background', 'rgb(214, 214, 214)');
      selections.css('margin-bottom', '10px');
      addArea.css('text-align', 'center');

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

      // Selection Box modification function

      /**
       * Add a selection input to the list
       */
      this.appendSelection = function(selectionSpec){
        var selection = $('<div></div>');
        var controls = $('<div></div>');
        var styleName = $('<div></div>');
        var selectionName = $('<div></div>');
        selectionName.text(selectionSpec.id);
        selection.append(selectionName);

        selection.append(controls);
        selection.append(controls);

        var editButton = new button(icons.edit, 20);
        var removeButton = new button(icons.remove, 20);
        var hideButton = new button(icons.list, 20);
        var showButton = new button(icons.style, 20);

        controls.append(editButton.ui);
        controls.append(removeButton.ui);
        controls.append(hideButton.ui);
        controls.append(showButton.ui);


        spec = selectionSpec.spec;
        var keys = Object.keys(spec);
        keys.forEach((key)=>{
          var text = `${key} : ${spec[key]}`
          var div = $('<div></div>');
          div.text(text);
          selection.append(div);
          console.log(spec[key])
        });

        selections.append(selection);
        selection.data('sel-id', selectionSpec.id);
        
        selectionObjects.push(selection);
      }

      plusButton.ui.on('click', ()=>{
        stateManager.addSelection();
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
        console.log('trying to change the background', play);
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
          console.log('Hiding');
        }, timeout);
      }
  
    }
  
    function DialogBox(config){
      config = config || {};
      var height = config.height || 400;
      var width = config.width || 400;
  
      var boundingBox = this.ui = $('<div></div>');
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
  
      boundingBox.on('mouseenter', ()=>{
        boundingBox.css('opacity', '0.8');
      });
  
      boundingBox.on('mouseleave', ()=>{
        boundingBox.css('opacity', '0.4');
      });
  
      done.ui.on('click', ()=>{
        boundingBox.hide();
        stateManager.finalize(this.form.getValues());
        // stateManager.
      });
  
      drop.ui.on('click', ()=>{
        boundingBox.hide();
        stateManager.cancel();
      });
  
      this.addForm = function(form){
        displayBox.empty();
        displayBox.append(form.ui);
        form.ui.css('margin', 'auto');
        boundingBox.show();

        this.form = form;
      }
  
    }
  
    // HTML Elements for the viewer
      /**
       * button - generates button with the given markup as contents
       * @param {string} SVG markup string that contains the content of the button
       * @param {number} Height of the content
       * @param {}
       * @return {object}  Returns the jquery object for button
       */
  
       function button(svg, height){
        // Button instance
        var button = this.ui = $('<div></div>');
  
        // CSS
        button.css('box-sizing', 'border-box');
        button.css('display', 'inline-block');
        button.css('margin-right', '5px');
        button.css('height', height);
        button.css('width', height);
        button.css('border-radius', (height/2) + 'px');
        button.css('padding', '3px');
        button.css('color', 'black');
        button.css('vertical-align', 'middle');
        button.css('text-align', 'center');
        button.css('background', 'rgb(177, 194, 203)');
  
        // Content Manipulation
        this.setSVG = function(svg){
          button.empty();
          var formatted_content = $(svg).height(height-6);
          button.append(formatted_content);
        }
  
        // this.setSVG(svg);
  
        // Hover
        button.hover(
          ()=>{
              button.css('box-shadow', '0px 0px 3px black');
          },
          ()=>{
            button.css('box-shadow', 'none');
          }
        )
  
        // click
        button.mousedown(()=>{
          button.css('box-shadow', '0px 0px 1px black');
        });
  
        button.mouseup(()=>{
          button.css('box-shadow', '0px 0px 3px black');
        });
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
    /**
     * parseSelection - generate markup using the selection list
     *
     * @param {list} Takes in the list of objects that are to selected
     * @return {Objects}  Returns list of selection after preprocessing the
     * SelectionSpec and StyleSpec. This is called `ParsedSelectionObject`
     */
    function parseSelection(Selections){
      return [
        { 'selection-name': 'Name 1',
          'properties':[
            { 'name':'property1', 'value':'value1', 'type':'input' },
            { 'name':'property2', 'value':'value2', 'type':'radio' },
            { 'name':'property3', 'value':'value3', 'type':'checkbox' },
            { 'name':'property4', 'value':'value4', 'type':'list' },
          ],
          'style' : { 'name':'StyleType1', 'value':'value1', 'type':'input' }
        },
        { 'selection-name': 'Name 2',
          'properties':[
            { 'name':'property1', 'value':'value1', 'type':'input' },
            { 'name':'property2', 'value':'value2', 'type':'radio' },
            { 'name':'property3', 'value':'value3', 'type':'checkbox' },
            { 'name':'property4', 'value':'value4', 'type':'list' },
          ],
          'style' : { 'name':'StyleType1', 'value':'value1', 'type':'input' }
        },
      ];
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
       console.log('Setting location', parent.offset(), child.offset(), parent.height(), parent.width(), child.height(), child.width());

       // p_ stands for parent
       var p_position = parent.position();
       var p_width = getWidth(parent);
       var p_height = getHeight(parent);

       // c_ stand for child
       var c_width = child.outerWidth(); // includes padding and margin
       var c_height = child.outerHeight(); // includes padding and margin

       var padding = parseInt(parent.css('padding').replace('px', ''));

       var p_top = getTop(parent) + parseInt(parent.css('margin-top').replace('px',''));
       var p_left = getLeft(parent) + parseInt(parent.css('margin-left').replace('px',''));


      // Setting position
      var c_position = {
        left: 0,
        top: 0
      };

      if(x_type == 'left'){
        console.log('left is called');
        c_position.left = p_left + padding + x_offset;
      }
      else if(x_type == 'center'){
        c_position.left = p_left + p_width/2 - c_width/2 + x_offset;
      }
      else if(x_type == 'right'){
        c_position.left = p_left + p_width  - c_width - padding + x_offset;
      }
      else {
        c_position.left = p_left + x_offset + padding;
      }

      if(y_type == 'top'){
        c_position.top = p_top + y_offset + padding;
      }
      else if(y_type == 'center'){
        c_position.top = p_top + p_height/2 - c_height/2 + y_offset;
      }
      else if(y_type == 'bottom'){
        c_position.top = p_top + p_height - c_height + y_offset - padding;
      }
      else {
        c_position.top = p_top + y_offset + padding;
      }

       // var c_position = {
       //   left: (x_type == 'left' )? p_left + padding + x_offset : (x_type == 'right')? p_width - c_width - padding + x_offset : x_offset,
       //   top: (y_type == 'top' )? p_top + padding + y_offset : (y_type == 'bottom')? p_height - c_height - padding + y_offset : y_offset
       // };

       console.log('setting parent child location',p_position,  c_position, padding, p_width, p_height, c_width, c_height);
       setPosition(child, c_position.left, c_position.top);
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
     function button(svg, height){
       // Button instance
       var button = this.ui = $('<div></div>');
       // CSS
       button.css('box-sizing', 'border-box');
       button.css('display', 'inline-block');
       button.css('margin-right', '5px');
       button.css('height', height);
       button.css('width', height);
       button.css('border-radius', (height/2) + 'px');
       button.css('padding', '3px');
       button.css('color', 'black');
       button.css('background', 'rgb(177, 194, 203)');

       // content
       this.setSVG = function(svg){
         button.empty();
         var formatted_content = $(svg).height(height-6);
         button.append(formatted_content);
       }

       this.setSVG(svg);

       // Hover
       button.hover(
         ()=>{
             button.css('box-shadow', '0px 0px 3px black');
         },
         ()=>{
           button.css('box-shadow', 'none');
         }
       )

       // click
       button.mousedown(()=>{
         button.css('box-shadow', '0px 0px 1px black');
       });

       button.mouseup(()=>{
         button.css('box-shadow', '0px 0px 3px black');
       });

     }

  }

  return UI;
})();
