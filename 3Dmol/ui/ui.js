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
  function UI(element, config){
    // Extract the viewer and then render it
    var container = $(element);


    console.log('Container for viewer', container);


    // Generates the necessary UI elements
    var icons =new $3Dmol.UI.Icons();
    generateUI();

    function generateUI(){
      var ui_overlay = new UI_Overlay();

      var padding = 10;
      var top = getTop(container) + padding;
      var left = getLeft(container) + padding;
      var right = getLeft(container) + getWidth(container) - 10;
      var bottom = getTop(container) + getHeight(container) - 10;
      console.log('Get position', left, top, right, bottom);

      $('body').append(ui_overlay);


      var topbar =new Toolbar();
      ui_overlay.append(topbar);

      var selectionBox = new SelectionBox();
      setLocation(ui_overlay, selectionBox, 'right', 'top', -15, 15);
      // selectionBox.css('position', 'absolute');
      // selectionBox.css('right', 0);
      // selectionBox.css('top', top);
      ui_overlay.append(selectionBox);
    }

    function UI_Overlay(config){
      config = config || {};

      var top = config.top || 10;
      var left = config.left || 10;
      var width = config.width || '400px';
      var height = config.height || '400px';

      var ui_overlay = $('<div></div>');
      ui_overlay.height(height);
      ui_overlay.width(width);
      setPosition(ui_overlay, left, top);
      ui_overlay.css('border', '1px solid black');

      return ui_overlay;
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
      var width = config.width || '300px';
      var padding = config.padding || '3px';
      // var height = config.height || '20';

      var toolbar = $('<div></div>');

      // Toolbar content

      var optionButton = new button(icons.option, 20)
      var moveButton = new button(icons.move, 20);
      var rotateButton = new button(icons.rotate, 20);
      toolbar.append(optionButton);
      toolbar.append(moveButton);
      toolbar.append(rotateButton);

      // Toolbar design
      toolbar.css('border', 'none');
      toolbar.css('padding', padding);
      toolbar.css('background', 'none');
      toolbar.css('color', 'white');
      toolbar.css('width', width);
      toolbar.css('box-sizing', 'border-box');

      // Button Action
      optionButton.click(()=>{
        console.log('option button clicked');
      });

      moveButton.click(()=>{
        console.log('move button clicked');
      });

      rotateButton.click(()=>{
        console.log('click button clicked');
      });

      return toolbar;
    }

    /**
     * SelectionBox - Draws the box where all the selections on atoms are listed
     * This will be used to modify style for specific set of selection
     *
     * @return {Object}  Jquery element of div
     */
    function SelectionBox() {
      var selectionBox = $('<div></div>')
      var selections = $('<div></div>');
      var addArea = $('<div></div>');
      var plusButton = new button(icons.plus, 20);
      // Content
      selectionBox.append(selections);
      addArea.append(plusButton);
      selectionBox.append(addArea);

      // CSS
      selectionBox.css('box-sizing', 'border-box');
      selectionBox.css('padding', '3px');
      selectionBox.css('width', '120px');
      selectionBox.css('height', '300px');
      selectionBox.css('border-radius', '6px');
      selectionBox.css('background', 'rgb(214, 214, 214)');
      selections.css('margin-bottom', '10px');
      addArea.css('text-align', 'center');

      // Selection Box modification function

      /**
       * Add a selection input to the list
       */
      this.addSelection = function(ParsedSelectionObject){
        var selectionContainer = $('<div></div>');

        var selectionNameContainer = $('<div></div>');
        var removeButton = new button(icons.minux, 20);
        var selectionName = $('<div></div>').text(ParsedSelectionObject['selection-name']);
        selectionName.css('display', 'inline-block');

        selectionNameContainer.append(removeButton);
        selectionNameContainer.append(selectionName);
        selectionContainer.append(selectionNameContainer);
        // Add event hander for selection remove

        var selectionProperties = $('<div></div>');
        // Add code for different types of input property
      }

      return selectionBox;
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
      var button = $('<div></div>');

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
      var formatted_content = $(svg).height(height-6);
      button.append(formatted_content);

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

      return button;
    }


    /**
     * position : Sets the css position property : absolute
     * @param {object} jquery html element
     * @param {number} left : css left property
     * @param {number} top : css top peroperty
     */
     function setPosition(ele, left, top){
       ele.css('position', 'absolute');
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
       var p_width = parent.width();
       var p_height = parent.height();

       // c_ stand for child
       var c_width = child.outerWidth();
       var c_height = child.outerHeight();

       var p_top = p_position.top + parseInt(parent.css('margin'));
       var p_left = p_position.left + parseInt(parent.css('margin')) ;

       var padding = parseInt(parent.css('padding').replace('px', ''));

      // Setting position
       var c_position = {
         left: (x_type == 'left' )? p_left + padding + x_offset : (x_type == 'right')? p_width - c_width - padding + x_offset : x_offset,
         top: (y_type == 'top' )? p_top + padding + y_offset : (y_type == 'bottom')? p_height - c_height - padding + y_offset : y_offset
       };

       console.log('setting parent child location', c_position, padding);
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
       var button = $('<div></div>');

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
       var formatted_content = $(svg).height(height-6);
       button.append(formatted_content);

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

       return button;
     }

  }

  return UI;
})();
