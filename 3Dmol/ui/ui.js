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


    // Generates the necessary UI elements
    var icons =new $3Dmol.UI.Icons();
    addToViewer();

    function addToViewer(){
      var topbar =new Toolbar();
      $('body').append(topbar);
      var selectionBox = new SelectionBox();
      $('body').append(selectionBox);
    }

    /**
    * function to make the html toolbar div
    * @param x position
    * @param y position
    * @param width
    * @param height
    **/
    function Toolbar(){
      var toolbar = $('<div></div>');

      // Toolbar content

      var optionButton = new button_toolbar(icons.option, 20)
      var moveButton = new button_toolbar(icons.move, 20);
      var rotateButton = new button_toolbar(icons.rotate, 20);
      toolbar.append(optionButton);
      toolbar.append(moveButton);
      toolbar.append(rotateButton);

      // Toolbar design
      toolbar.css('position', 'absolute');
      toolbar.css('border', '1px solid black');
      toolbar.css('padding', '3px');
      toolbar.css('background', 'black');
      toolbar.css('color', 'white');
      toolbar.css('width', '300px');

      return toolbar;
    }

    /**
     * button - generates button with the given markup as contents
     * @param {string} SVG markup string that contains the content of the button
     * @param {number} Height of the content
     * @param {}
     * @return {object}  Returns the jquery object for button
     */
    function button_toolbar(svg, height){
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

      return button;
    }


    /**
     * SelectionBox - Draws the box where all the selections on atoms are listed
     * This will be used to modify style for specific set of selection
     *
     * @return {Object}  Jquery element of div
     */
    function SelectionBox() {
      var selectionBox = $('<div></div>')

      // Content

      // CSS
      selectionBox.css('box-sizing', 'border-box');
      selectionBox.css('padding', '3px');
      selectionBox.css('width', '120px');
      selectionBox.css('height', '300px');
      selectionBox.css('border-radius', '6px');
      selectionBox.css('background', 'rgb(214, 214, 214)');

      // CSS Position
      selectionBox.css('position','absolute');
      // selectionBox.css('top', '0px');
      // selectionBox.css('left', '0px');

      // Selection Box modification function

      /**
       * Add a selection input to the list
       */
      this.addSelection = function(ParsedSelectionObject){

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



  }

  return UI;
})();
