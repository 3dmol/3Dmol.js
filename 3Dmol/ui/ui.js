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

    addToViewer();

    function addToViewer(){
      var topbar = toolbar();
      $('body').append(topbar);

      var icons =new $3Dmol.UI.Icons();
      $('body').append(icons.rotate);
    }
    // High level elements
    function preferences(){

    }


    /**
    * function to make the html toolbar div
    * @param x position
    * @param y position
    * @param width
    * @param height
    **/
    function toolbar(){
      var toolbar = $('<div></div>');

      // Toolbar content
      toolbar.text('This is cool');

      // Toolbar design
      toolbar.css('position', 'absolute');
      toolbar.css('border', '1px solid black');
      toolbar.css('padding', '3px');
      toolbar.css('background', 'black');
      toolbar.css('color', 'white');
      toolbar.css('width', '300px');



      return toolbar;
    }

    function toolbar


  }

  return UI;
})();
