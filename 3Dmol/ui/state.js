// @param {startState} object it store the initial values to initiliaze the
// UI state
// Everytime UI state is changed this state is updated so that this can
// preserve the changes in case error happens

// This will also help in migrating the controls from to different
$3Dmol.StateManager = (function(){

  function States(glviewer, config){
    var currentForm = '';


    console.log('State Manager Initiated For the viewer', glviewer);
    console.log('Viewer Configuration', );

    var canvas = $(glviewer.getRenderer().domElement);
    var height = canvas.height();
    var width = canvas.width();
    var offset = canvas.offset();

    var uiOverlayConfig = {
      height : height,
      width : width,
      offset : offset,
      ui : config.ui || undefined
    }

    var ui = new $3Dmol.UI(this, uiOverlayConfig);
    // UI changes
    $(window).on('resize', ()=>{
      var height = canvas.height();
      var width = canvas.width();
      var offset = canvas.offset();
      var config = {
        height : height,
        width : width,
        offset : offset,
      }
      ui.resize(config);
      console.log('Resizing');
    });

    // Selection Handlers
    var selections = [];
    this.addSelection = function(){
      // console.log('Add Selection Called');

      var form =new $3Dmol.UI.form($3Dmol.GLModel.validAtomSelectionSpecs, "Atom Selection");
      ui.tools.dialog.addForm(form);
      currentForm = "Selection";
    }

    // Style Handlers
    var styles = [];

    this.finalize = function(formSpecs){
      console.log('Generated Config File', currentForm, formSpecs);
      if( Object.keys(formSpecs).length != 0) {
        spec = { id : makeid(6), spec : formSpecs };
        ui.tools.selectionBox.appendSelection(spec);
      }
    }

    this.cancel = function(){
      console.log('Current Form Cancelled', currentForm);
    }

    this.addStyle = function(){
      console.log('Add Style Called');
    }

    function makeid(length) {
      var result           = '';
      var characters       = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
      var charactersLength = characters.length;
      for ( var i = 0; i < length; i++ ) {
        result += characters.charAt(Math.floor(Math.random() * charactersLength));
     }
     return result;
    }
  }

  return States;
})()
