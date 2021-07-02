// @param {startState} object it store the initial values to initiliaze the
// UI state
// Everytime UI state is changed this state is updated so that this can
// preserve the changes in case error happens

// This will also help in migrating the controls from to different
$3Dmol.StateManager = (function(){

  function States(glviewer, config){
    config = config || {};
    config.ui = config.ui || false;
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

    this.showUI = function(){
      this.ui = new $3Dmol.UI(this, uiOverlayConfig);  
      $(window).on('resize', ()=>{
        var height = canvas.height();
        var width = canvas.width();
        var offset = canvas.offset();
        var config = {
          height : height,
          width : width,
          offset : offset,
        }
        this.ui.resize(config);
        console.log('Resizing');
      });

    }

    if(config.ui == true){
     this.ui = new $3Dmol.UI(this, uiOverlayConfig);
    }
    // UI changes

    // Selection Handlers
    var selections = [];
    var currentSelection = null;
    var currentForm = '';
    var currentStyles = null;

    this.setCurrentSelection = function(selectionId){
      currentSelection = selections.find(  sel => sel.id == selectionId);
      currentStyles = currentSelection.styles;
      this.ui.tools.selectionBox.styleBox.updateStyles(currentSelection.styles);
      console.log('Updating Current Style', selections, currentSelection, currentStyles);
    }
    
    this.addSelection = function(){
      // console.log('Add Selection Called');
      var form =new $3Dmol.UI.Form($3Dmol.GLModel.validAtomSelectionSpecs, "Atom Selection");
      this.ui.tools.dialog.addForm(form);
      currentForm = "new_selection";
    }

    this.addStyle = function(){
      var formList = [];

      var validStyles = $3Dmol.GLModel.validAtomStyleSpecs;
      var forms = Object.keys(validStyles);
      forms.forEach( (form)=>{
        let f = new $3Dmol.UI.Form(validStyles[form].validItems, form);
        formList.push({formName: form, form : f});
      });


      this.ui.tools.dialog.addForm(formList, 'list');
      currentForm = "new_style";
    }

    // Style Handlers

    this.finalize = function(formSpecs){
      console.log('Generated Config File', currentForm, formSpecs);
      if(currentForm == 'new_selection'){
        if( Object.keys(formSpecs).length != 0) {
          selectionSpec = { id : makeid(6), spec : formSpecs , styles : []};
          this.ui.tools.selectionBox.appendSelection(selectionSpec);
          selections.push(selectionSpec);
          this.setCurrentSelection(selectionSpec.id);

          // console.log(selectionSpec, selections);
        }
      }
      if( currentForm == "new_style"){
        if( Object.keys(formSpecs).length != 0) {
          if( currentSelection != null ){
            styleSpec = { id : makeid(6), spec : formSpecs , styles : []};
            currentSelection.styles.push(styleSpec);
            this.ui.tools.selectionBox.styleBox.updateStyles(currentSelection.styles);
            // this.setCurrentSelection(selectionSpec.id);
          }
          console.log(currentSelection);
        }
      }
    }

    
    this.cancel = function(){
      console.log('Current Form Cancelled', currentForm);
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
