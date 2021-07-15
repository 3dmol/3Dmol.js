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

    this.currentFormValue = null;

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
    var tempStyle = null;

    this.setCurrentSelection = function(selectionId){
      currentSelection = selections.find(  sel => sel.id == selectionId);
      currentStyles = currentSelection.styles;
      this.ui.tools.selectionBox.styleBox.updateStyles(currentSelection.styles);
      // console.log('Updating Current Style', selections, currentSelection, currentStyles);
    }
    
    this.addSelection = function(){
      // console.log('Add Selection Called');
      this.currentFormValue = { value: null, type: 'form', key: "Atom Selection"}
      var form =new $3Dmol.UI.Form($3Dmol.GLModel.validAtomSelectionSpecs, this.currentFormValue);
      this.ui.tools.dialog.addForm(form);
      currentForm = "new_selection";
    }

    this.addStyle = function(){
      // var formList = [];

      var validStyles = $3Dmol.GLModel.validAtomStyleSpecs;
      this.currentFormValue = { value: null, type: 'form', key: "Atom Style"}
      // var forms = Object.keys(validStyles);
      // forms.forEach( (form)=>{
      //   let f = new $3Dmol.UI.Form(validStyles[form].validItems, form);
      //   formList.push({formName: form, form : f});
      // });
      let form = new $3Dmol.UI.Form(validStyles, this.currentFormValue);
      this.ui.tools.dialog.addForm(form);
      currentForm = "new_style";
    }

    // Style Handlers

    this.finalize = function(formSpecs, form){
      // console.log('Generated Config File', currentForm, formSpecs);
      if(currentForm == 'new_selection'){
        if( Object.keys(formSpecs).length != 0) {
          selectionSpec = { id : makeid(6), spec : formSpecs.value, hidden: false , styles : [], form : form};

          
          // Value Sanitization
          if(selectionSpec.spec.resi){
            selectionSpec.spec.resi = selectionSpec.spec.resi.split(",");
          }
          console.log('Form:Specs:Value', selectionSpec);
          
          this.ui.tools.selectionBox.appendSelection(selectionSpec);
          selections.push(selectionSpec);

          this.setCurrentSelection(selectionSpec.id);

          // console.log(selectionSpec, selections);
        }
      }
      else if( currentForm == "new_style"){
        if( Object.keys(formSpecs).length != 0) {
          if( currentSelection != null ){
            styleSpec = { id : makeid(6), spec : formSpecs.value, hidden : false , styles : [], form : form};
            currentSelection.styles.push(styleSpec);
            this.ui.tools.selectionBox.styleBox.updateStyles(currentSelection.styles);

            console.log(currentSelection.spec, styleSpec.spec);

            render();
            // this.setCurrentSelection(selectionSpec.id);
          }
          // console.log(currentSelection);
        }
      }

      else if(currentForm == "edit_selection") {
        currentSelection.spec = formSpecs.value;
        currentSelection.form = form;
        render();
      }
      else if(currentForm == "edit_style"){
        var currentStyle = tempStyle;
        var index = currentSelection.styles.indexOf(currentSelection.styles.find((e)=>{
          if(e == currentStyle){
            return e;
          }
        }));

        currentSelection.styles[index].spec = formSpecs.spec;
        
      }
    }

    this.removeStyle = function(styleSpec){
      var styleToRemove = currentSelection.styles.indexOf(styleSpec);
      currentSelection.styles.splice(styleToRemove,1);
      console.log('StateManager::removeStyle', styleToRemove, currentSelection.styles);

      render();
    }

    this.removeSelection = function(selectionSpec) {
      // glviewer.clear();

      var selectionToRemove = selections.indexOf(selections.find((sel)=>{
        console.log("index of selection", sel);
        if(sel.id == selectionSpec.id)
          return sel
        }));
      // glviewer.setStyle(currentSelection)

      selections.splice(selectionToRemove, 1);
      console.log('StateManager::removeSelection', selectionToRemove, selectionSpec, selections);

      // glviewer.clear();
      clear(selectionSpec);
      render();
    }
    
    this.cancel = function(){
      console.log('Current Form Cancelled', currentForm);
    }

    this.showSelection = function(){
      var index = selections.indexOf(currentSelection);
      selections[index].hidden = false;
      console.log('StateManger::showSelection', selections.map((s)=>{return s.hidden}), currentSelection);
      render();
    }

    this.hideSelection = function(){
      currentSelection.hidden = true;
      console.log('StateManger::hideSelection', selections.map((s)=>{return s.hidden}), currentSelection);
      render();
    }

    this.hideStyle = function(styleSpec){
      var index = currentSelection.styles.indexOf(styleSpec);
      currentSelection.styles[index].hidden = true;

      render();
    }

    this.showStyle = function(styleSpec){
      var index = currentSelection.styles.indexOf(styleSpec);
      currentSelection.styles[index].hidden = false;

      render();
    }

    this.editSelection = function(){
      this.ui.tools.dialog.addForm(currentSelection.form);
      currentForm = "edit_selection";
      console.log("StateManager::editSelection", currentSelection.form);
    }

    this.editStyle = function(style){
      this.ui.tools.dialog.addForm(style.form);
      tempStyle = style;
      currentForm = "edit_style";
      console.log("StateManager::editStyle", currentSelection.form);
    }

    function render(){
      console.log('RENDER::', selections.map((s)=>{return s.hidden}));
      // glviewer.();

      selections.forEach((sel)=>{
        if(!sel.hidden){
          clear(sel);

          // var renderStyle = {};
          sel.styles.forEach((style)=>{
            if(!style.hidden){
              glviewer.addStyle(sel.spec, style.spec);
            }
          })
        }
        else{
          glviewer.setStyle(sel.spec, {});
        }
      });

      glviewer.render();
    }

    function clear(selectionSpec){
      glviewer.setStyle(selectionSpec.spec, {});
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
