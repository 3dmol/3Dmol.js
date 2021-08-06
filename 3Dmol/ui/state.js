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
    var parentElement = glviewer.container;
    console.log('Container', parentElement);
    var height = canvas.height();
    var width = canvas.width();
    var offset = canvas.offset();

    var uiOverlayConfig = {
      height : height,
      width : width,
      offset : offset,
      ui : config.ui || undefined
    }
    // Selection Handlers
    var selections = [];
    var currentSelection = null;
    var currentForm = '';
    var currentStyles = null;
    var tempStyle = null;

    // Surface handlers
    var surfaces = {};

    // Label Handlers
    var labels = {};

    var atomLabel = {};

    this.setCurrentSelection = function(selectionId){
      currentSelection = selections.find(  sel => sel.id == selectionId);
      currentStyles = currentSelection.styles;
      this.ui.tools.selectionBox.styleBox.updateStyles(currentSelection.styles);
      // console.log('Updating Current Style', selections, currentSelection, currentStyles);
    }
    
    this.addSelection = function(){
      // console.log('Add Selection Called');
      var currentFormValue = { value: null, type: 'form', key: "Atom Selection"}
      var form =new $3Dmol.UI.Form($3Dmol.GLModel.validAtomSelectionSpecs, currentFormValue);
      this.ui.tools.dialog.addForm(form);
      currentForm = "new_selection";
    }

    this.addStyle = function(){
      // var formList = [];

      var validStyles = $3Dmol.GLModel.validAtomStyleSpecs;
      var currentFormValue = { value: null, type: 'form', key: "Atom Style"}
      // var forms = Object.keys(validStyles);
      // forms.forEach( (form)=>{
      //   let f = new $3Dmol.UI.Form(validStyles[form].validItems, form);
      //   formList.push({formName: form, form : f});
      // });
      let form = new $3Dmol.UI.Form(validStyles, currentFormValue);
      this.ui.tools.dialog.addForm(form);
      currentForm = "new_style";
    }

    // Style Handlers

    this.finalize = function(formSpecs, form){
      // console.log('Generated Config File', currentForm, formSpecs);
      if(currentForm == 'new_selection'){
        if( Object.keys(formSpecs).length != 0) {
          selectionSpec = { id : makeid(6), spec : formSpecs.value, hidden: false , styles : [], form : form, mouseInteraction: false, showResLabel: false};

          
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
        // currentSelection.form = form;
        this.ui.tools.selectionBox.editSelection(currentSelection);
        console.log("Edit_Selection", currentSelection.spec);
        render();
      }
      else if(currentForm == "edit_style"){
        var currentStyle = tempStyle;
        var index = currentSelection.styles.indexOf(currentSelection.styles.find((e)=>{
          console.log("Checking Style on Edit",e, currentStyle)
          if(e.id == currentStyle.id){
            return e;
          }
        }));

        currentSelection.styles[index].spec = formSpecs.value;
        render();
        
      }
    }

    this.removeStyle = function(styleSpec){
      var styleToRemove = currentSelection.styles.indexOf(styleSpec);
      if(styleToRemove != -1){
        currentSelection.styles.splice(styleToRemove,1);
        console.log('StateManager::removeStyle', styleToRemove, currentSelection.styles);
        this.ui.tools.selectionBox.styleBox.updateStyles(currentSelection.styles);
        render();
      }
    }

    // Context menu display test function 
    // $(document).on('click', { ui: this.ui.tools }, function(e){
    //   console.log('GLobal event', e);
    //   e.data.ui.contextMenu.show({ x: e.clientX, y:e.clientY });
    // });

    this.removeSelection = function(selectionSpec) {

      var selectionToRemove = selections.indexOf(selections.find((sel)=>{
        console.log("index of selection", sel);
        if(sel.id == selectionSpec.id)
          return sel
      }));

      if(selectionToRemove != -1){
        selections.splice(selectionToRemove, 1);
        console.log('StateManager::removeSelection', selectionToRemove, selectionSpec, selections);

        // glviewer.clear();
        clear(selectionSpec);
        render();
      }
    }
    
    this.cancel = function(){
      currentForm = "";
      console.log('Current Form Cancelled', currentForm);
    }

    this.showSelection = function(){
      var index = selections.indexOf(currentSelection);
      selections[index].hidden = false;
      // console.log('StateManger::showSelection', selections.map((s)=>{return s.hidden}), currentSelection);
      render();
    }

    this.hideSelection = function(){
      currentSelection.hidden = true;
      // console.log('StateManger::hideSelection', selections.map((s)=>{return s.hidden}), currentSelection);

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

    this.toggleMouseInteraction = function(){
      if(currentSelection.mouseInteraction){
        currentSelection.mouseInteraction = false;
        glviewer.setHoverable(currentSelection.spec, false);
        glviewer.setClickable(currentSelection.spec, false);
      }
      else{
        currentSelection.mouseInteraction = true;
        glviewer.setHoverable(currentSelection.spec, true, ()=>{
          console.log('Hovering Over Spheres', this);
        },
        ()=>{
          console.log('UnHovering Over Spheres', this);
        });

        glviewer.setClickable(currentSelection.spec, true, 
        (atom)=>{
          console.log('Clicking', atom, this);
        });

        
      }
      console.log("Toggle Mouse Interaction", currentSelection.spec, currentSelection.mouseInteraction);
      glviewer.render();
    }

    this.toggleResLabel = function(){
      if(currentSelection.showResLabel){
        currentSelection.showResLabel = false;
      }
      else {
        currentSelection.showResLabel = true;
      }

      // render();
    }

    this.editStyle = function(style){
      this.ui.tools.dialog.addForm(style.form);
      tempStyle = style;
      currentForm = "edit_style";
      console.log("StateManager::editStyle", currentSelection.form);
    }

    this.getSelectionList = function(){
      return selections.map((selection)=>{
        return { id: selection.id, spec: selection.spec};
      })
    }

    this.addSurface = function(property){
      var id = makeid(4);
      property.id = id;

      var style = property.surfaceStyle.value;
      if(style == null)
        style = {};

      var sel = selections.find((sel)=>{
        if(sel.id == property.surfaceOf.value)
          return true;
      });

      console.log('StateManager::addSurface', property, style);

      glviewer.addSurface(
        $3Dmol.SurfaceType[property.surfaceType.value],
        style,
        sel.spec
      ).then((surfId)=>{
        surfaces[id] = surfId;
        console.log("StateManager::Surfaces", surfaces);
      }, (err)=>{
        console.log('It failed', err);
      });

      return id;
    }

    this.removeSurface = function(id){
      
      glviewer.removeSurface(surfaces[id])
      delete surfaces[id];

    }

    this.editSurface = function(surfaceProperty){
      var style = surfaceProperty.surfaceStyle.value || {}
      var sel = selections.find((sel)=>{
        if(sel.id ==  surfaceProperty.surfaceOf.value)
          return true;
      })

      console.log('Surfaces edited', surfaceProperty.id, surfaces, surfaces[surfaceProperty.id]);
      glviewer.removeSurface(surfaces[surfaceProperty.id]);

      glviewer.addSurface(
        $3Dmol.SurfaceType[surfaceProperty.surfaceType.value],
        style,
        sel.spec
      ).then((surfId)=>{
        surfaces[surfaceProperty.id] = surfId;
      });

      console.log('StateManager::editSurface#Updating Surface', surfaceProperty)
    }

    this.openContextMenu = function(atom, x, y){
      console.log('Open Context Menu', atom, x, y);  
      var atomExist = false;

      if(atom){
        atomExist = Object.keys(atomLabel).find((i)=>{
          if (i == atom.index)
            return true;
          else 
            return false;
        });
  
        if(atomExist != undefined )
          atomExist = true;
        else 
          atomExist = false;
        
      }

      this.ui.tools.contextMenu.show(x, y, atom, atomExist);    
    }

    this.addLabel = function(labelValue){
      console.log('Label Added', labelValue);
      labels[labelValue.sel.value] = labels[labelValue.sel.value] || [];

      var labelProp = $3Dmol.labelStyles[labelValue.style.value];
      var selection = selections.find((sel)=>{
        if(sel.id == labelValue.sel.value)
          return true;
      });

      var offset = labels[labelValue.sel.value].length;
      labelProp['screenOffset'] = new $3Dmol.Vector2(0, -1*offset*35);

      labels[labelValue.sel.value].push(glviewer.addLabel(labelValue.text.value, labelProp, selection.spec));

      this.ui.tools.contextMenu.hide();
    }

    this.addAtomLabel = function(labelValue, atom){
      var atomExist = Object.keys(atomLabel).find((i)=>{
        if (i == atom.index)
          return true;
        else 
          return false;
      });

      if(atomExist != undefined )
        atomExist = true;
      else 
        atomExist = false;


      if(atomExist){
        this.removeAtomLabel(atom);
      }

      console.log('Add Atom Label Value', labelValue);
      
      atomLabel[atom.index] = atomLabel[atom.index] || null;
      
      var labelProp = $3Dmol.deepCopy($3Dmol.labelStyles['milk']);
      labelProp.position = {
        x : atom.x, y : atom.y, z : atom.z
      }

      var labelText = [];
      for ( key in labelValue){
        labelText.push(`${key} : ${labelValue[key]}`);
      }
      labelText = labelText.join('\n');

      atomLabel[atom.index] = glviewer.addLabel(labelText, labelProp);
      
      console.log('Getting Atom Label', labelText, labelProp);
    }

    this.exitContextMenu = function(){
      console.log('Unfinished Labeling');
      this.ui.tools.contextMenu.hide();
    }

    this.removeLabel = function(){
      // Add code to remove label 
      console.log('Remove Label')
      this.ui.tools.contextMenu.hide();
    }

    this.removeAtomLabel = function(atom){
      var label = atomLabel[atom.index];
      console.log("Stuff", label, atomLabel);
      glviewer.removeLabel(label);
      delete atomLabel[atom.index]; 
      
      console.log("Stuff After removal", label, atomLabel);
      this.ui.tools.contextMenu.hide();
    }

    canvas.on('click', ()=>{
      if(this.ui.tools.contextMenu.hidden == false){
        this.ui.tools.contextMenu.hide();
      }
    });

    
    // Setting up UI generation 
    this.showUI = function(){
      var ui = new $3Dmol.UI(this, uiOverlayConfig, parentElement);  

      return ui;
    }

    if(config.ui == true){
     this.ui = this.showUI();
    }
    // UI changes


    // console.log('GetSelectionList', this.getSelectionList());

    function render(){
      console.log('RENDER::', selections.map((s)=>{return s.hidden}));
      // glviewer.();
      glviewer.setStyle({});

      selections.forEach((sel)=>{
        if(!sel.hidden){
          // clear(sel);

          // var renderStyle = {};
          sel.styles.forEach((style)=>{
            if(!style.hidden){
              glviewer.addStyle(sel.spec, style.spec);
            }
          });

          glviewer.setClickable(sel.spec, true, ()=>{});
          glviewer.enableContextMenu(sel.spec, true);
        }
        else {
          glviewer.setClickable(sel.spec, false, ()=>{});
          glviewer.enableContextMenu(sel.spec, false);
        }

        console.log( "Context Menu Status", !sel.hidden);

      });

      glviewer.render();
    }

    function clear(selectionSpec){
      glviewer.enableContextMenu(sel.spec, false);
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
