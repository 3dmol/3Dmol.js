// @ts-ignore
import $ from 'jquery';
import SurfaceType from '../enum/SurfaceType';
import deepCopy from '../util/deepCopy';
import download from '../util/download';
import { Vector2 } from '../WebGL/math';
import { labelStyles } from './defaultValues';
import UI from './ui';
/**
 * StateManager - StateManager creates the space to preserve the state of the ui and sync it with the GLViewer
 * @constructor 
 * @param {import("../GLViewer").default} glviewer StateManager is required to have interaction between glviewer and the ui. 
 * @param {Object} config Loads the user defined parameters to generate the ui and handle state
 */
const StateManager = (function(){

  function States(glviewer, config){
    config = config || {};
    config.ui = config.ui || false;

    const canvas = $(glviewer.getRenderer().domElement);
    const parentElement = glviewer.container;

    const height = parentElement.height();
    const width = parentElement.width();
    const offset = canvas.offset();
    // @ts-ignore
    const stateManager = this;

    const uiOverlayConfig = {
      height,
      width,
      offset,
      ui : config.ui || undefined
    }

    // Selection Handlers
    let selections = {};

    // Surface handlers
    let surfaces = {};

    // Label Handlers
    let labels = {};

    let atomLabel = {};

    /**
     * Add Selection from the ui to glviewer
     * 
     * @function StateManager#addSelection
     * @param {import("../specs").AtomSelectionSpec} spec Object that contains the output from the form 
     * @param {string | null} sid If surface id being edited then sid is set to some string
     * @returns {string | null}
     */
    this.addSelection = function(spec, sid = null){

      const id = sid || makeid(4);

      const selectionSpec = {
        spec,
        styles : {},
        hidden : false
      };

      if(sid === null)
        selections[id] = selectionSpec;
      else 
        selections[id].spec = selectionSpec.spec;

      render();
      return id;
    }

    /**
     * Return true if the selections contain at least one atom
     * 
     * @function StateManager#checkAtoms
     * @param {import('../specs').AtomSelectionSpec} sel Atom selection spec
     * @returns Boolean
     */
    this.checkAtoms = function(sel){
      const atoms = glviewer.selectedAtoms(sel);
      if( atoms.length > 0)
        return true

      return false;
    }

    /**
     * Toggle the hidden property of the selection 
     * @function StateManager#toggleHide
     * @param {String} sid Selection id
     */
    this.toggleHide = function(sid){
      selections[sid].hidden = !selections[sid].hidden;
      render();
    }

    /**
     * Removes the selection
     * @param {String} id Selection id
     */
    this.removeSelection = function(id) {
      delete selections[id];
      render();
    }

    /**
     * Add style and renders it into the viewport
     * 
     * @function StateManager#addStyle 
     * @param {import('../specs').AtomStyleSpec} spec Output object of style form 
     * @param {string | null} sid Selection Id
     * @param {string | null} stid Style Id
     * @returns {string | null}
     */
    this.addStyle = function( spec, sid, stid = null){
      const selection = selections[sid];
      
      
      const styleSpec = {
        spec,
        hidden : false
      }
      
      let id = null; 
      
      if(stid === null) {
        id = makeid(4);
        selection.styles[id] = styleSpec
      }
      else {
        id = stid;
        selection.styles[id].spec = spec;
      }
      
      render();

      return id;
    }

    
    /**
     * Removes the style specified by stid
     * 
     * @function StateManager#removeStyle 
     * @param {String} sid Selection id
     * @param {String} stid Style Id
     */
    this.removeStyle = function(sid, stid){
      delete selections[sid].styles[stid];
      render();
    }


    /**
     * Toggle hidden property of a style 
     * 
     * @function StateManager#toggleHideStyle
     * @param {String} sid Selection Id
     * @param {String} stid Style Id 
     */
    this.toggleHideStyle = function(sid, stid){
      selections[sid].styles[stid].hidden = !selections[sid].styles[stid].hidden;
      render();
    }

    /**
     * Adds surface to the viewport
     * 
     * @function StateManager#addSurface
     * @param {Object} property Surface output object
     * @param {Function} callback callback
     * @returns String
     */
    this.addSurface = function(property, callback){
      const id = makeid(4);
      property.id = id;

      let style = property.surfaceStyle.value;
      if(style === null)
        style = {};

      const sel = (property.surfaceFor.value === 'all') ? { spec : {} } : selections[property.surfaceFor.value];

      const generatorAtom = (property.surfaceOf.value === 'self')? sel.spec : {};


      glviewer.addSurface(
        SurfaceType[property.surfaceType.value],
        style,
        sel.spec,
        generatorAtom
      ).then((surfParam)=>{
        surfaces[id] = surfParam[0];

        if(callback !== undefined)
          callback(id, surfParam[0]);
      // @ts-ignore
      }, (err)=>{

      });

      return id;
    }

    /**
     * Removes surface from the viewport 
     * @function StateManager#removeSurface
     * @param {String} id Surface Id
     */
    this.removeSurface = function(id){
      glviewer.removeSurface(surfaces[id])

      delete surfaces[id];

    }

    /**
     * Edit the exisiting surface in the viewport
     * 
     * @function StateManager#editSurface
     * @param {Object} surfaceProperty Surface Style
     */
    this.editSurface = function(surfaceProperty){
      const style = surfaceProperty.surfaceStyle.value || {}

      const sel = (surfaceProperty.surfaceFor.value === 'all') ? { spec : {} } : selections[surfaceProperty.surfaceFor.value];
      const generatorAtom = (surfaceProperty.surfaceOf.value === 'self')? sel.spec : {};

      glviewer.removeSurface(surfaces[surfaceProperty.id]);

      console.log(surfaceProperty);
      glviewer.addSurface(
        SurfaceType[surfaceProperty.surfaceType.value],
        style,
        sel.spec,
        generatorAtom
      ).then((surfId)=>{
        surfaces[surfaceProperty.id] = surfId[0];
      });
    }

    /**
     * Returns the list of ids of selections that are created so far
     * @function StateManager#getSelectionList
     * @returns <Array of selection ids>
     */
    this.getSelectionList = function(){
      return Object.keys(selections);
    }

    /**
     * Opens context menu when called from glviewer
     * 
     * @function StateManager#openContextMenu
     * @param {import('../specs').AtomSpec} atom Atom spec obtained from context menu event
     * @param {Number} x x coordinate of mouse on viewport
     * @param {Number} y y coordinate of mouse on the viewport
     */
    this.openContextMenu = function(atom, x, y){ 
      let atomExist = false;

      if(atom){
        atomExist = !!Object.keys(atomLabel).find((i)=>{
          // @ts-ignore ignore the type error and eqeqeq rule so that we can compare str and number
          // eslint-disable-next-line eqeqeq
          if (i == atom.index)
            return true;
          return false;
        });
  
        if(atomExist !== undefined )
          atomExist = true;
        else 
          atomExist = false;
        
      }

      if(this.ui) this.ui.tools.contextMenu.show(x, y, atom, atomExist);    
    }

    /**
     * Adds Label to the viewport specific to the selection
     * @function StateManager#addLabel
     * @param {Object} labelValue Output object from label form of Context Menu
     */
    this.addLabel = function(labelValue){
      labels[labelValue.sel.value] = labels[labelValue.sel.value] || [];

      const labelProp = labelStyles[labelValue.style.value];
      const selection = selections[labelValue.sel.value];

      const offset = labels[labelValue.sel.value].length;
      labelProp.screenOffset = new Vector2(0, -1*offset*35);

      labels[labelValue.sel.value].push(glviewer.addLabel(labelValue.text.value, labelProp, selection.spec));

      // @ts-ignore
      this.ui.tools.contextMenu.hide();
    }

    /**
     * Adds atom label to the viewport
     * 
     * @function StateManager#addAtomLabel
     * @param {Object} labelValue Output object from propertyMenu form of Context Menu
     * @param {import('../specs').AtomSpec} atom Atom spec that are to be added in the label 
     */
    this.addAtomLabel = function(labelValue, atom, styleName='milk'){
      let atomExist = !!Object.keys(atomLabel).find((i)=>{
        // @ts-ignore ignore type checks and eqeqeq rule so that we can compare str and number
        // eslint-disable-next-line eqeqeq
        if (i == atom.index)
          return true;
        return false;
      });

      if(atomExist !== undefined )
        atomExist = true;
      else 
        atomExist = false;


      if(atomExist){
        this.removeAtomLabel(atom);
      }

      
      atomLabel[atom.index] = atomLabel[atom.index] || null;
      
      const labelProp = deepCopy(labelStyles[styleName]);
      labelProp.position = {
        x : atom.x, y : atom.y, z : atom.z
      }

      /** @type {string | Array<string>} */
      let labelText = [];
      for (let key in labelValue){
        labelText.push(`${key} : ${labelValue[key]}`);
      }
      labelText = labelText.join('\n');

      atomLabel[atom.index] = glviewer.addLabel(labelText, labelProp);
      
    }

    /**
     * Executes hide context menu and process the label if needed
     * 
     * @function StateManager#exitContextMenu
     * @param {Boolean} processContextMenu Specify the need to process the values in the context menu
     */
    this.exitContextMenu = function(processContextMenu = false){
        if(this.ui) {
            this.ui.tools.contextMenu.hide(processContextMenu);
        }
    }

    /**
     * Removes the label specific to the selection 
     * 
     * (under development)
     */
    this.removeLabel = function(){
      // Add code to remove label 
      // @ts-ignore
      this.ui.tools.contextMenu.hide();
    }

    /**
     * Removes the atom label from the viewpoer 
     * @function StateManager#removeAtomLabel
     * @param {import('../specs').AtomSpec} atom Atom spec
     */
    this.removeAtomLabel = function(atom){
      const label = atomLabel[atom.index];
      glviewer.removeLabel(label);
      delete atomLabel[atom.index]; 
      
      // @ts-ignore
      this.ui.tools.contextMenu.hide();
    }

    /**
     * Add model to the viewport
     * @function StateManager#addModel
     * @param {Object} modelDesc Model Toolbar output
     */
    this.addModel = function(modelDesc){
      glviewer.removeAllModels();
      glviewer.removeAllSurfaces();
      glviewer.removeAllLabels();
      glviewer.removeAllShapes();

      const query = `${modelDesc.urlType.value  }:${  modelDesc.url.value}`;
      download(query, glviewer, {}, ()=>{
        // @ts-ignore
        this.ui.tools.modelToolBar.setModel(modelDesc.url.value.toUpperCase());
      });

      // Remove all Selections
      selections = {};
      surfaces = {};
      atomLabel = {};
      labels = {};

      // Reset UI
      // @ts-ignore
      this.ui.tools.selectionBox.empty();
      // @ts-ignore
      this.ui.tools.surfaceMenu.empty();
    }

    // State Management helper function 
    function findSelectionBySpec(spec){
      const ids = Object.keys(selections);
      let matchingObjectIds = null;
      for(let i = 0; i < ids.length; i++){
        const lookSelection = selections[ids[i]].spec;

        let match = true;
        
        // looking for same parameters length 
        const parameters = Object.keys(spec);

        if( Object.keys(lookSelection).length === parameters.length){
          for(let j = 0; j < parameters.length; j++){
            if( lookSelection[parameters[j]] !== spec[parameters[j]]){
              match = false;
              break;
            }
          }
        } else {
          match = false;
        }

        if(match){
          matchingObjectIds = ids[i];
          break;
        }
      }

      return matchingObjectIds;
    }

    // State managment function 

    /**
     * Updates the state variable for selections and styles and trigger ui to show the 
     * ui elements for these selections and styles.
     * 
     * @function StateManager#createSelectionAndStyle
     * @param {import('../specs').AtomSelectionSpec} selSpec Atom Selection Spec
     * @param {import('../specs').AtomStyleSpec} styleSpec Atom Style Spec
     */
    this.createSelectionAndStyle = function(selSpec, styleSpec){

      let selId = findSelectionBySpec(selSpec);

      if(selId === null){
        selId = this.addSelection(selSpec);
      }

      let styleId = null;

      if(Object.keys(styleSpec).length !== 0){
        styleId = this.addStyle(styleSpec, selId);
      }

      // @ts-ignore
      this.ui.tools.selectionBox.editSelection(selId, selSpec, styleId, styleSpec);
      
    };

    /**
     * Creates selection and add surface with reference to that selection 
     * and triggers updates in the ui
     * @function StateManager#createSurface
     * @param {string} surfaceType Type of surface to be created
     * @param {import('../specs').AtomSelectionSpec} sel Atom selection spec
     * @param {import('../specs').AtomStyleSpec} style Atom style spec
     * @param {string|null} sid selection id
     */
    this.createSurface = function(surfaceType, sel, style, sid){
      let selId = findSelectionBySpec(sel);
      
      if(selId === null){
        selId = this.addSelection(sel);

      }
      // @ts-ignore
      this.ui.tools.selectionBox.editSelection(selId, sel, null);

      surfaceType = Object.keys(style)[0];

      const surfaceInput = {
        surfaceType : {
          value : surfaceType
        },

        surfaceStyle : {
          value : style[surfaceType],
        },

        surfaceOf : {
          value : 'self'
        },

        surfaceFor : {
          value : selId
        }
      }

      const surfId = makeid(4);
      surfaces[surfId] = sid;

      // @ts-ignore
      this.ui.tools.surfaceMenu.addSurface(surfId, surfaceInput);

      // Create Surface UI
    };

    /**
     * Sets the value of title in ModelToolBar
     * @function StateManager#setModelTitle
     * @param {String} title Model title
     */
    this.setModelTitle = function(title){
      // @ts-ignore
      this.ui.tools.modelToolBar.setModel(title);
    }

    canvas.on('click', ()=>{
      if(this.ui && this.ui.tools.contextMenu.hidden === false){
        this.ui.tools.contextMenu.hide();
      }
    });
    
    // Setting up UI generation 
    /**
     * Generates the ui and returns its reference
     * @returns UI
     */
    this.showUI = function(){
      const ui = new UI(this, uiOverlayConfig, parentElement);  
      return ui;
    };

    if(config.ui === true){
     this.ui = this.showUI(); 
    };

    this.initiateUI = function(){
      this.ui = new UI(this, uiOverlayConfig, parentElement);
      render();
    }
    /**
     * Updates the UI on viewport change 
     * 
     * @function StateManager#updateUI
     */
    this.updateUI = function(){
      if(this.ui){
        this.ui.resize();
      }
    };
    
    // UI changes

    function render(){
      // glviewer.();
      glviewer.setStyle({});

      const selList = Object.keys(selections);

      selList.forEach( (selKey) =>{
        const sel = selections[selKey];

        if( !sel.hidden ) {
          const styleList = Object.keys(sel.styles);
          
          styleList.forEach((styleKey)=>{
            const style = sel.styles[styleKey];

            if( !style.hidden){
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

      })

      glviewer.render();
    }

    function makeid(length) {
      let result           = '';
      const characters       = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
      const charactersLength = characters.length;
      for ( let i = 0; i < length; i++ ) {
        result += characters.charAt(Math.floor(Math.random() * charactersLength));
     }
     return result;
    }
  }

  return States;
})()

export default StateManager