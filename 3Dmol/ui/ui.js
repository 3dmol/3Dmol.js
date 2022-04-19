import $ from 'jquery';
import SurfaceType from '../enum/SurfaceType';
import GLModel from '../GLModel';
import editingForm from '../singletons/editingForm';
import {labelStyles} from './defaultValues';
import Form from './Form';
import Icons from './Icons';

/**
 * UI - UI creates panels in the viewer to assist control of the viewport
 * @constructor
 * @param {import("./StateManager").default} stateManager StateManager is required to have interaction between glviewer and the ui.
 * @param {Object} config Loads the user defined parameters to generate the ui
 * @param {Object} parentElement Refers the parent division used to hold the canvas for 3Dmol.js
 */
const UI = (function () {
  function UI(stateManager, config, parentElement) {
    config = config || {};

    // Extract the viewer and then render it
    const icons = new Icons();
    const body = $('body');

    const mainParent = $(parentElement[0]);
    // Generates the necessary UI elements
    const HEIGHT = config.height;
    const WIDTH = config.width;
    const uiElements = (this.tools = generateUI(config));
    /**
     * Resize the panel with respect to the new viewport
     *
     * @function UI#resize
     */
    this.resize = function () {
      const selectionBox = this.tools.selectionBox;
      const surfaceMenu = this.tools.surfaceMenu;
      const modelToolBar = this.tools.modelToolBar;
      const HEIGHT = mainParent.height();

      setLocation(mainParent, modelToolBar.ui, 'left', 'top');
      // modelToolBar.updateInputLength();
      setLocation(mainParent, selectionBox.ui, 'left', 'top', 0, modelToolBar.ui.height() + 5);
      selectionBox.updateScrollBox(HEIGHT);
      setLocation(mainParent, surfaceMenu.ui, 'right', 'top', 0, modelToolBar.ui.height() + 5);
      surfaceMenu.updateScrollBox(HEIGHT);
    };

    /**
     * button - generates button with the given markup as contents
     * @param {String} svg SVG markup string that contains the content of the button
     * @param {Number} height Height of the content
     * @param {Object} [config] Various properties to define the button
     *
     */
    function button(svg, height, config) {
      config = config || {};
      const borderRadius = config.bfr * height || height / 4; // body radius factor
      const bgColor = config.backgroundColor || 'rgb(177, 194, 203)';
      const color = config.color || 'black';
      const hoverable = config.hoverable || 'true';
      const tooltipText = config.tooltip || null;

      // Button instance
      this.ui = $('<div></div>');
      const innerButton = $('<div></div>');
      this.ui.append(innerButton);

      // CSS
      this.ui.css('box-sizing', 'border-box');
      this.ui.css('display', 'inline-block');
      this.ui.css('margin', '3px');
      this.ui.css('height', height);
      this.ui.css('width', height);
      this.ui.css('border-radius', `${borderRadius}px`);

      //  button.css('padding', '3px');
      this.ui.css('color', color);
      this.ui.css('background', bgColor);

      innerButton.css('display', 'flex');
      innerButton.css('justify-content', 'center');
      innerButton.css('align-items', 'center');
      innerButton.css('padding', '2px');

      // content
      this.setSVG = function (svg) {
        innerButton.empty();
        const formatted_content = $(svg);
        innerButton.append(formatted_content);
      };

      this.setSVG(svg);

      // Hover

      // Setting up tool tip
      this.ui.css({
        position: 'relative',
      });

      // setting up tool tip
      if (tooltipText !== null) {
        this.ui.attr('title', tooltipText);
      }

      if (hoverable === 'true') {
        this.ui
          .on('mouseenter', () => {
            this.ui.css('box-shadow', '0px 0px 3px black');
          })
          .on('mouseleave', () => {
            this.ui.css('box-shadow', 'none');
          });

        const longPressTime = 0;
        const mouseX = 0;
        const mouseY = 0;

        // click
        this.ui.on('mousedown', e => {
          this.ui.css('box-shadow', '0px 0px 1px black');
        });

        this.ui.on('mouseup', () => {
          this.ui.css('box-shadow', '0px 0px 3px black');
        });

        this.ui.on('mousemove', e => {
          // mouseX = e.clientX;
          // mouseY = e.clientY;
        });
      }
    }

    /**
     * ModelToolbar is part of UI to edit or change the model loaded into the viewer
     *
     * @function ModelToolbar
     */
    function ModelToolbar() {
      const boundingBox = (this.ui = $('<div></div>'));

      boundingBox.css({
        position: 'relative',
        'min-width': '150px',
      });

      const modelButton = new button(Icons.molecule, 20, {tooltip: 'Toggle Model Selection Bar'});
      boundingBox.append(modelButton.ui);

      modelButton.ui.css({
        display: 'inline-block',
        top: '3px',
      });

      const control = {
        urlType: {
          active: true,
          value: null,
          key: 'Model type',
        },

        url: {
          active: true,
          value: null,
          key: 'Url',
        },
      };

      const surroundingBox = $('<div></div>');

      surroundingBox.css({
        display: 'inline-block',
        background: '#e4e4e4',
        padding: '2px',
        'border-radius': '3px',
        // 'width' : '90%'
      });

      boundingBox.append(surroundingBox);

      const currentModelBox = $('<div></div>');
      currentModelBox.css({});

      const currentModel = $('<div></div>');
      currentModel.css({
        display: 'inline-block',
        'font-family': 'Arial',
        'font-size': '12px',
        'font-weight': 'bold',
        // 'padding' : '3px'
      });

      currentModelBox.append(currentModel);

      const changeButton = new button(Icons.change, 16, {
        tooltip: 'Change Model',
        backgroundColor: 'white',
        bfr: 0.5,
      });
      changeButton.ui.css({
        display: 'inline-block',
        'margin-left': '4px',
      });
      currentModelBox.append(changeButton.ui);

      currentModelBox.hide();
      surroundingBox.append(currentModelBox);

      const formBox = $('<div></div>');
      surroundingBox.append(formBox);

      const dbs = 'pdb,mmtf,cid'.split(',');
      const list = (this.list = new Form.ListInput(control.urlType, dbs));
      list.showAlertBox = false;

      list.ui.css({
        display: 'inline-block',
      });

      formBox.append(list.ui);

      const input = (this.url = new Form.Input(control.url));
      formBox.append(input.ui);

      input.ui.css({
        display: 'inline-block',
        width: '125px',
      });

      // input.setWidth(125);

      const submitButton = new button(Icons.tick, 16, {
        bfr: 0.5,
        backgroundColor: 'lightgreen',
        tooltip: 'Add Model',
      });
      submitButton.ui.css({
        margin: '0px',
      });
      formBox.append(submitButton.ui);

      this.updateInputLength = function () {
        // var width = parentElement.width()*0.3;
        // boundingBox.width(width);
        // input.setWidth(width - 12);
      };

      modelButton.ui.on('click', () => {
        surroundingBox.toggle();
      });

      submitButton.ui.on('click', () => {
        const validateDb = list.validate();
        const validateId = input.validate();

        if (validateId && validateDb) {
          stateManager.addModel(control);
        }
      });

      /**
       * Sets the title in the ui with specified value
       *
       * @function ModelToolbar#setModel
       * @param {String} heading Name of the molecule that is to be displayed on the title
       */
      this.setModel = function (heading) {
        currentModel.text(heading);
        currentModelBox.show();
        formBox.hide();
      };

      changeButton.ui.on('click', () => {
        currentModelBox.hide();
        formBox.show();
        input.setValue('');
      });

      boundingBox.on('keypress', e => {
        if (e.key === 'Enter' || e.key === 'Return') {
          submitButton.ui.trigger('click');
        }
      });
    }

    /**
     * Selection box creates the UI panel to manipulate selections and style that are drawn
     * on the viewport
     *
     * @function SelectionBox
     * @param {string} icon takes the svg code for the icon that is to be used to display
     * selection box
     * @return {Object}  Jquery element of div
     */
    function SelectionBox(icon, side = 'left') {
      const selectionBox = (this.ui = $('<div></div>'));
      editingForm.isEditing = false;
      const selectionObjects = [];

      const selections = $('<div></div>');
      const scrollBox = $('<div></div>');

      selections.css('opacity', '0.9');

      const showArea = $('<div></div>');
      const addArea = $('<div></div>');
      const plusButton = new button(Icons.plus, 20, {tooltip: 'Add New Selection'});
      plusButton.ui.css('margin', '0px');

      const hideButton = new button(icon, 20, {tooltip: 'Toggle Selection Menu'});
      this.selectionObjects = [];

      // Content
      selectionBox.append(hideButton.ui);
      selectionBox.append(showArea);
      selectionBox.css('position', 'absolute');

      scrollBox.append(selections);
      showArea.append(scrollBox);
      addArea.append(plusButton.ui);

      const alertBox = new AlertBox();
      showArea.append(alertBox.ui);
      showArea.append(addArea);
      alertBox.ui.css('width', 162);

      // CSS
      if (side === 'left') {
        selectionBox.css('text-align', 'left');
      } else if (side === 'right') {
        selectionBox.css('text-align', 'right');
      } else {
        // Add alert box code
        selectionBox.css('text-align', 'right');
      }

      showArea.css('box-sizing', 'border-box');
      showArea.css('padding', '3px');
      // showArea.css('width', '162px');

      scrollBox.css('max-height', HEIGHT * 0.8);
      scrollBox.css('overflow-y', 'auto');
      scrollBox.css('overflow-x', 'visible');

      selections.css('box-sizing', 'content-box');

      this.updateScrollBox = function (height) {
        scrollBox.css('max-height', height * 0.8);
      };

      // Action
      let hidden = true;
      showArea.hide();

      hideButton.ui.click(toggleHide);

      function toggleHide() {
        if (hidden) {
          showArea.show(100);
        } else {
          showArea.hide(100);
        }
        hidden = !hidden;
      }

      /**
       * Card for manipulation of a selection form and related styles
       *
       * @function Selection
       */
      function Selection() {
        const boundingBox = (this.ui = $('<div></div>'));
        let sid = (this.id = null);
        selectionObjects.push(this);
        boundingBox.css({
          background: '#e8e8e8',
          padding: '4px 4px 2px 4px',
          'border-radius': '6px',
          'margin-bottom': '3px',
          position: 'relative',
          width: '156px',
        });

        const header = $('<div></div>');
        boundingBox.append(header);
        const heading = $('<div></div>');
        const controls = $('<div></div>');

        header.append(heading, controls);
        heading.css({
          'font-family': 'Arial',
          'font-weight': 'bold',
          'font-size': '12px',
          display: 'inline-block',
          width: '60px',
        });

        controls.css({
          display: 'inline-block',
        });

        header.hide();
        controls.editMode = false;

        const removeButton = new button(Icons.minus, 16, {
          bfr: 0.5,
          backgroundColor: '#f06f6f',
          tooltip: 'Remove Selection',
        });
        const editButton = new button(Icons.pencil, 16, {tooltip: 'Edit Selection'});
        const visibleButton = new button(Icons.visible, 16, {tooltip: 'Show / Hide Selection'});

        controls.append(removeButton.ui);
        controls.append(editButton.ui);
        controls.append(visibleButton.ui);

        const parameters = $('<div></div>');
        boundingBox.append(parameters);

        const styleHolder = $('<div></div>');

        removeButton.ui.on('click', () => {
          stateManager.removeSelection(sid);
          boundingBox.detach();
          // delete this;
        });

        editButton.ui.on('click', () => {
          parameters.toggle();
        });

        let hidden = false;
        visibleButton.ui.on('click', () => {
          stateManager.toggleHide(sid);
          if (hidden) {
            hidden = false;
            visibleButton.setSVG(Icons.visible);
          } else {
            hidden = true;
            visibleButton.setSVG(Icons.invisible);
          }
        });

        const styleBox = new StyleBox();

        let showStyle = false;
        styleHolder.append(styleBox.ui);
        styleBox.ui.css({
          position: 'static',
          // 'left' : '0',
          width: 'px',
          'border-radius': '4px',
        });

        styleBox.ui.hide();

        const allControl = (this.allSelector = {
          key: 'Select All Atom',
          value: null,
          active: true,
        });

        const allCheckBox = new Form.Checkbox(allControl);
        parameters.append(allCheckBox.ui);

        /** @type {{key: string; value: any; active:boolean}} */
        const selectionFormControl  = {
          key: 'Selection Spec',
          value: null,
          active: true,
        };

        this.selectionValue = selectionFormControl

        const selectionSpecForm = new Form(
          GLModel.validAtomSelectionSpecs,
          selectionFormControl
        );
        parameters.append(selectionSpecForm.ui);

        const submitControls = $('<div></div>');
        const submit = new button(Icons.tick, 16, {
          backgroundColor: 'lightgreen',
          tooltip: 'Submit',
        });
        const cancel = new button(Icons.cross, 16, {
          backgroundColor: 'lightcoral',
          tooltip: 'Cancel',
        });
        submitControls.append(submit.ui, cancel.ui);

        const alertBox = new AlertBox();
        parameters.append(alertBox.ui);

        parameters.append(submitControls);
        boundingBox.append(styleHolder);

        allCheckBox.update = function () {
          selectionSpecForm.ui.toggle();
        };

        function finalizeSelection(id) {
          header.show();
          controls.editMode = true;
          sid = this.id = id;
          heading.text(`Sel#${id}`);
          boundingBox.attr('data-id', id);
          parameters.hide();
          showStyle = true;
          styleBox.setSid(id);
          styleBox.ui.show();
        }

        function checkAndAddSelection(sid = null) {
          const validate = selectionSpecForm.validate();
          if (validate) {
            selectionSpecForm.getValue();
            const checkAtoms = stateManager.checkAtoms(selectionFormControl.value);

            if (Object.keys(selectionFormControl.value).length === 0) {
              alertBox.error('Please enter some input');
            } else {
              if (checkAtoms) {
                const id = stateManager.addSelection(selectionFormControl.value, sid);
                finalizeSelection(id);
                if (sid === null) editingForm.isEditing = false;
              } else {
                alertBox.error('No atom selected');
              }
            }
          } else {
            alertBox.error('Invalid Input');
          }
        }

        function removeSelf(selection) {
          const selectionToRemove = selectionObjects.find(sel => {
            if (selection === sel) {
              console.log('Selection found', selection);
              return true;
            }
          });

          // delete selectionToRemove;
        }

        submit.ui.on('click', () => {
          if (controls.editMode === false) {
            if (allControl.value) {
              const id = stateManager.addSelection({});
              finalizeSelection(id);
              editingForm.isEditing = false;
            } else {
              checkAndAddSelection();
            }
          } else {
            if (allControl.value) {
              const id = sid;
              stateManager.addSelection({}, id);
              finalizeSelection(id);
            } else {
              const id = sid;
              checkAndAddSelection(id);
            }
          }
        });

        const self = this;

        cancel.ui.on('click', () => {
          if (controls.editMode) {
            parameters.hide();
          } else {
            boundingBox.detach();
            removeSelf(self);
            editingForm.isEditing = false;
          }
        });

        boundingBox.on('keyup', e => {
          if (e.key === 'Enter') {
            submit.ui.trigger('click');
          }
        });

        /**
         * @function Selection#setProperty
         * @param {string|null} id Id of the selection created in StateManager
         * @param {Object} specs Defination of the selection that will be used to set default
         * values in the form
         */
        this.setProperty = function (id, specs) {
          // check for all selection
          if (Object.keys(specs).length === 0) {
            allCheckBox.setValue(true);
          } else {
            selectionSpecForm.setValue(specs);
          }

          // finalize the selection
          finalizeSelection(id);
        };

        /**
         * Adds style to the given selection
         *
         * @function Selection#addStyle
         * @param {string|null} selId Id of the selection to inititate the StyleBox
         * @param {string|null} styleId Id of the style that is created through StateManager
         * @param {import("../specs").AtomStyleSpec} styleSpecs
         */
        this.addStyle = function (selId, styleId, styleSpecs) {
          styleBox.addStyle(selId, styleId, styleSpecs);
        };
      }

      plusButton.ui.on('click', () => {
        if (!editingForm.isEditing) {
          const newSelection = new Selection();
          selections.append(newSelection.ui);
          editingForm.isEditing = true;
        } else {
          alertBox.warning('Please complete the previous form');
        }
      });

      /**
       * Remove all the selection card from the ui
       */
      this.empty = function () {
        selections.empty();
        editingForm.isEditing = false;
      };

      /**
       * Adds or create new selection card
       *
       * @function SelectionBox#editSelection
       * @param {string|null} id Id created in StateManager and passed down to this function during call
       * @param {import("../specs").AtomSelectionSpec} selSpec Selection spec that is used to generate the selection form
       * @param {string|null} styleId Id of style created in StateManager
       * @param {import("../specs").AtomStyleSpec} [styleSpec] Style spec if specified add the selection to the current selection
       */
      this.editSelection = function (id, selSpec, styleId, styleSpec) {
        // if selection does not exist create new

        // This thing works but I am not sure how!

        // Search selection with id
        const selectionUI = selections.children(`[data-id=${id}]`);
        let selection;
        if (!(selectionUI.length !== 0)) {
          selection = new Selection();
          selection.setProperty(id, selSpec);
          selections.append(selection.ui);
        }

        if (styleId !== null && selection && styleSpec) {
          selection.addStyle(id, styleId, styleSpec);
        }
      };
    }

    /**
     *
     * @param {Object} Jquery dom object
     */
    /**
     * Creates StyleBox for listing out different styles inside the selection
     *
     * @function StyleBox
     * @param {String} [selId] Id of the selection for which the style box is created
     * @param {String} [side] Alignment of text inside the box
     */
    function StyleBox(selId, side = 'left') {
      const styleBox = (this.ui = $('<div></div>'));
      editingForm.isEditing = false;
      let sid = (this.sid = selId); // selection id

      this.setSid = function (id) {
        sid = this.sid = id;
      };

      const styles = $('<div></div>');
      const scrollBox = $('<div></div>');

      styles.css('opacity', '0.9');

      const showArea = $('<div></div>');
      const addArea = $('<div></div>');
      addArea.css('text-align', 'center');
      const plusButton = new button(Icons.plus, 20, {tooltip: 'Add New Style'});
      plusButton.ui.css('margin', '0px');

      this.selectionObjects = [];

      // Content
      styleBox.append(showArea);
      styleBox.css('position', 'absolute');

      scrollBox.append(styles);
      showArea.append(scrollBox);

      const alertBox = new AlertBox();
      showArea.append(alertBox.ui);

      addArea.append(plusButton.ui);
      showArea.append(addArea);

      // CSS
      if (side === 'left') {
        styleBox.css('text-align', 'left');
      } else if (side === 'right') {
        styleBox.css('text-align', 'right');
      } else {
        // Add alert box code
        styleBox.css('text-align', 'right');
      }

      showArea.css('box-sizing', 'border-box');
      showArea.css('padding', '3px');
      // showArea.css('width', '162px');
      showArea.css('background-color', '#a4a4a4');
      showArea.css('border-radius', '4px');

      // scrollBox.css('max-height', HEIGHT*0.8);
      scrollBox.css('overflow', 'hidden');

      // styles.css('max-height', HEIGHT*0.8);
      // styles.css('overflow', 'auto');
      styles.css('box-sizing', 'content-box');

      /**
       * Style card to define the value of the style
       *
       * @param {string|null} sid Id of the selction for which the style box is created
       * and this stye will be added under that selection
       */
      function Style(sid) {
        const boundingBox = (this.ui = $('<div></div>'));
        let stid = (this.id = null); // style id
        boundingBox.css({
          background: '#e8e8e8',
          padding: '4px 4px 2px 4px',
          'border-radius': '6px',
          'margin-bottom': '3px',
          position: 'relative',
        });

        const header = $('<div></div>');
        boundingBox.append(header);
        const heading = $('<div></div>');
        const controls = $('<div></div>');

        header.append(heading, controls);
        heading.css({
          'font-family': 'Arial',
          'font-weight': 'bold',
          'font-size': '12px',
          display: 'inline-block',
          width: '60px',
        });

        controls.css({
          display: 'inline-block',
        });

        header.hide();
        controls.editMode = false;

        const removeButton = new button(Icons.minus, 16, {
          bfr: 0.5,
          backgroundColor: '#f06f6f',
          tooltip: 'Remove Style',
        });
        const editButton = new button(Icons.pencil, 16, {tooltip: 'Edit Style'});
        const visibleButton = new button(Icons.visible, 16, {tooltip: 'Show / Hide Style'});

        controls.append(removeButton.ui);
        controls.append(editButton.ui);
        controls.append(visibleButton.ui);

        const parameters = $('<div></div>');
        boundingBox.append(parameters);

        removeButton.ui.on('click', {parent: this, stid}, e => {
          stateManager.removeStyle(sid, stid);
          boundingBox.detach();
          // delete this;
        });

        editButton.ui.on('click', () => {
          parameters.toggle();
        });

        let hidden = false;
        visibleButton.ui.on('click', () => {
          stateManager.toggleHideStyle(sid, stid);
          if (hidden) {
            hidden = false;
            visibleButton.setSVG(Icons.visible);
          } else {
            hidden = true;
            visibleButton.setSVG(Icons.invisible);
          }
        });

        /** @type {{key: string; value: any; active:boolean}} */
        const styleFormControl = {
          key: 'Style Spec',
          value: null,
          active: true,
        };
        this.selectionValue = styleFormControl;

        const styleSpecForm = new Form(GLModel.validAtomStyleSpecs, styleFormControl);
        parameters.append(styleSpecForm.ui);

        const submitControls = $('<div></div>');
        const submit = new button(Icons.tick, 16, {
          backgroundColor: 'lightgreen',
          tooltip: 'Submit',
        });
        const cancel = new button(Icons.cross, 16, {
          backgroundColor: 'lightcoral',
          tooltip: 'Cancel',
        });
        submitControls.append(submit.ui, cancel.ui);

        const alertBox = new AlertBox();
        parameters.append(alertBox.ui);

        parameters.append(submitControls);

        function finalizeStyle(id) {
          header.show();
          controls.editMode = true;
          stid = id;
          heading.text(`Sty#${id}`);
          parameters.hide();
        }

        function checkAndAddStyle(stid = null) {
          const validate = styleSpecForm.validate();
          if (validate) {
            styleSpecForm.getValue();

            if (Object.keys(styleFormControl.value).length === 0) {
              alertBox.error('Please enter some value');
            } else {
              const id = stateManager.addStyle(styleFormControl.value, sid, stid);
              finalizeStyle(id);
              if (stid === null) editingForm.isEditing = false;
            }
          } else {
            alertBox.error('Invalid Input');
          }
        }

        submit.ui.on('click', () => {
          if (controls.editMode === false) {
            checkAndAddStyle();
          } else {
            const id = stid;
            styleSpecForm.getValue();

            if (Object.keys(styleFormControl.value).length === 0) {
              alertBox.error('Please enter some value');
            } else {
              checkAndAddStyle(id);
            }
          }
        });

        cancel.ui.on('click', () => {
          if (controls.editMode) {
            parameters.hide();
          } else {
            boundingBox.detach();
            // delete this;
          }
        });

        boundingBox.on('keyup', e => {
          if (e.key === 'Enter') {
            submit.ui.trigger('click');
          }
        });

        /**
         * @function Style#updateStyle
         * @param {string|null} styleId Id of the style created by StateManager
         * @param {import("../specs").AtomStyleSpec} styleSpec Specs for defining the style and setting default values
         */
        this.updateStyle = function (styleId, styleSpec) {
          styleSpecForm.setValue(styleSpec);

          finalizeStyle(styleId);
        };
      }

      plusButton.ui.on('click', () => {
        if (!editingForm.isEditing) {
          const newStyle = new Style((/** @type {string} */(sid)));
          styles.append(newStyle.ui);
          editingForm.isEditing = true;
        } else {
          alertBox.warning('Please complete editing the current form');
        }
      });

      /**
       * @function StyleBox#addStyle
       * @param {string|null} selectionId Id of the selection for which styles will be created
       * @param {string|null} styleId Id of the style part of the selection
       * @param {import("../specs").AtomStyleSpec} styleSpecs Style specs that will be used to create
       * style for the specified selection and set default values in the Style card
       */
      this.addStyle = function (selectionId, styleId, styleSpecs) {
        const style = new Style(selectionId);
        styles.append(style.ui);
        style.updateStyle(styleId, styleSpecs);
      };
    }

    /**
     * Add alert messages to different panels
     *
     * @function AlertBox
     * @param {Object} [config] Configuraiton for alert box display
     */
    function AlertBox(config) {
      const boundingBox = (this.ui = $('<div></div>'));
      config = config || {};
      const delay = config.delay || 5000;
      const autohide = config.autohide === undefined ? true : config.autohide;

      boundingBox.css({
        'font-family': 'Arial',
        'font-size': '14px',
        padding: '3px',
        'border-radius': '4px',
        'margin-top': '2px',
        'margin-bottm': '2px',
        'font-weight': 'bold',
        'text-align': 'center',
      });

      boundingBox.hide();

      function hide() {
        if (autohide) {
          setTimeout(() => {
            boundingBox.hide();
          }, delay);
        }
      }

      /**
       * Generate Internal alert message
       * @param {String} msg Error Message
       */
      this.error = function (msg) {
        boundingBox.css({
          background: 'lightcoral',
          color: 'darkred',
          border: '1px solid darkred',
        });

        boundingBox.text(msg);
        boundingBox.show();

        hide();
      };

      /**
       * Generates Internal warning message
       * @param {String} msg Warming message
       */
      this.warning = function (msg) {
        boundingBox.css({
          background: '#fff3cd',
          color: '#856409',
          border: '1px solid #856409',
        });

        boundingBox.text(msg);
        boundingBox.show();

        hide();
      };

      /**
       * Generates Internal Info message
       * @param {String} msg Info message
       */
      this.message = function (msg) {
        boundingBox.css({
          background: 'lightgreen',
          color: 'green',
          border: '1px solid green',
        });

        boundingBox.text(msg);
        boundingBox.show();

        hide();
      };
    }

    /**
     * Creates the panel for manipulation of labels on the viewport
     *
     * @function ContextMenu
     */
    function ContextMenu() {
      const boundingBox = (this.ui = $('<div></div>'));

      boundingBox.css('position', 'absolute');
      // boundingBox.css('border', '1px solid black');
      boundingBox.css('border-radius', '3px');
      boundingBox.css('background', '#f1f1f1');
      boundingBox.css('z-index', 99);
      const contentBox = $('<div></div>');
      contentBox.css('position', 'relative');
      boundingBox.css('opacity', '0.85');

      boundingBox.append(contentBox);
      contentBox.css({
        background: '#f1f1f1',
        'border-radius': '4px',
        padding: '4px',
        width: '140px',
      });
      // Context Box
      // Remove Label Button

      const labelMenuStyle = {
        background: '#d3e2ee',
        padding: '2px',
        'font-family': 'Arial',
        'font-weight': 'bold',
        'font-size': '12px',
        'border-radius': '2px',
        // 'margin-top':'3px'
      };

      const removeLabelMenu = $('<div></div>');
      removeLabelMenu.text('Remove Label');
      removeLabelMenu.css(labelMenuStyle);
      removeLabelMenu.css('margin-bottom', '3px');

      contentBox.append(removeLabelMenu);
      removeLabelMenu.hide();

      // Label Property List
      const propertyKeys = Object.keys(GLModel.validAtomSpecs);
      const propertyList = [];
      let propertyObjectList = [];

      propertyKeys.forEach(prop => {
        const propObj = GLModel.validAtomSpecs;
        if (propObj[prop].prop === true) {
          propertyList.push(prop);
        }
      });

      // Property Menu
      const propertyMenu = $('<div></div>');
      contentBox.append(propertyMenu);

      /**
       * Property object used in property menu
       *
       * @function Property
       * @param {String} key Name of the atom property
       * @param {*} value Value of the property
       */
      function Property(key, value) {
        this.row = $('<tr></tr>');
        const propLabelValue = (this.control = {
          key: '',
          value: null,
          active: true,
          name: key,
        });

        this.key = key;
        this.value = value;

        const checkbox = new Form.Checkbox(propLabelValue);
        const checkboxHolder = $('<td></td>');
        checkboxHolder.append(checkbox.ui);
        const keyHolder = $('<td></td>');
        const separatorHolder = $('<td></td>').text(':');
        const valueHolder = $('<td></td>');

        this.row.append(checkboxHolder, keyHolder, separatorHolder, valueHolder);

        keyHolder.text(key);

        if (typeof value === 'number') {
          valueHolder.text(value.toFixed(2));
        } else {
          valueHolder.text(value.replace(/\^/g, ''));
        }

        console.log('Type of value', typeof value, value);
      }

      /**
       * @param {import("../specs").AtomSpec} atom Value of different property of the atom, if the atom has prop : true
       * then that option is made visible in the context menu
       */
      function setProperties(atom) {
        propertyMenu.empty();
        propertyObjectList = [];

        const propertyTable = $('<table></table>');

        propertyList.forEach(prop => {
          const propObj = new Property(prop, atom[prop]);
          propertyTable.append(propObj.row);
          propertyObjectList.push(propObj);
        });

        propertyMenu.append(propertyTable);

        let labelStyle = {
          value: null,
          key: 'Atom Label Style',
        };

        const labelStyleHolder = $('<div><div>');

        labelStyle = $('<div><div>');
        labelStyle.text('Style');
        labelStyle.css({
          display: 'inline-block',
          'font-family': 'Arial',
          'font-size': '14px',
          'margin-right': '6px',
          'margin-left': '6px',
        });

        const stylesForLabel = new Form.ListInput(labelStyle, Object.keys(labelStyles));
        stylesForLabel.ui.css({
          display: 'inline-block',
        });

        stylesForLabel.setValue('milk');

        labelStyleHolder.append(labelStyle, stylesForLabel.ui);
        propertyMenu.append(labelStyleHolder);

        const submit = new button(Icons.tick, 18, {
          backgroundColor: 'lightgreen',
          tooltip: 'Submit',
        });
        const cancel = new button(Icons.cross, 18, {
          backgroundColor: 'lightcoral',
          tooltip: 'Cancel',
        });

        const controlButtons = $('<div></div>');
        controlButtons.append(submit.ui, cancel.ui);
        // controlButtons.css('text-align', 'center');

        const alertBox = new AlertBox();
        propertyMenu.append(alertBox.ui);

        propertyMenu.append(controlButtons);

        submit.ui.on('click', () => {
          const props = processPropertyList();
          const labelStyleValidation = stylesForLabel.validate();

          if (props !== null) {
            if (labelStyleValidation) {
              stateManager.addAtomLabel(props, atom, stylesForLabel.getValue().value);
              stateManager.exitContextMenu(false);
            } else {
              alertBox.error('Select style for label');
            }
          } else {
            alertBox.error('No value selected for label');
          }
        });

        cancel.ui.on('click', () => {
          stateManager.exitContextMenu();
        });
      }

      // Previous Labels
      const labelHolder = $('<div></div>');
      contentBox.append(labelHolder);

      // Add Menu
      const addMenu = $('<div></div>');
      contentBox.append(addMenu);
      addMenu.css('width', '100%');

      const addLabelMenu = $('<div></div>');
      addMenu.append(addLabelMenu);

      addLabelMenu.text('Add Label');
      addLabelMenu.css(labelMenuStyle);
      addLabelMenu.css('margin-bottom', '3px');
      addLabelMenu.hide();

      // Edit Menu
      const editMenu = $('<div></div>');
      contentBox.append(editMenu);

      contentBox.css({
        position: 'relative',
      });

      editMenu.css({
        background: '#dfdfdf',
        'border-radius': '3px',
        'font-family': 'Arial',
        'font-weight': 'bold',
        'font-size': '12px',
        // 'position': 'absolute',
        // 'left' : '105%',
        // 'top' : '0',,
        'box-sizing': 'border-box',
        width: '100%',
      });
      editMenu.hide();

      const alertBox = new AlertBox({autohide: false});
      contentBox.append(alertBox.ui);

      // Add Label Inputs

      /**
       * Generate input elements that are used as form values in the context menu under addLabelForm
       * @returns {Object} that holds different input elements
       */
      function generateAddLabelForm() {
        const addLabelForm = $('<div></div>');

        const addLabelValue = {
          text: {
            key: 'Label Text',
            value: null,
            active: true,
          },

          style: {
            key: 'Style',
            value: null,
            active: true,
          },

          sel: {
            key: 'Selection',
            value: null,
            active: true,
          },
        };
        const formModifierControl = $('<div></div>');
        const removeButton = new button(Icons.minus, 16);
        const tick = new button(Icons.tick, 16, {backgroundColor: 'lightgreen', tooltip: 'Submit'});
        const cross = new button(Icons.cross, 16, {
          backgroundColor: 'lightcoral',
          tooltip: 'Cancel',
        });
        formModifierControl.append(removeButton.ui, tick.ui, cross.ui);
        removeButton.ui.hide();
        addLabelForm.append(formModifierControl);

        const addLabelTextBox = $('<div></div>');
        const lt = $('<div></div>').text('Label Text');
        const addLabelTextInput = new Form.Input(addLabelValue.text);
        addLabelTextBox.append(lt, addLabelTextInput.ui);
        const width = 126; // editMenu.innerWidth()*0.8;
        addLabelTextInput.setWidth(width);
        addLabelForm.append(addLabelTextBox);

        const addLabelStyleBox = $('<div></div>');
        const ls = $('<div></div>').text('Label Style');
        const addLabelStyleInput = new Form.ListInput(
          addLabelValue.style,
          Object.keys(labelStyles)
        );
        addLabelStyleInput.setValue('milk');
        addLabelStyleBox.append(ls, addLabelStyleInput.ui);
        addLabelForm.append(addLabelStyleBox);

        const selectionList = stateManager.getSelectionList();

        const addLabelSelectionBox = $('<div></div>');
        const lsl = $('<div></div>').text('Label Selection');
        const addLabelSelectionInput = new Form.ListInput(addLabelValue.sel, selectionList);
        addLabelSelectionBox.append(lsl, addLabelSelectionInput.ui);
        addLabelForm.append(addLabelSelectionBox);

        // CSS
        addLabelForm.css({
          padding: '2px',
        });

        tick.ui.on('click', () => {
          let validate = true;

          if (!addLabelStyleInput.validate()) validate = false;

          if (!addLabelTextInput.validate()) validate = false;

          if (!addLabelSelectionInput.validate()) validate = false;

          if (validate) {
            stateManager.addLabel(addLabelValue);
          }
        });

        cross.ui.on('click', () => {
          stateManager.exitContextMenu();
        });

        removeButton.ui.on('click', () => {
          stateManager.removeLabel();
        });

        addLabelForm.on('keyup', e => {
          if (e.key === 'Enter') {
            tick.ui.trigger('click');
          }
        });

        return {
          boundingBox: addLabelForm,
          text: addLabelTextInput,
          style: addLabelStyleInput,
          selection: addLabelSelectionInput,
          editMode() {
            removeButton.ui.show();
          },
        };
      }

      function processPropertyList() {
        const propsForLabel = {};

        propertyObjectList.forEach(propObj => {
          if (propObj.control.value === true) {
            propsForLabel[propObj.key] = propObj.value;
          }
        });

        if (Object.keys(propsForLabel).length !== 0) {
          return propsForLabel;
        }
        return null;
      }

      // Context Menu UI Funciton
      boundingBox.hide();
      this.hidden = true;
      this.atom = null;

      removeLabelMenu.on('click', {atom: this.atom}, e => {
        stateManager.removeAtomLabel(removeLabelMenu.atom);
      });

      /**
       * Shows the context menu
       *
       * @function ContextMenu#show
       *
       * @param {Number} x x coordinate of the mouse
       * @param {Number} y y coordinate of the mouse in the viewport in pixels
       * @param {import("../specs").AtomSpec} atom Value of the atoms that is selected
       * @param {Boolean} atomExist if atom label is previously added it is set true else false
       */
      this.show = function (x, y, atom, atomExist) {
        if (atomExist) {
          removeLabelMenu.show();
          removeLabelMenu.atom = atom;
        } else {
          removeLabelMenu.hide();
          removeLabelMenu.atom = null;
        }

        alertBox.ui.hide();
        addLabelMenu.hide();

        if (stateManager.getSelectionList().length === 0) {
          alertBox.message('Please create selections before adding label');
        } else {
          addLabelMenu.show();
        }

        unsetForm();
        setPosition(boundingBox, x, y);
        boundingBox.show();
        this.hidden = false;

        if (atom) {
          setProperties(atom);
          this.atom = atom;
        } else {
          propertyMenu.empty();
        }
      };

      /**
       * Hides the context menu and if needed process the propertyMenu
       *
       * @function ContextMenu#hide
       * @param {Boolean} [processContextMenu] If true then submission of the property to add label is executed
       */

      this.hide = (processContextMenu) => {
        if (processContextMenu) {
          const propsForLabel = processPropertyList();
          if (propsForLabel !== null) {
            stateManager.addAtomLabel(propsForLabel, this.atom);
          }
        }

        boundingBox.hide();
        this.hidden = true;
        unsetForm();
      };

      addLabelMenu.on('click', () => {
        const addLabelMenuForm = generateAddLabelForm();
        setForm(addLabelMenuForm);
      });

      function setForm(form) {
        editMenu.children().detach();
        editMenu.append(form.boundingBox);
        editMenu.show();
      }

      function unsetForm() {
        editMenu.children().detach();
        editMenu.hide();
      }
    }

    /**
     * Creates UI panel for surface manipulations
     *
     * @function SurfaceMenu
     */
    function SurfaceMenu() {
      const boundingBox = (this.ui = $('<div></div>'));
      editingForm.isEditing = false;
      // Selection Layout

      boundingBox.css({
        position: 'absolute',
        width: '140px',
      });

      const surfaceButton = new button(Icons.surface, 20, {tooltip: 'Toggle Surface Menu'});

      boundingBox.append(surfaceButton.ui);

      const displayBox = $('<div></div>');
      boundingBox.append(displayBox);

      // Overflow fix
      boundingBox.css({
        overflow: 'visible',
      });

      const newSurfaceSpace = $('<div></div>');
      newSurfaceSpace.css({
        'max-height': HEIGHT * 0.8,
        'overflow-y': 'auto',
        'overflow-x': 'hidden',
      });

      this.updateScrollBox = function (height) {
        newSurfaceSpace.css('max-height', height * 0.8);
      };
      // newSurfaceSpace.append(controlButton);
      // controlButton.hide();

      displayBox.append(newSurfaceSpace);

      const alertBox = new AlertBox();
      displayBox.append(alertBox.ui);

      const addArea = $('<div></div>');
      const addButton = new button(Icons.plus, 20, {tooltip: 'Add New Surface'});
      addArea.append(addButton.ui);
      displayBox.append(addArea);
      displayBox.hide();

      const surfaces = (this.surfaces = []);

      /**
       * Creates cards for manipulation of surface
       *
       * @function Surface
       */
      function Surface() {
        const control = {
          surfaceType: {
            key: 'Surface Type',
            value: null,
          },
          surfaceStyle: {
            key: 'Surface Style',
            value: null,
          },
          surfaceFor: {
            key: 'Selection Atoms',
            value: null,
          },
          surfaceOf: {
            key: 'Surface Generating Atoms',
            value: null,
          },
        };

        const surfaceBox = (this.ui = $('<div></div>'));
        surfaceBox.css({
          'margin-top': '3px',
          padding: '6px',
          'border-radius': '3px',
          'background-color': '#e8e8e8',
          // 'position':'relative',
          width: '100%',
          'box-sizing': 'border-box',
          // 'left': "-100%",
          opacity: 0.9,
          'text-align': 'left',
        });

        const heading = (this.heading = $('<div></div>'));
        const header = $('<div></div>');

        header.css({
          'text-align': 'right',
        });

        // Control Buttons
        const toolButtons = $('<div></div>');

        const editButton = new button(Icons.pencil, 16, {tooltip: 'Edit Surface'});
        const removeButton = new button(Icons.minus, 16, {bfr: 0.5, backgroundColor: '#f06f6f'});

        toolButtons.append(removeButton.ui);
        toolButtons.append(editButton.ui);

        toolButtons.editButton = editButton;
        toolButtons.removeButton = removeButton;
        toolButtons.editMode = false;

        const defaultTextStyle = {
          'font-weight': 'bold',
          'font-family': 'Arial',
          'font-size': '12px',
        };

        heading.css('display', 'inline-block');
        heading.css(defaultTextStyle);

        toolButtons.css('display', 'inline-block');
        header.hide();

        header.append(heading, toolButtons);
        surfaceBox.append(header);

        // toolButtons.hide();
        const surfacePropertyBox = $('<div></div>');
        surfaceBox.append(surfacePropertyBox);

        // Surface Type
        const surfaceType = $('<div></div>');

        const labelSurfaceType = $('<div></div>');
        labelSurfaceType.text('Surface Type');
        labelSurfaceType.css(defaultTextStyle);

        const listSurfaceType = new Form.ListInput(control.surfaceType, Object.keys(SurfaceType));

        surfaceType.append(labelSurfaceType, listSurfaceType.ui);
        surfacePropertyBox.append(surfaceType);

        listSurfaceType.setValue(Object.keys(SurfaceType)[0]);
        // Surface Style
        const surfaceStyle = $('<div></div>');

        const labelSurfaceStyle = $('<div></div>');
        // labelSurfaceStyle.text('Surface Style');

        const formSurfaceStyle = new Form(GLModel.validSurfaceSpecs, control.surfaceStyle);

        surfaceStyle.append(labelSurfaceStyle, formSurfaceStyle.ui);
        surfacePropertyBox.append(surfaceStyle);

        // Surface Of
        const surfaceOf = $('<div></div>');

        const labelSurfaceOf = $('<div></div>');
        labelSurfaceOf.text('Surface Atoms');
        labelSurfaceOf.css(defaultTextStyle);

        const surfaceGeneratorAtomType = ['self', 'all'];
        const surfaceGeneratorDesc = {
          self: 'Atoms in the selections will be used to generate the surface',
          all: 'All the atoms will be used to generate the surface',
        };

        const listSurfaceOf = new Form.ListInput(control.surfaceOf, surfaceGeneratorAtomType);

        const hintbox = $('<div></div>');
        hintbox.css({
          'background-color': '#e4e4e4',
          border: '1px solid grey',
          color: 'grey',
          padding: '2px',
          'border-radius': '3px',
          'font-family': 'Arial',
          'font-size': '12px',
          'font-weight': 'bold',
          'margin-top': '3px',
        });

        hintbox.hide();

        listSurfaceOf.update = function (control) {
          if (control.value === 'self') {
            hintbox.show();
            hintbox.text(surfaceGeneratorDesc.self);
          } else if (control.value === 'all') {
            hintbox.show();
            hintbox.text(surfaceGeneratorDesc.all);
          } else {
            hintbox.hide();
          }
        };

        listSurfaceOf.setValue('all');

        surfaceOf.append(labelSurfaceOf, listSurfaceOf.ui, hintbox);
        surfacePropertyBox.append(surfaceOf);

        // Surface For
        const selectionListElement = ['all'].concat(stateManager.getSelectionList());
        const surfaceFor = $('<div></div>');

        const labelSurfaceFor = $('<div></div>');
        labelSurfaceFor.text('Show Atoms');
        labelSurfaceFor.css(defaultTextStyle);

        const listSurfaceFor = new Form.ListInput(control.surfaceFor, selectionListElement);
        listSurfaceFor.setValue('all');

        surfaceFor.append(labelSurfaceFor, listSurfaceFor.ui);
        surfacePropertyBox.append(surfaceFor);

        const alertBox = new AlertBox();
        surfacePropertyBox.append(alertBox.ui);

        // Control Button
        const controlButton = $('<div></div>');
        const submit = new button(Icons.tick, 16, {
          backgroundColor: 'lightgreen',
          tooltip: 'Submit',
        });
        const cancel = new button(Icons.cross, 16, {
          backgroundColor: 'lightcoral',
          tooltip: 'Cancel',
        });
        controlButton.append(submit.ui);
        controlButton.append(cancel.ui);
        surfacePropertyBox.append(controlButton);

        // Functionality
        removeButton.ui.on('click', {surfaceBox}, e => {
          const id = e.data.surfaceBox.data('surf-id');
          surfaceBox.remove();
          stateManager.removeSurface(id);
        });

        editButton.ui.on('click', () => {
          surfacePropertyBox.toggle();

          // After creation of the surface box all the changes will be edit to the surfaces so on first submit toolButtons.editMode === true;
        });

        // Form Validation

        const validateInput = (this.validateInput = function () {
          let validated = true;

          if (!listSurfaceFor.validate()) {
            validated = false;
          }

          if (!listSurfaceOf.validate()) {
            validated = false;
          }

          if (!listSurfaceType.validate()) {
            validated = false;
          }

          if (!formSurfaceStyle.validate()) {
            validated = false;
          }

          return validated;
        });

        // edit this code to add on edit selection option to work
        // boundingBox.on('mouseenter', function(){
        //   selections = stateManager.getSelectionList();
        //   selectionListElement = selections.map( (m)=>{
        //     return m.id;
        //   });
        //   listSurfaceFor.updateList(selectionListElement);
        //   listSurfaceOf.updateList(selectionListElement);
        // });

        function finalize(id) {
          // element properties
          surfaceBox.data('surf-id', id);
          heading.text(`surf#${id}`);

          header.show();
          toolButtons.editMode = true;
          surfacePropertyBox.hide();
        }

        // Submit
        submit.ui.on('click', {}, function () {
          listSurfaceFor.getValue();
          listSurfaceOf.getValue();
          listSurfaceType.getValue();
          formSurfaceStyle.getValue();

          if (validateInput()) {
            if (toolButtons.editMode === false) {
              const id = stateManager.addSurface(control);
              control.id = id;

              finalize(id);

              // @ts-ignore ignore this conversion as it is abusing the jquery this type
              surfaces.push(this);
              editingForm.isEditing = false;
            } else {
              formSurfaceStyle.getValue();
              control.id = surfaceBox.data('surf-id');
              console.log('Edit surface called');
              stateManager.editSurface(control); // -> add updateSurface funciton to surfaceMenu
              surfacePropertyBox.hide();
            }
          } else {
            alertBox.error('Invalid Input');
          }
        });

        // Cancel Edit
        cancel.ui.on('click', {}, () => {
          if (toolButtons.editMode === false) {
            surfaceBox.detach();
            surfaceBox.remove();
            editingForm.isEditing = false;
          } else {
            surfacePropertyBox.hide();
            toolButtons.editMode = false;
          }
        });

        surfaceBox.on('keyup', e => {
          if (e.key === 'Enter') {
            submit.ui.trigger('click');
          }
        });

        /**
         * Finalizes the surface card with value specified in the surfaceSpec
         *
         * @function Surface#editSurface
         * @param {String} id Id of the surface generated by StateManager
         * @param {Object} surfaceSpec Different value of the surface menu
         */
        this.editSurface = function (id, surfaceSpec) {
          finalize(id);
          listSurfaceType.setValue(surfaceSpec.surfaceType.value);
          formSurfaceStyle.setValue(surfaceSpec.surfaceStyle.value);
          listSurfaceOf.setValue(surfaceSpec.surfaceOf.value);
          listSurfaceFor.setValue(surfaceSpec.surfaceFor.value);

          listSurfaceFor.getValue();
          listSurfaceOf.getValue();
          listSurfaceType.getValue();
          formSurfaceStyle.getValue();
        };
      }

      // Functionality

      // Surface addition

      addButton.ui.on('click', {surfaces: this}, e => {
        if (!editingForm.isEditing) {
          const newSurface = new Surface();
          newSurfaceSpace.append(newSurface.ui);
          editingForm.isEditing = true;
        } else {
          alertBox.warning('Please complete the previous form first');
        }
      });

      surfaceButton.ui.on('click', () => {
        displayBox.toggle();
      });

      /**
       * Clear all the surface cards
       * @function SurfaceMenu#empty
       */

      this.empty = function () {
        newSurfaceSpace.empty();
        editingForm.isEditing = false;
      };

      /**
       * Add Surface in the Surface Menu
       *
       * @function SurfaceMenu#addSurface
       * @param {String} id Id of the surface generated in the StateManager
       * @param {Object} surfaceSpec Values of different property required for setting values in surface menu
       */
      this.addSurface = function (id, surfaceSpec) {
        const newSurface = new Surface();
        newSurfaceSpace.append(newSurface.ui);

        newSurface.editSurface(id, surfaceSpec);
      };
    }

    /**
     * Creates all the jquery object of different UI features
     * @param  {Object} config
     */
    function generateUI(config) {
      const modelToolBar = new ModelToolbar();
      mainParent.append(modelToolBar.ui);
      setLocation(mainParent, modelToolBar.ui, 'left', 'top');
      // modelToolBar.updateInputLength();

      const contextMenu = new ContextMenu();
      mainParent.append(contextMenu.ui);
      setPosition(contextMenu.ui, 100, 100);

      const surfaceMenu = new SurfaceMenu();
      mainParent.append(surfaceMenu.ui);
      setLocation(mainParent, surfaceMenu.ui, 'right', 'top', 0, modelToolBar.ui.height() + 5);

      const selectionBox = new SelectionBox(Icons.select);
      mainParent.append(selectionBox.ui);
      setLocation(mainParent, selectionBox.ui, 'left', 'top', 0, modelToolBar.ui.height() + 5);

      // Fixing Context Menu Behaviour
      selectionBox.ui.on('mousedown', () => {
        stateManager.exitContextMenu();
      });

      surfaceMenu.ui.on('mousedown', () => {
        stateManager.exitContextMenu();
      });

      return {
        modelToolBar,
        selectionBox,
        contextMenu,
        surfaceMenu,
      };
    }

    /**
     * Sets the css position property left and top for the element
     *
     * @function setPosition
     *
     * @param {object} ele - jquery html element
     * @param {number} left - css left property
     * @param {number} top - css top peroperty
     */
    function setPosition(ele, left, top) {
      ele.css('left', left);
      ele.css('top', top);
    }

    /**
     * Sets the location of the element relative to the parseInt
     * as per position types
     * @function setLocation
     *
     * @param  {Object} parent jquery object
     * @param  {Object} child  jquery object
     * @param  {String} x_type 'left|right'
     * @param  {String} y_type 'top|bottom'
     * @param  {Number} x_offset Offset x values in pixels
     * @param  {Number} y_offset Offset y values in pixels
     */
    function setLocation(
      parent,
      child,
      x_type = 'left',
      y_type = 'top',
      x_offset = 0,
      y_offset = 0
    ) {
      // p_ stands for parent
      child.css('z-index', 99);

      const p_position = parent.position();
      const p_width = getWidth(parent);
      const p_height = getHeight(parent);

      // c_ stand for child
      const c_width = child.outerWidth(); // includes padding and margin
      const c_height = child.outerHeight(); // includes padding and margin

      let padding = parseInt(parent.css('padding').replace('px', ''));
      padding = padding || 0;
      const p_top = getTop(parent) + parseInt(parent.css('margin-top').replace('px', ''));
      const p_left = getLeft(parent) + parseInt(parent.css('margin-left').replace('px', ''));

      // Setting position
      const c_position = {
        left: 0,
        top: 0,
      };

      if (x_type === 'left') {
        c_position.left = padding + x_offset;
      } else if (x_type === 'center') {
        c_position.left = p_width / 2 - c_width / 2 + x_offset;
      } else if (x_type === 'right') {
        c_position.left = p_width - c_width - padding + x_offset;
      } else {
        c_position.left = x_offset + padding;
      }

      if (y_type === 'top') {
        c_position.top = y_offset + padding;
      } else if (y_type === 'center') {
        c_position.top = p_height / 2 - c_height / 2 + y_offset;
      } else if (y_type === 'bottom') {
        c_position.top = p_height - c_height - y_offset - padding;
      } else {
        c_position.top = y_offset + padding;
      }

      setPosition(child, c_position.left, c_position.top);
    }

    // Copied from glviewer.js
    function getRect(container) {
      const div = container[0];
      let rect = div.getBoundingClientRect();
      if (rect.width === 0 && rect.height === 0 && div.style.display === 'none') {
        const oldpos = div.style.position;
        const oldvis = div.style.visibility;
        div.style.display = 'block';
        div.style.visibility = 'hidden';
        div.style.position = 'absolute';
        rect = div.getBoundingClientRect();
        div.style.display = 'none';
        div.style.visibility = oldvis;
        div.style.position = oldpos;
      }
      return rect;
    }

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
  }

  UI.Form = Form;

  return UI;
})();

export default UI;
