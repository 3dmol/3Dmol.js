// @param {startState} object it store the initial values to initiliaze the
// UI state
// Everytime UI state is changed this state is updated so that this can
// preserve the changes in case error happens

// This will also help in migrating the controls from to different

$3Dmol.UI.States = (function(){
  function States(startState){
    this.position.x ;
    this.position.y ; //=
    this.height; // =
    this.width ; //=
    this.viewerState ; //=
    this.viewToolbar ; //=
    this.viewSidebar ; //=
    this.viewDialogBox ; //=
    this.viewAnimationControl ; //=
    this.viewSelectionBox ; //=
    this.viewStyleBox ; //=
    this.controlTyle ; //=
    this.toolTip ; //=

  }

  return States;
})()
