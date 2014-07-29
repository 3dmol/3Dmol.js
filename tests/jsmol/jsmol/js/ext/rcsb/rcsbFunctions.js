// jsmol/ext/JSmolRCSB.js   
//
// a JSmol extension for RCSB-related functionality
// 
// The Jmol.rcsb extension allows easy creation of Java or HTML5 applets on a 
// page that reference PDB 4-characater and 3-character (ligand) IDs.
// This applet's id will be "jmol_" + Info._id
//
// Jmol.rcsb can be used for any URL, not just those from RCSB.
// Options include a label in the upper left-hand corner and 
// one or more buttons below the applet.
//
// Jmol.rcsb.getHtml3D(Info)  // includes _id, _url, _label, _buttonOptions
// Jmol.rcsb.toggleHydrogens(applet)
// Jmol.rcsb.toggleLabels(applet)
// Jmol.rcsb,toggleMep(applet)  // (ligand) molecular electrostatic potential
// Jmol.rcsb.toggleSpin(applet)
//
// Bob Hanson 11/8/2013 6:32:35 AM  hansonr@stolaf.edu
//

;(function(Jmol) {

Jmol.rcsb = {

DefaultInfo : {
  // rcsb-related:
  _id: "1crn",// 1crn, HEM, etc.
  _url: null, // partial or full url
  _label: null,// top left
  _buttonOptions: "",//  hydrogens labels mep spin debug
  
  // Jmol:
  
  script: "",
	width: "100%",
	height: "100%",
	debug: false,
	color: "white",
	use: "HTML5",
  addSelectionOptions: false, // true for debugging
  disableInitialConsole: true,
  j2sPath: "j2s",
	isSigned: true
},

getHtml3D : function(Info) {
  // returns HTML for the applet
  Info || (Info = {});
  Jmol._addDefaultInfo(Info, Jmol.rcsb.DefaultInfo);
  var id = Info._id;
  var loadID = Info._url || (id.length == 3 ? "==" + id : id.length == 4 ? "=" + id : id);  // standard PDB ID for a ligand or protein
  Info.script = "set zoomlarge false;load \"" + loadID + "\";frank off;set antialiasDisplay;" + (Info.script || "");
  if (Info._label)
    Info.script += ";set echo top left;font echo 10;echo " + Info._label + ";color echo blue;";

  var options = (Info._buttonOptions || "");
  if (options && Info.height == "100%")
    Info.height = "90%";
  var appID = "jmol_" + id;
  var html = Jmol.getAppletHtml(appID, Info);
  // based on options specified, create one or more anchors
  var opts = options.split(" ");
  for (var i = 0, o; o=opts[i++];) {
    switch(o){
    case "hydrogens":
      html += '<a href="javascript:Jmol.rcsb.toggleHydrogens(' + appID + ')">hydrogens</a> ';
      break;
    case "labels":
      html += '<a href="javascript:Jmol.rcsb.toggleLabels(' + appID + ')">labels</a> ';
      break;
    case "mep":
      html += '<a href="javascript:Jmol.rcsb.toggleMep(' + appID + ')" title="molecular electrostatic potential">mep</a> ';
      break;
    case "spin":
      html += '<a href="javascript:Jmol.rcsb.toggleSpin(' + appID + ')">spin</a> ';
      break;
    case "debug":
      html += '<a href="javascript:' + appID + '._showInfo(true)">debug</a> ';
      break;
    }
  }
  return html;     
},

// note to Kyle: [16:41:26.611] Jmol.evaluate(jmol_HEM, "script('show _?')")

toggleSpin: function(applet) {
  var s = "if(_spinning){spin off}else{spin on}";
  Jmol.script(applet, s);
},

toggleLabels: function(applet) {
  // turn on labels, but selectively if protein or nucleic
  var s = "font label 10;color labels black;if ({protein|nucleic}){if ({*.CA|*.P}[1].label){labels off} else {select *.CA|*.P;label %n %r;select !protein&!nucleic;label %a} } else if ({*}[1].label){labels off}else{labels %a}";
  Jmol.script(applet, s);
},

toggleMep: function(applet) {
  // a molecule electrostatic potential
  var s = "if (!script('isosurface list').find('visible')){set echo bottom left;echo working...;refresh;select ligand; if ({selected}.partialCharge.max == {selected}.partialCharge.min){calculate partialcharge};isosurface select {ligand} only vdw map mep translucent;echo}else if(script('isosurface list').find('visible:true')){isosurface off}else{isosurface on}";
  Jmol.script(applet, s);
},

toggleHydrogens: function(applet) {
  // we calculate hydrogen positions if this is a protein and there are no H atoms
  var s = 'if({_H}==0&&{*.CA|*.P}){set pdbAddHydrogens true;save orientation o1;load "";restore orientation o1;set pdbAddHydrogens false;}else{showHydrogens=!showHydrogens}';
  Jmol.script(applet, s);
}

} // rcsb

})(Jmol); // closure

