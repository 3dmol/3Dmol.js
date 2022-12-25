
//put the global $3Dmol object into a module
if (typeof module === "object" && typeof module.exports === "object") {
	//for node.js exporting
	module.exports = window.$3Dmol;
}