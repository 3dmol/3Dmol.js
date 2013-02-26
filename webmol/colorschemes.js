//color scheme mappings
var WebMol = WebMol || {};

//red to white to blue, for charges
WebMol.RWB = function(min, max) {
	
	//map value to hex color, range is provided
	this.valueToHex = function(val, range) {
		var lo = range[0];
		var hi = range[1];
		
		if(typeof(val) == "undefined")
			return 0xffffff;
		
		if(val < lo) val = lo;
		if(val > hi) val = hi;
		
		var middle = (hi+lo)/2;
		
		//scale bottom from red to white
		if(val <= middle) {
			var scale = Math.floor(255*(val-lo)/(middle-lo));
			var color = 0xff0000 + 0x100*scale + scale;
			return color;
		}
		else { //form white to blue
			var scale = 255-Math.floor(255*(val-middle)/(hi-middle));
			var color =  0x10000*scale+0x100*scale+0xff;
			return color;
		}
	}
	
	this.jmolID = function() {
		return "rwb";
	};

	//return range used for color mapping, null if none set
	this.range = function() {
		if(typeof(min) != "undefined" && typeof(max) != "undefined") {
			return [min,max];
		}
		return null;
	};

};