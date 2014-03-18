//color scheme mappings
var WebMol = WebMol || {};

//red to white to blue, for charges
WebMol.RWB = function(min, max) {
	
	//map value to hex color, range is provided
	this.valueToHex = function(val, range) {
		var lo, hi;
		if(range) {
			lo = range[0];
			hi = range[1];
		}
		else {
			lo = min;
			hi = max;
		}
	
		if(val === undefined)
			return 0xffffff;
		
		if(val < lo) val = lo;
		if(val > hi) val = hi;
		
		var middle = (hi+lo)/2;
		
		//scale bottom from red to white
		if(val <= middle) {
			var scale = Math.floor(255*Math.sqrt((val-lo)/(middle-lo)));
			var color = 0xff0000 + 0x100*scale + scale;
			return color;
		}
		else { //form white to blue
			var scale = Math.floor(255*Math.sqrt((1-(val-middle)/(hi-middle))));
			var color =  0x10000*scale+0x100*scale+0xff;
			return color;
		}
	};
	
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

//rainbow gradient, but without purple to match jmol
WebMol.ROYGB = function(min, max) {
	
	//map value to hex color, range is provided
	this.valueToHex = function(val, range) {
		var lo, hi;
		if(range) {
			lo = range[0];
			hi = range[1];
		}
		else {
			lo = min;
			hi = max;
		}
	
		if(typeof(val) == "undefined")
			return 0xffffff;
		
		if(val < lo) val = lo;
		if(val > hi) val = hi;
		
		var mid = (lo+hi)/2;
		var q1 = (lo+mid)/2;
		var q3 = (mid+hi)/2;
		
		if(val < q1) { //scale green up, red up, blue down
			var scale = Math.floor(255*Math.sqrt((val-lo)/(q1-lo)));
			var color = 0xff0000 + 0x100*scale + 0;
			return color;
		}
		else if(val < mid) { //scale red down, green up, blue down
			var scale = Math.floor(255*Math.sqrt((1-(val-q1)/(mid-q1))));
			var color =  0x010000*scale+0xff00+0x0;
			return color;
		}
		else if(val < q3) { //scale blue up, red down, green up
			var scale = Math.floor(255*Math.sqrt((val-mid)/(q3-mid)));
			var color = 0x000000 + 0xff00 + 0x1*scale;
			return color;
		}
		else { //scale green down, blue up, red down
			var scale = Math.floor(255*Math.sqrt((1-(val-q3)/(hi-q3))));
			var color =  0x000000+0x0100*scale+0xff;
			return color;
		}		
	};
	
	this.jmolID = function() {
		return "roygb";
	};

	//return range used for color mapping, null if none set
	this.range = function() {
		if(typeof(min) != "undefined" && typeof(max) != "undefined") {
			return [min,max];
		}
		return null;
	};

};

//rainbow gradient with constant saturation, all the way to purple!
WebMol.Sinebow = function(min, max) {
	
	//map value to hex color, range is provided
	this.valueToHex = function(val, range) {
		var lo, hi;
		if(range) {
			lo = range[0];
			hi = range[1];
		}
		else {
			lo = min;
			hi = max;
		}
	
		if(typeof(val) == "undefined")
			return 0xffffff;
		
		if(val < lo) val = lo;
		if(val > hi) val = hi;
		
		var scale = (val-lo)/(hi-lo);
		var h = (5*scale/6.0+.5);
		var r = Math.sin(Math.PI*h);
		r *= r*255;
		var g = Math.sin(Math.PI*(h+1/3.0));
		g *= g*255;
		var b = Math.sin(Math.PI*(h+2/3.0));
		b *= b*255;
		
		return 0x10000*Math.floor(r)+0x100*Math.floor(b)+0x1*Math.floor(g);
	};
	
	this.jmolID = function() {
		return "sinebow";
	};

	//return range used for color mapping, null if none set
	this.range = function() {
		if(typeof(min) != "undefined" && typeof(max) != "undefined") {
			return [min,max];
		}
		return null;
	};

};
