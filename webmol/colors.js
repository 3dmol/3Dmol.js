
WebMol.elementColors = WebMol.elementColors || {};

WebMol.elementColors.defaultColor = 0xff1493;

WebMol.elementColors.Jmol = {
		'H': 0xFFFFFF,
		'He': 0xD9FFFF,
		'HE': 0xD9FFFF,
		'Li': 0xCC80FF,
		'LI': 0xCC80FF,
		'B': 0xFFB5B5,
		'C': 0x909090,
		'N': 0x3050F8,
		'O': 0xFF0D0D,
		'F': 0x90E050,
		'Na': 0xAB5CF2,
		'NA': 0xAB5CF2,
		'Mg': 0x8AFF00,
		'MG': 0x8AFF00,
		'Al': 0xBFA6A6,
		'AL': 0xBFA6A6,
		'Si': 0xF0C8A0,
		'SI': 0xF0C8A0,
		'P': 0xFF8000,
		'S': 0xFFFF30,
		'Cl': 0x1FF01F,
		'CL': 0x1FF01F,
		'Ca': 0x3DFF00,
		'CA': 0x3DFF00,
		'Ti': 0xBFC2C7,
		'TI': 0xBFC2C7,
		'Cr': 0x8A99C7,
		'CR': 0x8A99C7,
		'Mn': 0x9C7AC7,
		'MN': 0x9C7AC7,
		'Fe': 0xE06633,
		'FE': 0xE06633,
		'Ni': 0x50D050,
		'NI': 0x50D050,
		'Cu': 0xC88033,
		'CU': 0xC88033,
		'Zn': 0x7D80B0,
		'ZN': 0x7D80B0,
		'Br': 0xA62929,
		'BR': 0xA62929,
		'Ag': 0xC0C0C0,
		'AG': 0xC0C0C0,
		'I': 0x940094,
		'Ba': 0x00C900,
		'BA': 0x00C900,
		'Au': 0xFFD123,
		'AU': 0xFFD123
};

WebMol.elementColors.rasmol = {
		'H': 0xFFFFFF,
		'He': 0xFFC0CB,
		'HE': 0xFFC0CB,
		'Li': 0xB22222,
		'LI': 0xB22222,
		'B': 0x00FF00,
		'C': 0xC8C8C8,
		'N': 0x8F8FFF,
		'O': 0xF00000,
		'F': 0xDAA520,
		'Na': 0x0000FF,
		'NA': 0x0000FF,
		'Mg': 0x228B22,
		'MG': 0x228B22,
		'Al': 0x808090,
		'AL': 0x808090,
		'Si': 0xDAA520,
		'SI': 0xDAA520,
		'P': 0xFFA500,
		'S': 0xFFC832,
		'Cl': 0x00FF00,
		'CL': 0x00FF00,
		'Ca': 0x808090,
		'CA': 0x808090,
		'Ti': 0x808090,
		'TI': 0x808090,
		'Cr': 0x808090,
		'CR': 0x808090,
		'Mn': 0x808090,
		'MN': 0x808090,
		'Fe': 0xFFA500,
		'FE': 0xFFA500,
		'Ni': 0xA52A2A,
		'NI': 0xA52A2A,
		'Cu': 0xA52A2A,
		'CU': 0xA52A2A,
		'Zn': 0xA52A2A,
		'ZN': 0xA52A2A,
		'Br': 0xA52A2A,
		'BR': 0xA52A2A,
		'Ag': 0x808090,
		'AG': 0x808090,
		'I': 0xA020F0,
		'Ba': 0xFFA500,
		'BA': 0xFFA500,
		'Au': 0xDAA520,
		'AU': 0xDAA520	
};

WebMol.elementColors.defaultColors = WebMol.elementColors.rasmol;

WebMol.elementColors.greenCarbon = $.extend({},WebMol.elementColors.defaultColors);
WebMol.elementColors.greenCarbon['C'] = 0x00ff00;

WebMol.elementColors.cyanCarbon =  $.extend({},WebMol.elementColors.defaultColors);
WebMol.elementColors.cyanCarbon['C'] = 0x00ffff;

WebMol.elementColors.magentaCarbon =  $.extend({},WebMol.elementColors.defaultColors);
WebMol.elementColors.magentaCarbon['C'] = 0xff00ff;

WebMol.elementColors.yellowCarbon =  $.extend({},WebMol.elementColors.defaultColors);
WebMol.elementColors.yellowCarbon['C'] = 0xffff00;

WebMol.elementColors.whiteCarbon =  $.extend({},WebMol.elementColors.defaultColors);
WebMol.elementColors.whiteCarbon['C'] = 0xffffff;

WebMol.elementColors.orangeCarbon =  $.extend({},WebMol.elementColors.defaultColors);
WebMol.elementColors.orangeCarbon['C'] = 0xff6600;

WebMol.elementColors.purpleCarbon =  $.extend({},WebMol.elementColors.defaultColors);
WebMol.elementColors.purpleCarbon['C'] = 0x800080;
