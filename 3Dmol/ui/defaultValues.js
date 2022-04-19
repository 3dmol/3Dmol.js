import { htmlColors } from "../colors";

export const labelStyles = {
  purple : {
    backgroundColor: 0x800080, 
    backgroundOpacity: 0.8
  },

  milk : {
    font : 'Arial',
    fontSize: 12,
    fontColor: htmlColors.black,
    borderThickness: 1,
    borderColor: htmlColors.azure,
    backgroundColor: htmlColors.aliceblue,
    backgroundOpacity: 0.9
  },

  blue : {
    font : 'Arial',
    fontSize: 12,
    fontColor: htmlColors.aliceblue,
    borderThickness: 1,
    borderColor: htmlColors.darkviolet,
    backgroundColor: htmlColors.darkblue,
    backgroundOpacity: 0.9
  },

  chocolate : {
    font : 'Arial',
    fontSize: 12,
    fontColor: htmlColors.aliceblue,
    borderThickness: 1,
    borderColor: htmlColors.brown,
    backgroundColor: htmlColors.chocolate,
    backgroundOpacity: 0.9
  },

  lime : {
    font : 'Arial',
    fontSize: 12,
    fontColor: htmlColors.black,
    borderThickness: 1,
    borderColor: htmlColors.lightgreen,
    backgroundColor: htmlColors.lime,
    backgroundOpacity: 0.9
  },

  rose : {
    font : 'Arial',
    fontSize: 12,
    fontColor: htmlColors.black,
    borderThickness: 1,
    borderColor: htmlColors.mintcream,
    backgroundColor: htmlColors.mistyrose,
    backgroundOpacity: 0.9
  },

  yellow : {

    font : 'Arial',
    fontSize: 12,
    fontColor: htmlColors.black,
    borderThickness: 1,
    borderColor: htmlColors.orange,
    backgroundColor: htmlColors.yellow,
    backgroundOpacity: 0.9
  },

};

export const longPressDuration = 1500;

export const defaultCSS = {
  ListInput : {
    select : {
      'width' : 'auto',
      'border' : 'none',
      'margin' : '0px',
      'font-family' : 'Arial',
      'padding' : '0px',
      'height' : '20px',
      'border-radius' : '4px',
      'box-sizing' : 'border-box'
    }
  },
  Input : {
    input : {
      'margin-bottom' : '0px',
      'padding' : '0px',
      'border-radius' : '4px',
      'font-family' : 'Arial',
      'width' : '96%'
    },

    boundingBox : {
      'margin-left' : '4px',
      'margin-right' : '',
    }
  },
  Checkbox : {},
  Slider : {},
  Color : {},
  TextDefault : {
    'font-family' : 'Arial',
    'margin-left' : '4px'
  },

  LinkBreak : {
    'height' : '3px',
    'border-bottom' : '1px solid #687193'
  }

}