import {
  SpriteAlignment,
  Texture,
  SpriteMaterial,
  Sprite,
  Vector3,
  Vector2,
} from "./WebGL";
import { Gradient } from "./Gradient";
import { Color, CC, ColorSpec } from "./colors";
import {XYZ} from "./WebGL/math"

//Adapted from the text sprite example from http://stemkoski.github.io/Three.js/index.html

export let LabelCount = 0;

// function for drawing rounded rectangles - for Label drawing
function roundRect(ctx, x, y, w, h, r, drawBorder) {
  ctx.beginPath();
  ctx.moveTo(x + r, y);
  ctx.lineTo(x + w - r, y);
  ctx.quadraticCurveTo(x + w, y, x + w, y + r);
  ctx.lineTo(x + w, y + h - r);
  ctx.quadraticCurveTo(x + w, y + h, x + w - r, y + h);
  ctx.lineTo(x + r, y + h);
  ctx.quadraticCurveTo(x, y + h, x, y + h - r);
  ctx.lineTo(x, y + r);
  ctx.quadraticCurveTo(x, y, x + r, y);
  ctx.closePath();
  ctx.fill();
  if (drawBorder) ctx.stroke();
}

//do all the checks to figure out what color is desired
function getColor(style, stylealpha?: any, init?: any) {
  var ret = init;
  if (typeof style != "undefined") {
    //convet regular colors
    if (style instanceof Color) ret = style.scaled();
    else {
      //hex or name
      ret = CC.color(style);
      if (typeof ret.scaled != "undefined") {
        ret = ret.scaled(); //not already scaled to 255
      }
    }
  }
  if (typeof stylealpha != "undefined") {
    ret.a = parseFloat(stylealpha);
  }
  return ret;
}

/** Label style specification */
export interface LabelSpec {
/** font name, default sans-serif */
  font?: string;
  /** height of text, default 18 */
  fontSize?: number;
  /** font color, default white */
  fontColor?: ColorSpec;
  /** font opacity, default 1 */
  fontOpacity?: number;
  /** line width of border around label, default 0  */
  borderThickness?: number;
  /** color of border, default backgroundColor */
  borderColor?: ColorSpec;
  /** opacity of border */
  borderOpacity?: number;
  /** color of background, default black */
  backgroundColor?: ColorSpec;
  /** opacity of background, default 1.0 */
  backgroundOpacity?: number;
  /** coordinates for label */
  position?: XYZ;
  /** x,y pixel offset of label from position */
  screenOffset?: Vector2;
  /** always put labels in front of model */
  inFront?: boolean;
  /** show background rounded rectangle, default true */
  showBackground?: boolean;
  /** position is in screen (not model) coordinates which are pixel offsets from the upper left corner */
  useScreen?: boolean;
  /** An elment to draw into the label. Any CanvasImageSource is allowed.  Label is resized to size of image */
  backgroundImage?: any;
  /** how to orient the label w/respect to position: "topLeft" (default),
   * "topCenter", "topRight", "centerLeft", "center", "centerRight",
   * "bottomLeft", "bottomCenter", "bottomRight", or an arbitrary offset */
  alignment?: string | Vector2;
  /** if set, only display in this frame of an animation */
  frame?: number;
}

/**
 * Renderable labels
 * @constructor $3Dmol.Label
 * @param {string} tag - Label text
 * @param {LabelSpec} parameters Label style and font specifications
 */
export class Label {
  id: number;
  stylespec: any;
  canvas: HTMLCanvasElement;
  context: any;
  sprite: any;
  text: any;
  frame: any;
  constructor(text, parameters) {
    this.id = LabelCount++;
    this.stylespec = parameters || {};

    this.canvas = document.createElement("canvas");
    //todo: implement resizing canvas..
    this.canvas.width = 134;
    this.canvas.height = 35;
    this.context = this.canvas.getContext("2d");
    this.sprite = new Sprite();
    this.text = text;
    this.frame = this.stylespec.frame;
  }

  getStyle() {
    return this.stylespec;
  }

  setContext() {
    var style = this.stylespec;
    var useScreen =
      typeof style.useScreen == "undefined" ? false : style.useScreen;

    var showBackground = style.showBackground;
    if (showBackground === "0" || showBackground === "false")
      showBackground = false;
    if (typeof showBackground == "undefined") showBackground = true; //default
    var font = style.font ? style.font : "sans-serif";

    var fontSize = parseInt(style.fontSize) ? parseInt(style.fontSize) : 18;

    var fontColor = getColor(style.fontColor, style.fontOpacity, {
      r: 255,
      g: 255,
      b: 255,
      a: 1.0,
    });

    var padding = style.padding ? style.padding : 4;
    var borderThickness = style.borderThickness ? style.borderThickness : 0;

    var backgroundColor = getColor(
      style.backgroundColor,
      style.backgroundOpacity,
      {
        r: 0,
        g: 0,
        b: 0,
        a: 1.0,
      }
    );

    var borderColor = getColor(
      style.borderColor,
      style.borderOpacity,
      backgroundColor
    );

    var position = style.position
      ? style.position
      : {
        x: -10,
        y: 1,
        z: 1,
      };

    // Should labels always be in front of model?
    var inFront = style.inFront !== undefined ? style.inFront : true;
    if (inFront === "false" || inFront === "0") inFront = false;

    // clear canvas

    var spriteAlignment = style.alignment || SpriteAlignment.topLeft;
    if (
      typeof spriteAlignment == "string" &&
      spriteAlignment in SpriteAlignment
    ) {
      spriteAlignment = SpriteAlignment[spriteAlignment];
    }

    var bold = "";
    if (style.bold) bold = "bold ";
    this.context.font = bold + fontSize + "px  " + font;

    var metrics = this.context.measureText(this.text);
    var textWidth = metrics.width;

    if (!showBackground) borderThickness = 0;

    var width = textWidth + 2.5 * borderThickness + 2 * padding;
    var height = fontSize * 1.25 + 2 * borderThickness + 2 * padding; // 1.25 is extra height factor for text below baseline: g,j,p,q.

    if (style.backgroundImage) {
      //resize label to image
      var img = style.backgroundImage;
      var w = style.backgroundWidth ? style.backgroundWidth : img.width;
      var h = style.backgroundHeight ? style.backgroundHeight : img.height;
      if (w > width) width = w;
      if (h > height) height = h;
    }

    this.canvas.width = width;
    this.canvas.height = height;
    this.context.clearRect(0, 0, this.canvas.width, this.canvas.height);

    bold = "";
    if (style.bold) bold = "bold ";
    this.context.font = bold + fontSize + "px  " + font;

    // background color
    this.context.fillStyle =
      "rgba(" +
      backgroundColor.r +
      "," +
      backgroundColor.g +
      "," +
      backgroundColor.b +
      "," +
      backgroundColor.a +
      ")";
    // border color
    this.context.strokeStyle =
      "rgba(" +
      borderColor.r +
      "," +
      borderColor.g +
      "," +
      borderColor.b +
      "," +
      borderColor.a +
      ")";

    if (style.backgroundGradient) {
      let gradient = this.context.createLinearGradient(
        0,
        height / 2,
        width,
        height / 2
      );
      let g = Gradient.getGradient(style.backgroundGradient);
      let minmax = g.range();
      let min = -1;
      let max = 1;
      if (minmax) {
        min = minmax[0];
        max = minmax[1];
      }
      let d = max - min;
      for (let i = 0; i < 1.01; i += 0.1) {
        let c = getColor(g.valueToHex(min + d * i));
        let cname = "rgba(" + c.r + "," + c.g + "," + c.b + "," + c.a + ")";
        gradient.addColorStop(i, cname);
      }
      this.context.fillStyle = gradient;
    }

    this.context.lineWidth = borderThickness;
    if (showBackground) {
      roundRect(
        this.context,
        borderThickness,
        borderThickness,
        width - 2 * borderThickness,
        height - 2 * borderThickness,
        6,
        borderThickness > 0
      );
    }

    if (style.backgroundImage) {
      this.context.drawImage(img, 0, 0, width, height);
    }

    // text color
    this.context.fillStyle =
      "rgba(" +
      fontColor.r +
      "," +
      fontColor.g +
      "," +
      fontColor.b +
      "," +
      fontColor.a +
      ")";

    this.context.fillText(
      this.text,
      borderThickness + padding,
      fontSize + borderThickness + padding,
      textWidth
    );

    // canvas contents will be used for a texture
    var texture = new Texture(this.canvas);
    texture.needsUpdate = true;
    this.sprite.material = new SpriteMaterial({
      map: texture,
      useScreenCoordinates: useScreen,
      alignment: spriteAlignment,
      depthTest: !inFront,
      screenOffset: style.screenOffset || null,
    });

    this.sprite.scale.set(1, 1, 1);

    this.sprite.position.set(position.x, position.y, position.z);
  }

  // clean up material and texture
  dispose() {
    if (this.sprite.material.map !== undefined)
      this.sprite.material.map.dispose();
    if (this.sprite.material !== undefined) this.sprite.material.dispose();
  }
}
