// Adapted from the text sprite example from http://stemkoski.github.io/Three.js/index.html

import {CC} from './colors';
import Gradient from './Gradient';
import Color from './WebGL/core/Color';
import {SpriteAlignment, SpriteMaterial, Texture} from './WebGL/materials';
import {Sprite} from './WebGL/objects';

// eslint-disable-next-line import/no-mutable-exports
export let LabelCount = 0;

// do all the checks to figure out what color is desired
function getColor(style, stylealpha, init) {
  let ret = init;
  if (typeof style != 'undefined') {
    // convet regular colors
    if (style instanceof Color) ret = style.scaled();
    else {
      // hex or name
      ret = CC.color(style);
      if (typeof ret.scaled != 'undefined') {
        ret = ret.scaled(); // not already scaled to 255
      }
    }
  }
  if (typeof stylealpha != 'undefined') {
    ret.a = parseFloat(stylealpha);
  }
  return ret;
}

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

/**
 * Renderable labels
 */
export default class Label {
  /**
   * @param {string} text - Label text
   * @param {import("./specs").LabelSpec} parameters Label style and font specifications
   */
  constructor(text, parameters) {
    this.id = LabelCount++;
    this.stylespec = parameters || {};

    this.canvas = document.createElement('canvas');
    // todo: implement resizing canvas..
    this.canvas.width = 134;
    this.canvas.height = 35;
    this.context = this.canvas.getContext('2d');
    this.sprite = new Sprite();
    this.text = text;
    this.frame = this.stylespec.frame;

    if (!this.context) throw new Error('Failed to create canvas context');
  }

  getStyle() {
    return this.stylespec;
  }

  setContext() {
    const style = this.stylespec;
    const useScreen = typeof style.useScreen == 'undefined' ? false : style.useScreen;

    let {showBackground} = style;
    if (showBackground === '0' || showBackground === 'false') showBackground = false;
    if (typeof showBackground == 'undefined') showBackground = true; // default
    const font = style.font ? style.font : 'sans-serif';

    // eslint-disable-next-line no-nested-ternary
    const fontSize = style.fontSize
      ? typeof style.fontSize === 'number'
        ? style.fontSize
        : parseInt(style.fontSize, 10)
      : 18;

    const fontColor = getColor(style.fontColor, style.fontOpacity, {
      r: 255,
      g: 255,
      b: 255,
      a: 1.0,
    });

    const padding = style.padding ? style.padding : 4;
    let borderThickness = style.borderThickness ? style.borderThickness : 0;

    const backgroundColor = getColor(style.backgroundColor, style.backgroundOpacity, {
      r: 0,
      g: 0,
      b: 0,
      a: 1.0,
    });

    const borderColor = getColor(style.borderColor, style.borderOpacity, backgroundColor);

    const position = style.position
      ? style.position
      : {
          x: -10,
          y: 1,
          z: 1,
        };

    // Should labels always be in front of model?
    let inFront = style.inFront !== undefined ? style.inFront : true;
    if (inFront === 'false' || inFront === '0') inFront = false;

    // clear canvas

    let spriteAlignment = style.alignment || SpriteAlignment.topLeft;
    if (typeof spriteAlignment == 'string' && spriteAlignment in SpriteAlignment) {
      spriteAlignment = SpriteAlignment[spriteAlignment];
    }

    let bold = '';
    if (style.bold) bold = 'bold ';
    this.context.font = `${bold + fontSize}px  ${font}`;

    const metrics = this.context.measureText(this.text);
    const textWidth = metrics.width;

    if (!showBackground) borderThickness = 0;

    let width = textWidth + 2.5 * borderThickness + 2 * padding;
    let height = fontSize * 1.25 + 2 * borderThickness + 2 * padding; // 1.25 is extra height factor for text below baseline: g,j,p,q.

    if (style.backgroundImage) {
      const img = style.backgroundImage;
      const w = style.backgroundWidth ? style.backgroundWidth : (/** @type {number} */(img.width));
      const h = style.backgroundHeight ? style.backgroundHeight : (/** @type {number} */(img.height));
      if (w > width) width = w;
      if (h > height) height = h;
    }

    this.canvas.width = width;
    this.canvas.height = height;
    this.context.clearRect(0, 0, this.canvas.width, this.canvas.height);

    bold = '';
    if (style.bold) bold = 'bold ';
    this.context.font = `${bold + fontSize}px  ${font}`;

    // background color
    this.context.fillStyle = `rgba(${backgroundColor.r},${backgroundColor.g},${backgroundColor.b},${backgroundColor.a})`;
    // border color
    this.context.strokeStyle = `rgba(${borderColor.r},${borderColor.g},${borderColor.b},${borderColor.a})`;

    if (style.backgroundGradient) {
      const gradient = this.context.createLinearGradient(0, height / 2, width, height / 2);
      const g = Gradient.getGradient(style.backgroundGradient);
      const minmax = g.range();
      let min = -1;
      let max = 1;
      if (minmax) {
        [min, max] = minmax;
      }
      const d = max - min;
      for (let i = 0; i < 1.01; i += 0.1) {
        const c = getColor(g.valueToHex(min + d * i));
        const cname = `rgba(${c.r},${c.g},${c.b},${c.a})`;
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
      const img = style.backgroundImage;
      const w = style.backgroundWidth ? style.backgroundWidth : (/** @type {number} */(img.width));
      const h = style.backgroundHeight ? style.backgroundHeight : (/** @type {number} */(img.height));
      this.context.drawImage(img, 0, 0, w, h);
    }

    // text color
    this.context.fillStyle = `rgba(${fontColor.r},${fontColor.g},${fontColor.b},${fontColor.a})`;

    this.context.fillText(
      this.text,
      borderThickness + padding,
      fontSize + borderThickness + padding,
      textWidth
    );

    // canvas contents will be used for a texture
    const texture = new Texture(this.canvas);
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
    if (this.sprite.material.map !== undefined) this.sprite.material.map.dispose();
    if (this.sprite.material !== undefined) this.sprite.material.dispose();
  }
}
