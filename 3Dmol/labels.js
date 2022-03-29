// Adapted from the text sprite example from http://stemkoski.github.io/Three.js/index.html

$3Dmol.LabelCount = 0;

/**
 * Renderable labels
 * @constructor $3Dmol.Label
 * @param {string} tag - Label text
 * @param {LabelSpec} parameters Label style and font specifications
 */
$3Dmol.Label = function Label(text, parameters) {
  this.id = $3Dmol.LabelCount++;
  this.stylespec = parameters || {};

  this.canvas = document.createElement('canvas');
  // todo: implement resizing canvas..
  this.canvas.width = 134;
  this.canvas.height = 35;
  this.context = this.canvas.getContext('2d');
  this.sprite = new $3Dmol.Sprite();
  this.text = text;
  this.frame = this.stylespec.frame;
};

$3Dmol.Label.prototype = {
  constructor: $3Dmol.Label,

  getStyle() {
    return this.stylespec;
  },

  setContext: (function setContextConstructor() {
    // function for drawing rounded rectangles - for Label drawing
    const roundRect = function roundRect(ctx, x, y, w, h, r, drawBorder) {
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
    };

    // do all the checks to figure out what color is desired
    const getColor = function getColor(style, stylealpha, init) {
      let ret = init;
      if (typeof style != 'undefined') {
        // convet regular colors
        if (style instanceof $3Dmol.Color) ret = style.scaled();
        else {
          // hex or name
          ret = $3Dmol.CC.color(style);
          if (typeof ret.scaled != 'undefined') {
            ret = ret.scaled(); // not already scaled to 255
          }
        }
      }
      if (typeof stylealpha != 'undefined') {
        ret.a = parseFloat(stylealpha);
      }
      return ret;
    };

    /**
     * Label type specification
     * @typedef LabelSpec
     * @struct
     * @prop {string} font - font name, default sans-serif
     * @prop {number} fontSize - height of text, default 18
     * @prop {ColorSpec} fontColor - font color, default white
     * @prop {number} fontOpacity - font opacity, default 1
     * @prop {number} borderThickness - line width of border around label, default 0
     * @prop {ColorSpec} borderColor - color of border, default backgroundColor
     * @prop {string} borderOpacity - color of border
     * @prop {ColorSpec} backgroundColor - color of background, default black
     * @prop {string} backgroundOpacity - opacity of background, default 1
     * @prop {$3Dmol.Vector3} position - x,y,z coordinates for label
     * @prop {$3Dmol.Vector2} screenOffset - x,y _pixel_ offset of label from position
     * @prop {boolean} inFront - always put labels in from of model
     * @prop {boolean} showBackground - show background rounded rectangle, default true
     * @prop {boolean} fixed - sets the label to change with the model when zooming
     * @prop {boolean} useScreen - position is in screen (not model) coordinates which are pixel offsets from upper left corner.
     * @prop {Object} backgroundImage - An element to draw into the label.  Any CanvasImageSource is allowed.
     * @prop {string} alignment - how to orient the label w/respect to position: topLeft (default), topCenter, topRight, centerLeft, center, centerRight, bottomLeft, bottomCenter, bottomRight
     * @prop {number} frame - if set, only display in this frame of an animation
     */
    return function setContext() {
      const style = this.stylespec;
      const useScreen = typeof style.useScreen == 'undefined' ? false : style.useScreen;

      let {showBackground} = style;
      if (showBackground === '0' || showBackground === 'false') showBackground = false;
      if (typeof showBackground == 'undefined') showBackground = true; // default
      const font = style.font ? style.font : 'sans-serif';

      const fontSize = parseInt(style.fontSize, 10) ? parseInt(style.fontSize, 10) : 18;

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

      let spriteAlignment = style.alignment || $3Dmol.SpriteAlignment.topLeft;
      if (typeof spriteAlignment == 'string' && spriteAlignment in $3Dmol.SpriteAlignment) {
        spriteAlignment = $3Dmol.SpriteAlignment[spriteAlignment];
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
        const w = style.backgroundWidth ? style.backgroundWidth : img.width;
        const h = style.backgroundHeight ? style.backgroundHeight : img.height;
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
        const g = $3Dmol.Gradient.getGradient(style.backgroundGradient);
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
        const w = style.backgroundWidth ? style.backgroundWidth : img.width;
        const h = style.backgroundHeight ? style.backgroundHeight : img.height;
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
      const texture = new $3Dmol.Texture(this.canvas);
      texture.needsUpdate = true;
      this.sprite.material = new $3Dmol.SpriteMaterial({
        map: texture,
        useScreenCoordinates: useScreen,
        alignment: spriteAlignment,
        depthTest: !inFront,
        screenOffset: style.screenOffset || null,
      });

      this.sprite.scale.set(1, 1, 1);

      this.sprite.position.set(position.x, position.y, position.z);
    };
  })(),

  // clean up material and texture
  dispose() {
    if (this.sprite.material.map !== undefined) this.sprite.material.map.dispose();
    if (this.sprite.material !== undefined) this.sprite.material.dispose();
  },
};
