
import createViewerGrid from './createViewerGrid';

export default function createStereoViewer(element) {
  if (typeof element === 'string') element = $(`#${element}`);
  if (!element) return;

  const viewers = createViewerGrid(element, {rows: 1, cols: 2, controlAll: true});
  if (!viewers) return;
  /** @type {import("../GLViewer").default} */
  // @ts-ignore nullishness
  this.glviewer1 = viewers[0][0];
  /** @type {import("../GLViewer").default} */
  // @ts-ignore
  this.glviewer2 = viewers[0][1];

  this.glviewer1.setAutoEyeSeparation(false);
  this.glviewer2.setAutoEyeSeparation(true);

  this.glviewer1.linkViewer(this.glviewer2);
  this.glviewer2.linkViewer(this.glviewer1);

  const methods = Object.getOwnPropertyNames(this.glviewer1) // get all methods of glviewer object
    .filter(property => this.glviewer1 && typeof this.glviewer1[property] == 'function');

  for (let i = 0; i < methods.length; i++) {
    // create methods of the same name
    this[methods[i]] = 
      ((method) =>
      (...args) =>
        [
          // eslint-disable-next-line prefer-spread
          this.glviewer1[method].apply(this.glviewer1, ...args),
          // eslint-disable-next-line prefer-spread
          this.glviewer2[method].apply(this.glviewer2, ...args),
        ]
    // @ts-ignore
    )(methods[i]);
  }

  // special cased methods
  this.setCoordinates = function (models, data, format) {
    // for setting the coordinates of the models
    for (let i = 0; i < models.length; i++) {
      models[i].setCoordinates(data, format);
    }
  };

  this.surfacesFinished = function () {
    return this.glviewer1.surfacesFinished() && this.glviewer2.surfacesFinished();
  };

  this.isAnimated = function () {
    return this.glviewer1.isAnimated() || this.glviewer2.isAnimated();
  };

  this.render = function (callback) {
    this.glviewer1.render();
    this.glviewer2.render();
    if (callback) {
      callback(this); // call only once
    }
  };

  this.getCanvas = function () {
    return this.glviewer1.getCanvas(); // same for both
  };
}
