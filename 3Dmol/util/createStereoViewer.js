import GLViewer from '../GLViewer';
import createViewerGrid from './createViewerGrid';

export default function createStereoViewer(element) {
  if (typeof element === 'string') element = $(`#${element}`);
  if (!element) throw new Error('No element provided to stereo viewer constructor');
  const viewers = createViewerGrid(element, {rows: 1, cols: 2, controlAll: true});
  if (!viewers) throw new Error('Could not create stereo viewer, createViewerGrid failed');
  const [glviewer1, glviewer2] = viewers[0];
  if (!glviewer1 || !glviewer2)
    throw new Error('Could not create stereo viewer, missing one or more viewers');
  const allV1PropNames = [
    ...Object.getOwnPropertyNames(glviewer1),
    ...Object.getOwnPropertyNames(GLViewer.prototype),
  ];
  const methods = allV1PropNames.filter(
    property => glviewer1 && typeof glviewer1[property] == 'function'
  );
  const properties = allV1PropNames.filter(
    property => glviewer1 && typeof glviewer1[property] != 'function'
  );

  const proxyTarget = {
    glviewer1,
    glviewer2,
    // eslint-disable-next-line class-methods-use-this
    setCoordinates(models, data, format) {
      // for setting the coordinates of the models
      for (let i = 0; i < models.length; i++) {
        models[i].setCoordinates(data, format);
      }
    },

    surfacesFinished() {
      return this.glviewer1.surfacesFinished() && this.glviewer2.surfacesFinished();
    },

    isAnimated() {
      return this.glviewer1.isAnimated() || this.glviewer2.isAnimated();
    },

    render(callback) {
      this.glviewer1.render();
      this.glviewer2.render();
      if (callback) {
        callback(this); // call only once
      }
    },

    getCanvas() {
      return this.glviewer1.getCanvas(); // same for both
    },
  };

  const proxyTargetProperties = Object.getOwnPropertyNames(proxyTarget);
  const stereoViewer = new Proxy(proxyTarget, {
    get(target, property) {
      if (proxyTargetProperties.includes(/** @type {string} */ (property))) return target[property];

      if (methods.includes(/** @type {string} */ (property))) {
        return (...args) => [
          target.glviewer1[property](...args),
          target.glviewer2[property](...args),
        ];
      }

      if (properties.includes(/** @type {string} */ (property))) {
        return target.glviewer1[property];
      }

      return undefined;
    },
    set(target, property, value) {
      if (proxyTargetProperties.includes(/** @type {string} */ (property)))
        target[property] = value;

      if (properties.includes(/** @type {string} */ (property))) {
        target.glviewer1[property] = value;
        target.glviewer2[property] = value;
        return true;
      }
      return false;
    },
  });

  glviewer1.setAutoEyeSeparation(false);
  glviewer2.setAutoEyeSeparation(true);
  glviewer1.linkViewer(glviewer2);
  glviewer2.linkViewer(glviewer1);
  return stereoViewer;
}
