/* jshint -W117:false */
(function (QUnit, env) {
  if (env.__quit_once_initialized) {
    return;
  }
  env.__quit_once_initialized = true;

  if (typeof QUnit !== 'object') {
    throw new Error('undefined QUnit object');
  }

  var _module = QUnit.module;
  if (typeof _module !== 'function') {
    throw new Error('QUnit.module should be a function');
  }

  QUnit.module = function (name, config) {
    if (typeof config !== 'object') {
      return;
    }

    (function addSetupOnce() {
      if (QUnit.supports &&
        QUnit.supports.setupOnce) {
        return;
      }

      if (typeof config.setupOnce === 'function') {
        var _setupOnceRan = false;
        var _setup = typeof config.setup === 'function' ?
          config.setup : null;

        config.setup = function () {
          if (!_setupOnceRan) {
            config.setupOnce();
            _setupOnceRan = true;
          }

          if (_setup) {
            _setup.call(config);
          }
        };
      }
    }());

    (function addTeardownOnce() {

      if (QUnit.supports &&
        QUnit.supports.teardownOnce) {
        return;
      }

      function isLastTestInModule() {
        if (QUnit.config && Array.isArray(QUnit.config.queue)) {
          return QUnit.config.queue.length === 1;
        } else {
          // we cannot determine if the test is the last one in this module
          return false;
        }
      }

      if (typeof config.teardownOnce === 'function') {
        var _teardown = typeof config.teardown === 'function' ?
          config.teardown : null;

        config.teardown = function () {
          if (_teardown) {
            _teardown.call(config);
          }

          if (isLastTestInModule()) {
            config.teardownOnce();
          }
        };
      }
    }());

    _module.call(QUnit, name, config);
  };
}(QUnit, typeof global === 'object' ? global : window));

