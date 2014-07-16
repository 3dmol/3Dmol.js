/*
 * grunt-inspector
 * https://github.com/ChrisWren/grunt-inspector
 *
 * Copyright (c) 2013 Chris Wren
 * Licensed under the MIT license.
 */
/*jshint node: true*/
module.exports = function (grunt) {
  'use strict';
  grunt.registerMultiTask('node-inspector', 'Runs node-inspector to debug your node.js JavaScripts', function () {
    var options = this.options();
    var done = this.async();
    var args = [require.resolve('node-inspector/bin/inspector')];
    [
      'web-port',
      'web-host',
      'debug-port',
      'save-live-edit',
      'readTimeout',
      'stack-trace-limit',
      'no-preload',
      'hidden'
    ].forEach(function (option) {
      if (option in options) {
        args.push('--' + option);
        if (option === 'hidden') {
          args.push(JSON.stringify(options[option]));
          return;
        }
        args.push(options[option]);
      }
    });

    grunt.util.spawn({
      cmd: 'node',
      args: args,
      opts: {
        stdio: 'inherit'
      }
    },
    function (error) {
      if (error) {
        grunt.fail.fatal(error);
      }
      done();
    });
  });
};
