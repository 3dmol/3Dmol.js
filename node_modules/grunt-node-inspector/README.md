# grunt-node-inspector
> Run [node-inspector](https://github.com/node-inspector/node-inspector) as a grunt task for easy configuration and integration with the rest of your workflow

[![NPM version](https://badge.fury.io/js/grunt-node-inspector.png)](http://badge.fury.io/js/grunt-node-inspector)
[![Dependency Status](https://david-dm.org/ChrisWren/grunt-node-inspector.png)](https://david-dm.org/ChrisWren/grunt-node-inspector)

## Getting Started
If you haven't used grunt before, be sure to check out the [Getting Started](http://gruntjs.com/getting-started) guide, as it explains how to create a gruntfile as well as install and use grunt plugins. Once you're familiar with that process, install this plugin with this command:

```bash
npm install grunt-node-inspector --save-dev
```

Then add this line to your project's `Gruntfile.js` gruntfile:

```javascript
grunt.loadNpmTasks('grunt-node-inspector');
```

## Documentation

### Usage

The minimal usage of node-inspector runs with no options:

```js
'node-inspector': {
  dev: {}
}
```

When this is run, node-inspector will be available at `0.0.0.0:8080`.

Here is a config that uses all of the available options for node-inspector:

```js
'node-inspector': {
  custom: {
    options: {
      'web-port': 1337,
      'web-host': 'localhost',
      'debug-port': 5857,
      'save-live-edit': true,
      'no-preload': true,
      'stack-trace-limit': 4,
      'hidden': ['node_modules']
    }
  }
}
```

### Options

#### web-port

Type: `Number` Default: 8080

Port to host the inspector.

#### web-host

Type: `String` Default: '0.0.0.0'

Host to listen on.

#### debug-port

Type: `Number` Default: 5858

Port to connect to the debugging app.

#### save-live-edit

Type: `Boolean` Default: false

Save live edit changes to disk.

#### no-preload

Type: `Boolean` Default: false

Disables preloading *.js to speed up startup

#### stack-trace-limit

Type: `Number` Default: 50

Number of stack frames to show on a breakpoint.

#### hidden

Type: `Array` Default: []

Array of files to hide from the UI (breakpoints in these files will be ignored).

# Changelog

**0.1.4** - Added `--hidden` option for hiding certain files/directories.

**0.1.3** - Bumped node-inspector version to ~0.7.0, adding `--no-preload` option for faster loading.

**0.1.2** - Bumped node-inspector version to ~0.6.0, adding the new `--stack-trace-limit` option. Allowed node-inspector to be listed as a dependency in a project's package.json instead of forcing it to be in grunt-node-inspector's node_modules folder.

**0.1.1** - Bumped node-inspector version to ~0.5.0.

**0.1.0** - Added debug-port and save-live-edit options. Renamed port to web-port and host to web-host to match node-inspector cli naming.

**Breaking changes:**

options.host is now options['web-host'] and options.port is now options['web-port'].

**0.0.1** - Initial release
