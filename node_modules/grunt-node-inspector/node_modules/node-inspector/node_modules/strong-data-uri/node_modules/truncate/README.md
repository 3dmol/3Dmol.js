Truncate [![Build Status](https://drone.io/github.com/FGRibreau/node-truncate/status.png)](https://drone.io/github.com/FGRibreau/node-truncate/latest)
==================

Truncate text and keeps urls safe.

## NPM
Install the module with: `npm install truncate`

## Usage

```javascript
// Browser
String.truncate("1234 http://google.com hey :)", 2) === "12..."
```

```javascript
// NodeJS
> truncate = require('truncate');
> truncate("1234 http://google.com hey :)", 4);
"1234..."
> truncate("1234 http://google.com hey :)", 4, {ellipsis:null}); // or ellipsis:''
"1234"
> truncate("1234 http://google.com hey :)", 6);
"1234 http://google.com..."
> truncate("1234 http://google.com hey :)", 100);
"1234 http://google.com hey :)"
```

## Release History
v1.0.0 - Initial commit (5 apr. 2012)

## License
Copyright (c) 2013 Francois-Guillaume Ribreau
Licensed under the MIT license.
