const { series, parallel, src, watch, dest } = require("gulp");

const jshint = require("gulp-jshint");
const concat = require("gulp-concat");
const uglify = require("gulp-terser");
const rename = require("gulp-rename");
const merge = require("gulp-merge");
const shell = require("gulp-shell");
const del = require("del");
const jsdoc = require("gulp-jsdoc3");

tssrc = ["tmp/index.js"];

extsrc = [
  "vendor/disable_amd.js",
  "vendor/mmtf.js",
  "../../node_modules/pako/dist/pako.js",
  "../../node_modules/netcdfjs/dist/netcdfjs.js",
  "../../node_modules/upng-js/UPNG.js",
];

jqsrc = ["../../node_modules/jquery/dist/jquery.js"];

function clean(cb) {
  del("build/*.js");
  del(["doc/*"]);
  cb();
}

function doc(cb) {
  var config = require("./jsdoc.conf.json");
  return src(["src/*.js", "src/ui/*", "doc.md"], { read: false }).pipe(
    jsdoc(config, cb)
  );
}

function domin(srcs, name) {
  return src(srcs)
    .pipe(concat(name + ".js"))
    .pipe(dest("build"))
    .pipe(rename(name + "-min.js"))
    .pipe(
      uglify().on("error", function (e) {
        console.log(e);
      })
    )
    .pipe(dest("build"));
}

function minify() {
  return domin(
    jqsrc.concat(extsrc).concat(tssrc),
    "3Dmol"
  );
}

function minify_nojquery() {
  return domin(
    extsrc.concat(tssrc),
    "3Dmol-nojquery"
  );
}

function tests(cb) {
  src("tests/auto/generate_tests.py", { read: false }).pipe(
    shell("python3 <%= file.path %>")
  );
  cb();
}

function build_quick() {
  //nomin
  return src(jqsrc.concat(extsrc).concat(tssrc))
    .pipe(concat("3Dmol.js"))
    .pipe(dest("build"));
}

//const webpack = shell.task("webpack --config webpack.config.js");
//todo: separate dev/prod
const webpack = require('webpack')
const webpackDev = require('./webpack.dev.js')
const webpackProd = require('./webpack.prod.js')
function assets(cb) {
    return new Promise((resolve, reject) => {
        webpack(webpackDev, (err, stats) => {
            if (err) {
                return reject(err)
            }
            if (stats.hasErrors()) {
                return reject(new Error(stats.compilation.errors.join('\n')))
            }
            resolve()
        })
    })
}
exports.build = series(assets, parallel(tests, minify, minify_nojquery));
exports.default = series(clean, parallel(exports.build, doc));
exports.build_quick = series(assets, parallel(build_quick, tests));
exports.clean = clean;
exports.doc = doc
