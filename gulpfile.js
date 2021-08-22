const { series, parallel, src, watch, dest} = require('gulp');
jshint = require('gulp-jshint'),
concat = require('gulp-concat'),
uglify = require('gulp-terser');
rename = require('gulp-rename');
merge = require('gulp-merge');
shell = require('gulp-shell');
del = require('del');
jsdoc = require('gulp-jsdoc3');

coresrc = ['3Dmol/3dmol.js','3Dmol/WebGL/math.js','3Dmol/WebGL/shapes.js','3Dmol/WebGL/core.js','3Dmol/WebGL/**.js','3Dmol/**.js','!3Dmol/SurfaceWorker.js','3Dmol/SurfaceWorker.js'];
extsrc = ['js/mmtf.js','node_modules/pako/dist/pako_inflate.js','node_modules/netcdfjs/dist/netcdfjs.js'];
uisrc = ['3Dmol/ui/ui.js', '3Dmol/ui/state.js', '3Dmol/ui/icon.js', '3Dmol/ui/form.js', '3Dmol/ui/defaultValues.js'];
jqsrc = ['node_modules/jquery/dist/jquery.js'];

function clean(cb) {

    del('build/*.js');
    del(['doc/*']);
    cb();
}

function doc(cb) {
    var config = require('./jsdoc.conf.json');
    return src(['3Dmol/*.js', '3Dmol/ui/*', 'doc.md'], {read: false})
        .pipe(jsdoc(config, cb));
}

function check() {
    return src(coresrc).pipe(jshint({latedef:'nofunc',  esversion:6, laxbreak:true, undef:true, unused:true,
	    globals: {"$3Dmol":true,
		    'console':true, //set in webworker
		    'document':false,
		    '$':false,
		    'window':false,
		    'self': true, //for passing data to workers?
		    'module':false,
		    'Blob':false,
		    'pako':false,
		    'netcdfjs':false,
		    'XMLHttpRequest':false,
		    'alert':false,
		    'setTimeout':false,
		    'clearTimeout':false,
		    'setInterval':false,
		    'clearInterval': false,
		    'Worker': false,
		    'MMTF': false,
		    'TextDecoder': false,
	            'define':false}}))
    .pipe(jshint.reporter('default'));
}

function domin(srcs, name) {
  return src(srcs).pipe(concat(name+'.js'))
      .pipe(dest('build'))
      .pipe(rename(name+'-min.js'))
      .pipe(uglify().on('error', function(e) { console.log(e);}))
      .pipe(dest('build'));
}

function minify() {
   return domin(jqsrc.concat(extsrc).concat(coresrc).concat(uisrc), '3Dmol');
}


function minify_nojquery() {
   return domin(extsrc.concat(coresrc).concat(uisrc), '3Dmol-nojquery');
}

function tests(cb) {
  src('tests/auto/generate_tests.py',{read:false}).pipe(shell('python3 <%= file.path %>'));
  cb();
}

function build_quick() { //nomin
	return src(jqsrc.concat(extsrc).concat(coresrc).concat(uisrc)).pipe(concat('3Dmol.js')).pipe(dest('build'));
}
exports.build = series(check, parallel(tests,
                    minify, minify_nojquery));
exports.default = series(clean, parallel(exports.build, doc));
exports.build_quick = parallel(build_quick,tests);
exports.clean = clean;
exports.doc = doc;
