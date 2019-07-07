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
extsrc = ['js/mmtf.js','js/pako_inflate.js','js/netcdfjs.js'];
jqsrc = ['js/jquery-3.2.1.js'];

function clean(cb) {

    del('build/*.js');
    del(['doc/*']);
    cb();
}

function doc(cb) {
    var config = require('./jsdoc.conf.json');
    src(['3Dmol/*.js', 'doc.md'], {read: false})
        .pipe(jsdoc(config, cb));
}

function check(cb) {
    src(coresrc).pipe(jshint({latedef:'nofunc',  esversion:6, laxbreak:true, undef:true, unused:true, 
	    globals: {"$3Dmol":true,
		    'console':true, //set in webworker
		    'document':false,
		    '$':false,
		    'window':false,
		    'module':false,
		    'Blob':false,
		    'XMLHttpRequest':false,
		    'alert':false,
	            'define':false}}))
    .pipe(jshint.reporter('default'));
    cb();
}

function domin(srcs, name) {
	src(srcs).pipe(concat(name+'.js'))
      .pipe(dest('build'))
      .pipe(rename(name+'-min.js'))
      .pipe(uglify().on('error', function(e) { console.log(e);}))
      .pipe(dest('build'));
}

function minify(cb) {
   domin(jqsrc.concat(extsrc).concat(coresrc), '3Dmol');
   cb();
}


function minify_nojquery(cb) {
   domin(extsrc.concat(coresrc), '3Dmol-nojquery');
   cb();
}

function tests(cb) {
  src('tests/auto/generate_tests.py',{read:false}).pipe(shell('python <%= file.path %>'));
  cb();
}

function build_quick(cb) { //nomin
	src(jqsrc.concat(extsrc).concat(coresrc)).pipe(concat('3Dmol.js')).pipe(dest('build'));
	cb();
}
exports.build = series(check, parallel(tests,
                    minify, minify_nojquery, doc));
exports.default = series(clean, exports.build);
exports.build_quick = parallel(build_quick,tests);
exports.clean = clean;
exports.doc = doc;
