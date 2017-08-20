
module.exports = function(grunt) {
    'use strict';
    
    grunt.initConfig({    
          
        pkg: grunt.file.readJSON('package.json'),
        
        'node-inspector': {
            dev: {}
        },
        
        clean : {
            doc: ['doc/*','!doc/example.html'],
            build: ['build'],
            tmp: ['build/*pre.js'],
            release: ['release']
        },
        
        jshint : {
            all : {
                src : ['Gruntfile.js', '3Dmol/**.js', '!3Dmol/jmolmodel.js', '!3Dmol/jmolviewer.js']
            },
            main : {
                src : ['Gruntfile.js', '3Dmol/3Dmol.js', '3Dmol/glcartoon.js', '3Dmol/glmodel.js', '3Dmol/glviewer.js', '3Dmol/glshape.js', '3Dmol/gldraw.js']    
            },
            aux : {
                src : ['Gruntfile.js', '3Dmol/*.js', '!3Dmol/glcartoon.js', '!3Dmol/glmodel.js', '!3Dmol/glviewer.js', '!3Dmol/glshape.js',
                       '!3Dmol/jmolmodel.js', '!3Dmol/jmolviewer.js']    
            }
            
        },
        
        concat : {
            options : {
                separator : ''
            },
            
            pre : {
                src : ['3Dmol/3dmol.js','3Dmol/WebGL/math.js','3Dmol/WebGL/shapes.js',
                '3Dmol/WebGL/core.js','3Dmol/WebGL/*.js','3Dmol/**.js','!3Dmol/SurfaceWorker.js','3Dmol/SurfaceWorker.js'],
                dest : 'build/3Dmol-pre.js'            
            },            
            
            big : {
                src : ['js/jquery-3.2.1.js','js/mmtf.js','js/pako_inflate.js','js/netcdfjs.js','build/3Dmol-pre.js'],
                dest : 'build/3Dmol.js'
            },
            
            bignojquery : {
                src : ['js/mmtf.js','js/pako_inflate.js','js/netcdfjs.js','build/3Dmol-pre.js'],
                dest : 'build/3Dmol-nojquery.js'
            },
            
            closure : {
                src : ['build/jquery-3.2.1-min-pre.js','build/mmtf-min-pre.js','build/pako_inflate-min-pre.js','build/netcdfjs-min-pre.js','build/3Dmol-min-pre.js'],
                dest : 'build/3Dmol-min.js'
            },
            closurenojquery: {
                src : ['build/mmtf-min-pre.js','build/pako_inflate-min-pre.js','build/netcdfjs-min-pre.js','build/3Dmol-min-pre.js'],
                dest : 'build/3Dmol-nojquery-min.js'
            }, 
            append : {
                src : ['build/3Dmol-min.js', 'append.js'],
                dest : 'build/3Dmol-min.js'
            }
        },
        
        uglify : {
            options : {
                mangle : false
            },
            $3Dmol : {
                src : ['build/3Dmol-pre.js'],
                dest : 'build/3Dmol-min-pre.js'
            },
            jquery : {
                src : ['js/jquery-3.2.1.js'],
                dest : 'build/jquery-3.2.1-min-pre.js'
            },
            mmtf : {
                src : ['js/mmtf.js'],
                dest : 'build/mmtf-min-pre.js'
            },
	    pako : {
		src : ['js/pako_inflate.js'],
                dest : 'build/pako_inflate-min-pre.js'
            },
            netcdfjs : {
                src : ['js/netcdfjs.js'],
                dest : 'build/netcdfjs-min-pre.js'
            },
        },
        
        'closure-compiler' : {
            
            $3Dmol : {
                closurePath : 'lib/closure_compiler',
                js : ['build/3Dmol-pre.js'],
                jsOutputFile : 'build/3Dmol-min-pre.js',
                noreport : true,
		maxBuffer:500,
                options : {
                    'compilation_level': 'SIMPLE_OPTIMIZATIONS',
                    'warning_level': 'DEFAULT',
                    'language_in': 'ECMASCRIPT6',
                    'create_source_map': 'script.map'                 
                }
            },            
            jquery : {
                closurePath : 'lib/closure_compiler',
                js : ['js/jquery-3.2.1.js'],
                jsOutputFile : 'build/jquery-3.2.1-min-pre.js',
                noreport : true,
		maxBuffer:500,
                options : {
                    'compilation_level': 'SIMPLE_OPTIMIZATIONS',
                    'warning_level': 'DEFAULT',
                    'language_in': 'ECMASCRIPT6'
                }
            },
            mmtf : {
                closurePath : 'lib/closure_compiler',
                js : ['js/mmtf.js'],
                jsOutputFile : 'build/mmtf-min-pre.js',
                noreport : true,
		maxBuffer:500,
                options : {
                    'compilation_level': 'SIMPLE_OPTIMIZATIONS',
                    'warning_level': 'DEFAULT',
                    'language_in': 'ECMASCRIPT6'
                }
            },
            pako : {
                closurePath : 'lib/closure_compiler',
                js : ['js/pako_inflate.js'],
                jsOutputFile : 'build/pako_inflate-min-pre.js',
                noreport : true,
		maxBuffer:500,
                options : {
                    'compilation_level': 'SIMPLE_OPTIMIZATIONS',
                    'warning_level': 'DEFAULT',
                    'language_in': 'ECMASCRIPT6'
                }
            },           
        },
        
        shell : {
            
            doc : {
                options : {
                    stdout: true
                },
                command: "node node_modules/jsdoc/jsdoc.js 3Dmol/*.js doc.md -c jsdoc.conf.json -t 3Dmol-doc-template -u tutorials/ -d doc/"
            },
             pythonServer: {
            options: {
                stdout: true
            },
            command: 'python tests/auto/generate_tests.py'
        }
        },

        copy : {
            release : {
                expand : true,
                src : 'build/*.js',
                dest : 'release/',
                flatten : 'true'
            }
        }
        
    });
    
    grunt.registerTask('doc', ['clean:doc', 'shell:doc']);
    grunt.registerTask('concat_pre_build', ['concat:pre']);
    grunt.registerTask('concat_post_build', ['concat:big', 'concat:bignojquery', 'concat:closure', 'concat:closurenojquery']);
    
    grunt.registerTask('test', ['shell:pythonServer']);
    grunt.registerTask('test_closure', ['clean:build', 'concat_pre_build', 'closure-compiler', 'concat_post_build', 'concat:append']);
    
    grunt.registerTask('build', ['clean:build', 'clean:doc', 'concat_pre_build', 'closure-compiler', 'concat_post_build', 'shell:doc', 'clean:tmp','shell:pythonServer']);
    grunt.registerTask('build-quick', ['clean:build', 'concat_pre_build', 'concat_post_build', 'clean:tmp']);
    grunt.registerTask('build-noclean', ['concat_pre_build', 'closure-compiler', 'concat_post_build', 'shell:doc', 'clean:tmp','shell:pythonServer']);

    grunt.registerTask('release-update', ['clean:release', 'build', 'copy:release']);
    
    grunt.registerTask('debug-doc', ['shell:default']);
    
    grunt.loadNpmTasks('grunt-contrib-jshint');
    grunt.loadNpmTasks('grunt-contrib-clean');
    grunt.loadNpmTasks('grunt-contrib-concat');
    grunt.loadNpmTasks('grunt-contrib-uglify');
    grunt.loadNpmTasks('grunt-closure-compiler');
    grunt.loadNpmTasks('grunt-shell');
    grunt.loadNpmTasks('grunt-node-inspector');
    grunt.loadNpmTasks('grunt-contrib-copy');
    
};
