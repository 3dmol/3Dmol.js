
module.exports = function(grunt) {
    'use strict';
    
    grunt.initConfig({    
          
        pkg: grunt.file.readJSON('package.json'),
        
        'node-inspector': {
            dev: {}
        },
        
        clean : {
            doc: ['doc'],
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
            },
            webgl : {
                src : ['3Dmol/WebGL/*.js']
            }
            
        },
        
        concat : {
            options : {
                separator : ''
            },
            
            pre : {
                src : ['3Dmol/3Dmol.js', '3Dmol/marchingcube.js', '3Dmol/ProteinSurface4.js', '3Dmol/**.js', '!3Dmol/WebGL/*.js',
                       '!3Dmol/MarchingCubeData.js', '!3Dmol/jmolmodel.js', '!3Dmol/jmolviewer.js'],
                dest : 'build/3Dmol-pre.js'            
            },
            
            webGL : {
                src : ['js/jquery-1.9.1.js', '3Dmol/WebGL/math.js', '3Dmol/WebGL/shapes.js', 
                       '3Dmol/WebGL/core.js', '3Dmol/WebGL/*.js', '3Dmol/properties.js'],
                dest : 'build/webGL-pre.js'
            },
            
            test : {
            	src : ['js/jquery-1.9.1.js', '3Dmol/WebGL/math.js', '3Dmol/WebGL/shapes.js', 
                       '3Dmol/WebGL/core.js', '3Dmol/WebGL/*.js', '3Dmol/properties.js',
                       '3Dmol/3Dmol.js', '3Dmol/marchingcube.js', '3Dmol/ProteinSurface4.js', '3Dmol/**.js',
                       '!3Dmol/MarchingCubeData.js', '!3Dmol/jmolmodel.js', '!3Dmol/jmolviewer.js'],
                dest : 'build/3Dmol-pre.js'
            },
            
            big : {
                src : ['build/webGL-pre.js', 'build/3Dmol-pre.js'],
                dest : 'build/3Dmol.js'
            },
            
            closure : {
                src : ['build/webGL-min-pre.js', 'build/3Dmol-min-pre.js'],
                dest : 'build/3Dmol-min.js'
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
            webGL : {
                src : ['build/webGL-pre.js'],
                dest : 'build/webGL-min-pre.js'
            }
        },
        
        'closure-compiler' : {
            
            $3Dmol : {
                closurePath : 'lib/closure_compiler',
                js : ['build/3Dmol-pre.js'],
                jsOutputFile : 'build/3Dmol-min-pre.js',
                noreport : true,
                options : {
                    'compilation_level': 'ADVANCED_OPTIMIZATIONS',
                    'warning_level': 'DEFAULT',
                    'language_in': 'ECMASCRIPT5',
                    'externs': ['externs/jquery.js', 'externs/3Dmol.js', 'externs/webGL.js'],
                    'create_source_map': 'script.map'                 
                }
            },            
            webgl : {
                closurePath : 'lib/closure_compiler',
                js : ['build/webGL-pre.js'],
                jsOutputFile : 'build/webGL-min-pre.js',
                noreport : true,
                options : {
                    'compilation_level': 'SIMPLE_OPTIMIZATIONS',
                    'warning_level': 'DEFAULT',
                    'language_in': 'ECMASCRIPT5'
                }
            }
        },
        
        shell : {
            
            doc : {
                options : {
                    stdout: true
                },
                command: "node node_modules/jsdoc/jsdoc.js externs/3Dmol.js README.md -c jsdoc.conf.json -t 3Dmol-doc-template -u tutorials/ -d doc/"
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
    grunt.registerTask('concat_pre_build', ['concat:pre', 'concat:webGL']);
    grunt.registerTask('concat_post_build', ['concat:big', 'concat:closure']);
    
    grunt.registerTask('test', ['clean:build', 'concat:test', 'closure-compiler:test', 'concat:append']);
    grunt.registerTask('test_closure', ['clean:build', 'concat_pre_build', 'closure-compiler', 'concat_post_build', 'concat:append']);
    
    grunt.registerTask('build', ['clean:build', 'concat_pre_build', 'closure-compiler', 'concat_post_build', 'clean:tmp']);
    grunt.registerTask('build-quick', ['clean:build', 'concat_pre_build', 'concat_post_build', 'clean:tmp']);

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
