
module.exports = function(grunt) {
    'use strict';
    
    grunt.initConfig({    
          
        pkg: grunt.file.readJSON('package.json'),
        
        clean : {
            doc: ['doc'],
            build: ['build'],
            tmp: ['build/*pre.js']
        },
        
        jsdoc : {  
            dist : {
                src : ['externs/webmol.js', 'README.md'],
                options : {
                    destination: 'doc',
                    configure: 'jsdoc.conf.json',
                    template: 'node_modules/grunt-jsdoc/node_modules/ink-docstrap/template'
                }                
            }                   
        },
        
        jshint : {
            all : {
                src : ['Gruntfile.js', 'webmol/**.js', '!webmol/jmolmodel.js', '!webmol/jmolviewer.js']
            },
            main : {
                src : ['Gruntfile.js', 'webmol/webmol.js', 'webmol/glcartoon.js', 'webmol/glmodel.js', 'webmol/glviewer.js', 'webmol/glshape.js']    
            },
            aux : {
                src : ['Gruntfile.js', 'webmol/*.js', '!webmol/glcartoon.js', '!webmol/glmodel.js', '!webmol/glviewer.js', '!webmol/glshape.js',
                       '!webmol/jmolmodel.js', '!webmol/jmolviewer.js']    
            },
            webgl : {
                src : ['webmol/WebGL/*.js']
            }
            
        },
        
        concat : {
            options : {
                separator : ''
            },
            
            pre : {
                src : ['webmol/webmol.js', 'webmol/marchingcube.js', 'webmol/ProteinSurface4.js', 'webmol/**.js', '!webmol/WebGL/*.js',
                       '!webmol/MarchingCubeData.js', '!webmol/jmolmodel.js', '!webmol/jmolviewer.js'],
                dest : 'build/webmol-pre.js'            
            },
            
            webGL : {
                src : ['js/jquery-1.9.1.js', 'webmol/WebGL/math.js', 'webmol/WebGL/shapes.js', 
                       'webmol/WebGL/core.js', 'webmol/WebGL/*.js', 'webmol/properties.js'],
                dest : 'build/webGL-pre.js'
            },
            
            big : {
                src : ['build/webGL-pre.js', 'build/webmol-pre.js'],
                dest : 'build/webmol.js'
            },
            
            closure : {
                src : ['build/webGL-min-pre.js', 'build/webmol-min-pre.js'],
                dest : 'build/webmol-min.js'
            },
            
            append : {
                src : ['js/jquery-1.11.1.min.js', 'build/webmol-closure-pre.js', 'append.js'],
                dest : 'build/webmol-closure.js'
            }
        },
        
        uglify : {
            options : {
                mangle : false
            },
            webmol : {
                src : ['build/webmol-pre.js'],
                dest : 'build/webmol-min-pre.js'
            },
            webGL : {
                src : ['build/webGL-pre.js'],
                dest : 'build/webGL-min-pre.js'
            }
        },
        
        'closure-compiler' : {
            
            webmol : {
                closurePath : 'lib/closure_compiler',
                js : ['build/webmol-pre.js'],
                jsOutputFile : 'build/webmol-min-pre.js',
                //noreport : true,
                options : {
                    'compilation_level': 'ADVANCED_OPTIMIZATIONS',
                    'warning_level': 'DEFAULT',
                    'language_in': 'ECMASCRIPT5',
                    'externs': ['externs/jquery.js', 'externs/webmol.js', 'externs/webGL.js'],
                    'create_source_map': 'script.map'                 
                }
            },            
            webgl : {
                closurePath : 'lib/closure_compiler',
                js : ['build/webGL-pre.js'],
                jsOutputFile : 'build/webgl-min-pre.js',
                //noreport : true,
                options : {
                    'compilation_level': 'SIMPLE_OPTIMIZATIONS',
                    'warning_level': 'DEFAULT',
                    'language_in': 'ECMASCRIPT5'
                }
            }
        }   
        
    });
    
    grunt.registerTask('doc', ['clean:doc', 'jsdoc']);
    grunt.registerTask('concat_pre_build', ['concat:pre', 'concat:webGL']);
    grunt.registerTask('concat_post_build', ['concat:big', 'concat:closure']);
    
    grunt.registerTask('test', ['clean:build', 'concat_pre_build', 'closure-compiler:test', 'concat:append']);
    grunt.registerTask('test_closure', ['clean:build', 'concat_pre_build', 'closure-compiler:build', 'concat:append']);
    
    grunt.registerTask('build', ['clean:build', 'concat_pre_build', 'closure-compiler', 'concat_post_build', 'clean:tmp']);
    
    grunt.loadNpmTasks('grunt-jsdoc');
    grunt.loadNpmTasks('grunt-contrib-jshint');
    grunt.loadNpmTasks('grunt-contrib-clean');
    grunt.loadNpmTasks('grunt-contrib-concat');
    grunt.loadNpmTasks('grunt-contrib-uglify');
    grunt.loadNpmTasks('grunt-closure-compiler');
    
};
