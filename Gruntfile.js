
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
                src : ['webmol/**.js', 'README.md'],
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
                src : ['webmol/webmol.js', 'webmol/marchingcube.js', 'webmol/ProteinSurface4.js', 'webmol/WebGL/math.js', 'webmol/WebGL/shapes.js', 
                       'webmol/WebGL/core.js', 'webmol/WebGL/*.js', 'webmol/**.js',
                       '!webmol/MarchingCubeData.js', '!webmol/jmolmodel.js', '!webmol/jmolviewer.js'],
                dest : 'build/webmol-pre.js'            
            },
            
            big : {
                src : ['js/jquery-1.11.1.min.js', 'build/webmol-pre.js'],
                dest : 'build/webmol.js'
            },
            
            min : {
                src : ['js/jquery-1.11.1.min.js', 'build/webmol-min-pre.js'],
                dest : 'build/webmol-min.js'
            },
            
            closure : {
                src : ['js/jquery-1.11.1.min.js', 'build/webmol-min-closure-pre.js'],
                dest : 'build/webmol-min-closure.js'
            }
        },
        
        uglify : {
            options : {
                mangle : false
            },
            build : {
                src : ['build/webmol-pre.js'],
                dest : 'build/webmol-min-pre.js'
            }
        },
        
        'closure-compiler' : {
            build : {
                closurePath : 'lib/closure_compiler',
                js : ['build/webmol-pre.js'],
                jsOutputFile : 'build/webmol-min-closure-pre.js',
                noreport : true,
                options : {
                    'compilation_level': 'ADVANCED_OPTIMIZATIONS',
                    'warning_level': 'DEFAULT',
                    'language_in': 'ECMASCRIPT5',
                    'externs': 'externs/externs.js'
                }
            }
        }   
        
    });
    
    grunt.registerTask('doc', ['clean:doc', 'jsdoc']);
    grunt.registerTask('concat_build', ['concat:big', 'concat:min', 'concat:closure']);
    grunt.registerTask('build', ['clean:build', 'concat:pre', 'uglify', 'closure-compiler', 'concat_build', 'clean:tmp']);
    grunt.registerTask('compile', ['clean:build', 'concat', 'closure-compiler']);
    
    grunt.loadNpmTasks('grunt-jsdoc');
    grunt.loadNpmTasks('grunt-contrib-jshint');
    grunt.loadNpmTasks('grunt-contrib-clean');
    grunt.loadNpmTasks('grunt-contrib-concat');
    grunt.loadNpmTasks('grunt-contrib-uglify');
    grunt.loadNpmTasks('grunt-closure-compiler');
    
};
