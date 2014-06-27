
module.exports = function(grunt) {
    'use strict';
    
    grunt.initConfig({    
          
        pkg: grunt.file.readJSON('package.json'),
        
        clean : {
            doc: ['doc'],
            build: ['build'],
            tmp: ['build/tmp.js']
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
            
            dist : {
                src : ['webmol/webmol.js', 'webmol/marchingcube.js', 'webmol/ProteinSurface4.js', 'webmol/WebGL/math.js', 'webmol/WebGL/shapes.js', 
                       'webmol/WebGL/core.js', 'webmol/WebGL/*.js', 'webmol/**.js', 
                       '!webmol/MarchingCubeData.js', '!webmol/jmolmodel.js', '!webmol/jmolviewer.js'],
                dest : 'build/webmol.js'
            }   
        },
        
        uglify : {
            options : {
                mangle : false
            },
            build : {
                src : ['js/jquery-1.9.1.js', 'build/webmol.js'],
                dest : 'build/webmol-min.js'
            }
        },
        
        'closure-compiler' : {
            build : {
                js : ['build/webmol.js'],
                jsOutputFile : 'build/webmol-min-closure.js',
                noreport : true,
                options : {
                    compilation_levl: 'SIMPLE_OPTIMIZATIONS',
                    warning_level: 'DEFAULT'
                }
            }
        }   
        
    });
    
    grunt.registerTask('doc', ['clean:doc', 'jsdoc']);
    grunt.registerTask('build', ['clean:build', 'concat', 'uglify:build']);
    grunt.registerTask('compile', ['clean:build', 'concat', 'closure-compiler']);
    
    grunt.loadNpmTasks('grunt-jsdoc');
    grunt.loadNpmTasks('grunt-contrib-jshint');
    grunt.loadNpmTasks('grunt-contrib-clean');
    grunt.loadNpmTasks('grunt-contrib-concat');
    grunt.loadNpmTasks('grunt-contrib-uglify');
    grunt.loadNpmTasks('grunt-closure-compiler');
    
};
