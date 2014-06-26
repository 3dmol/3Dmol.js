
module.exports = function(grunt) {
    
    grunt.initConfig({    
          
        pkg: grunt.file.readJSON('package.json'),
        
        clean : {
            doc: ['doc'],
            build: ['build']
        },
        
        jsdoc : {  
            dist : {
                src : ['webmol/**.js', 'README.md'],
                options : {
                    destination: 'doc',
                    configure: '.jsdocrc',
                    template: 'node_modules/grunt-jsdoc/node_modules/ink-docstrap/template'
                }                
            }                   
        },
        
        jshint : {
            main : {
                src : ['Gruntfile.js', 'webmol/*.js']    
            },
            webgl : {
                src : ['webmol/WebGL/*.js', '!webmol/WebGL/renderer.js']
            },
            renderer : {
                src : ['webmol/WebGL/renderer.js']
            }
            
        },
        
        concat : {
            options : {
                separator : ''
            },
            
            dist : {
                src : ['js/jquery-1.9.1.js', 'webmol/webmol.js', 'webmol/WebGL/math.js', 'webmol/WebGL/shapes.js', 
                       'webmol/WebGL/core.js', 'webmol/WebGL/*.js', 'webmol/**.js', 
                       '!webmol/MarchingCubeData.js', '!webmol/jmolmodel.js', '!webmol/jmolviewer.js'],
                dest : 'build/webmol.js'
            }   
        },
        
        uglify : {
            build : {
                src : ['build/webmol.js'],
                dest : 'build/webmol-min.js',
                options : {
                    mangle: false
                }
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
    grunt.registerTask('build', ['clean:build', 'concat', 'uglify']);
    grunt.registerTask('compile', ['clean:build', 'concat', 'closure-compiler']);
    
    grunt.loadNpmTasks('grunt-jsdoc');
    grunt.loadNpmTasks('grunt-contrib-jshint');
    grunt.loadNpmTasks('grunt-contrib-clean');
    grunt.loadNpmTasks('grunt-contrib-concat');
    grunt.loadNpmTasks('grunt-contrib-uglify');
    grunt.loadNpmTasks('grunt-closure-compiler');
    
};
