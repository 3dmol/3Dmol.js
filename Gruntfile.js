
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
        }   
        
    });
    
    grunt.registerTask('doc', ['clean:doc', 'jsdoc']);
    grunt.registerTask('build', ['clean:build', 'concat', 'uglify']);
    
    grunt.loadNpmTasks('grunt-jsdoc');
    grunt.loadNpmTasks('grunt-contrib-clean');
    grunt.loadNpmTasks('grunt-contrib-concat');
    grunt.loadNpmTasks('grunt-contrib-uglify');
    
};
