
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
        
        uglify : {
            build : {
                src: ['webmol/**.js'],
                dest: 'build/webmol-min.js',
                options: {
                    mangle: false
                }
            }
        }   
    });
    
    grunt.registerTask('doc', ['clean:doc', 'jsdoc']);
    grunt.registerTask('build', ['clean:build', 'uglify']);
    
    grunt.loadNpmTasks('grunt-jsdoc');
    grunt.loadNpmTasks('grunt-contrib-clean');
    grunt.loadNpmTasks('grunt-contrib-uglify');
    
};
