uniform sampler2D colormap;
varying highp vec2 vTexCoords;
uniform vec2 dimensions;
//DEFINEFRAGCOLOR
void main (void) {
   gl_FragColor = texture2D(colormap, vTexCoords);

   //gl_FragColor.g = gl_FragColor.b =  gl_FragColor.r; //debug shading 
}
        