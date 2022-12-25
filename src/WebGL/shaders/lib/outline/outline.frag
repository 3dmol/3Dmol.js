

uniform float opacity;
uniform vec3 outlineColor;
uniform vec3 fogColor;
uniform float fogNear;
uniform float fogFar;
//DEFINEFRAGCOLOR

void main() {
    gl_FragColor = vec4( outlineColor, 1 );
}


