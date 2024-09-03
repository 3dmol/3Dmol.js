

uniform float opacity;
uniform vec3 outlineColor;
uniform vec3 fogColor;
uniform float fogNear;
uniform float fogFar;
//DEFINEFRAGCOLOR

void main() {
    gl_FragColor = vec4( outlineColor, 1 );

    float depth = gl_FragCoord.z / gl_FragCoord.w;
    float fogFactor = smoothstep( fogNear, fogFar, depth );
    gl_FragColor = mix( gl_FragColor, vec4( fogColor, gl_FragColor.w ), fogFactor );    
}


