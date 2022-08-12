uniform mat4 viewMatrix;
uniform float opacity;

uniform vec3 fogColor;
uniform float fogNear;
uniform float fogFar;

varying vec3 vLightFront;
varying vec3 vColor;
//DEFINEFRAGCOLOR

void main() {

    gl_FragColor = vec4( vec3 ( 1.0 ), opacity );

    #ifndef WIREFRAME
    gl_FragColor.xyz *= vLightFront;
    #endif

    gl_FragColor = gl_FragColor * vec4( vColor, opacity );
    float depth = gl_FragCoord.z / gl_FragCoord.w;

    float fogFactor = smoothstep( fogNear, fogFar, depth );

    gl_FragColor = mix( gl_FragColor, vec4( fogColor, gl_FragColor.w ), fogFactor );

}