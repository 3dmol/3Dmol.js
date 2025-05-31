uniform mat4 viewMatrix;
uniform float opacity;

uniform vec3 fogColor;
uniform float fogNear;
uniform float fogFar;
#ifdef SHADED
uniform highp sampler2D shading;
#endif
varying vec3 vLightFront;
varying vec3 vColor;
varying vec4 mvPosition;


//DEFINEFRAGCOLOR

void main() {

    gl_FragColor = vec4( vec3 ( 1.0 ), opacity );

    #ifndef WIREFRAME
    gl_FragColor.xyz *= vLightFront;
    #endif
#ifdef SHADED
    ivec2 dim = textureSize(shading,0);
    float shadowFactor = texture2D(shading,vec2(gl_FragCoord.x/float(dim.x),gl_FragCoord.y/float(dim.y))).r;
    vColor *= shadowFactor;
#endif
    gl_FragColor = gl_FragColor * vec4( vColor, opacity );

    if(fogNear != fogFar) {
        float depth = -mvPosition.z;
        float fogFactor = smoothstep( fogNear, fogFar, depth );
        gl_FragColor = mix( gl_FragColor, vec4( fogColor, gl_FragColor.w ), fogFactor );
    }

}