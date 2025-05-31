

uniform vec3 color;
uniform sampler2D map;
uniform float opacity;

uniform vec3 fogColor;
uniform float fogNear;
uniform float fogFar;
uniform float alphaTest;

varying vec2 vUV;
//DEFINEFRAGCOLOR

void main() {

    vec4 texture = texture2D(map, vUV);

    if (texture.a <= alphaTest) discard;

    gl_FragColor = vec4(color * texture.xyz, texture.a * opacity);

    if (fogNear != fogFar) {

        float depth = gl_FragCoord.z / gl_FragCoord.w; //probably wrong
        float fogFactor = smoothstep(fogNear, fogFar, depth);        
        gl_FragColor = mix(gl_FragColor, vec4(fogColor, gl_FragColor.w), fogFactor);
    }
}

