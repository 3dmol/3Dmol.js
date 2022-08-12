

uniform vec3 color;
uniform sampler2D map;
uniform float opacity;

uniform int fogType;
uniform vec3 fogColor;
uniform float fogDensity;
uniform float fogNear;
uniform float fogFar;
uniform float alphaTest;

varying vec2 vUV;
//DEFINEFRAGCOLOR

void main() {

    vec4 texture = texture2D(map, vUV);

    if (texture.a < alphaTest) discard;

    gl_FragColor = vec4(color * texture.xyz, texture.a * opacity);

    if (fogType > 0) {

        float depth = gl_FragCoord.z / gl_FragCoord.w;
        float fogFactor = 0.0;

        if (fogType == 1) {
            fogFactor = smoothstep(fogNear, fogFar, depth);
        }

        else {
            const float LOG2 = 1.442695;
            float fogFactor = exp2(- fogDensity * fogDensity * depth * depth * LOG2);
            fogFactor = 1.0 - clamp(fogFactor, 0.0, 1.0);
        }

        gl_FragColor = mix(gl_FragColor, vec4(fogColor, gl_FragColor.w), fogFactor);

    }
}

