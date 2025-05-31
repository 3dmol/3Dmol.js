    float dotProduct = dot( norm, vLight );
    vec3 light = vec3( max( dotProduct, 0.0 ) );
    color *= light;
#ifdef SHADED
    ivec2 dim = textureSize(shading,0);
    float shadowFactor = texture2D(shading,vec2(gl_FragCoord.x/float(dim.x),gl_FragCoord.y/float(dim.y))).r;
    color *= shadowFactor;
#endif    
    gl_FragColor = vec4(color, opacity*opacity );
    if(fogNear != fogFar) {
        float depth = -qi.z;
        float fogFactor = smoothstep( fogNear, fogFar, depth );
        gl_FragColor = mix( gl_FragColor, vec4( fogColor, gl_FragColor.w ), fogFactor );
    }
}