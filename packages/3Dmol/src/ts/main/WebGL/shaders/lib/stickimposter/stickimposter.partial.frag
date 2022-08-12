    float dotProduct = dot( norm, vLight );
    vec3 light = vec3( max( dotProduct, 0.0 ) );
    gl_FragColor = vec4(light*color, opacity*opacity );
    float fogFactor = smoothstep( fogNear, fogFar, depth );
    gl_FragColor = mix( gl_FragColor, vec4( fogColor, gl_FragColor.w ), fogFactor );
}