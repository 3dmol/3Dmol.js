
    gl_FragColor = vec4(color, opacity*opacity );
    if(fogNear != fogFar) {
        float depth = -qi.z;
        float fogFactor = smoothstep( fogNear, fogFar, depth );
        gl_FragColor = mix( gl_FragColor, vec4( fogColor, gl_FragColor.w ), fogFactor );
    }
}