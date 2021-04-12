function basic_shader(){
  // Send the resolution of sketch into the shader
  theShader.setUniform('u_resolution', [width, height]);

  // shader() sets the active shader with our shaders
  shader(theShader);

  //rect gives us some geometry on the screen
  rect(0, 0, width, height);

}
