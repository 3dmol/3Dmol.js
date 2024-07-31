uniform sampler2D depthmap;
varying highp vec2 vTexCoords;
uniform float vWidth;
uniform float vHeight;
uniform float total_strength;
uniform float radius;
uniform mat4 projinv;
uniform mat4 projectionMatrix;

//code for pseudorandom vector from chatgpt
float hash(vec3 p) {
    p = fract(p * vec3(0.1031, 0.1030, 0.0973));
    p += dot(p, p.yzx + 33.33);
    return fract((p.x + p.y) * p.z);
}

// Generate a pseudorandom vec3 from a seed vec3
vec3 pseudorandomVec3(vec3 seed) {
    vec3 randomVec;
    randomVec.x = hash(seed);
    randomVec.y = hash(seed + vec3(1.0, 0.0, 17.1));
    randomVec.z = hash(seed + vec3(0.0, 13.23, 0.0));
    return randomVec;
}

void main(void)
{   
   const float base = 0.2;
   const float area = 0.75;
   const int cycles = 1;

   const int samples = 64;
   vec3 sample_sphere[64] = vec3[](
      vec3(0.091258,-0.510164,0.000000),
      vec3(-0.204347,-0.872967,0.187199),
      vec3(0.009690,-0.263696,-0.110414),
      vec3(0.175208,-0.563987,0.228527),
      vec3(-0.001824,-0.003113,-0.000323),
      vec3(0.411134,-0.719869,-0.261530),
      vec3(-0.074272,-0.377368,0.276290),
      vec3(-0.147773,-0.381587,-0.284529),
      vec3(0.173317,-0.199635,0.063295),
      vec3(-0.186452,-0.199460,0.076965),
      vec3(0.143985,-0.308160,-0.307687),
      vec3(0.053194,-0.148286,0.169589),
      vec3(-0.547656,-0.486476,-0.317378),
      vec3(0.020804,-0.015092,-0.004574),
      vec3(-0.038006,-0.043165,0.054059),
      vec3(-0.094795,-0.443908,-0.731525),
      vec3(0.547552,-0.396466,0.461477),
      vec3(-0.176886,-0.089989,0.007315),
      vec3(0.074401,-0.048840,-0.074039),
      vec3(-0.008240,-0.075697,0.178197),
      vec3(-0.307880,-0.185053,-0.368943),
      vec3(0.309520,-0.108483,0.041646),
      vec3(-0.773478,-0.292946,0.538166),
      vec3(0.184487,-0.231594,-0.820065),
      vec3(0.207318,-0.100531,0.361797),
      vec3(-0.173306,-0.037737,-0.055289),
      vec3(0.548102,-0.105342,-0.253237),
      vec3(-0.119342,-0.043907,0.285162),
      vec3(-0.270247,-0.087861,-0.751357),
      vec3(0.449312,-0.039777,0.236146),
      vec3(-0.743773,-0.036095,0.196056),
      vec3(0.148819,-0.004300,-0.231448),
      vec3(0.008773,0.000809,0.051047),
      vec3(-0.461467,0.027390,-0.357386),
      vec3(0.169626,0.013338,-0.014053),
      vec3(-0.043786,0.007095,0.047331),
      vec3(0.004821,0.140371,-0.988260),
      vec3(0.092402,0.023994,0.101860),
      vec3(-0.295335,0.061530,-0.027372),
      vec3(0.024903,0.007537,-0.018901),
      vec3(-0.081463,0.125402,0.447794),
      vec3(-0.397119,0.231805,-0.631062),
      vec3(0.163853,0.059014,0.044905),
      vec3(-0.495220,0.214357,0.254139),
      vec3(0.306123,0.373687,-0.825715),
      vec3(0.021665,0.026737,0.053220),
      vec3(-0.208231,0.117129,-0.098685),
      vec3(0.139749,0.080968,-0.043086),
      vec3(-0.153599,0.182814,0.262090),
      vec3(-0.159673,0.496777,-0.743568),
      vec3(0.134797,0.117152,0.095753),
      vec3(-0.155626,0.120533,0.019395),
      vec3(0.042311,0.054462,-0.049709),
      vec3(0.001257,0.031288,0.034468),
      vec3(-0.002271,0.003199,-0.002304),
      vec3(0.662104,0.717307,0.033854),
      vec3(-0.373100,0.576021,0.308274),
      vec3(0.024233,0.231316,-0.173688),
      vec3(0.161311,0.420217,0.234273),
      vec3(-0.045248,0.078031,-0.010411),
      vec3(0.167453,0.376942,-0.094872),
      vec3(-0.056194,0.433247,0.173218),
      vec3(-0.016224,0.123149,-0.035569),
      vec3(0.067127,0.407641,0.028479)
   );

   float depth = texture2D(depthmap, vTexCoords).r;
   if(depth == 1.0) {
      discard;
   }

   vec4 screen_position = vec4(vTexCoords, depth,1.0);
   vec4 pos = projinv*screen_position;
   pos /= pos.w;
   vec3 position = pos.xyz;

   //approximate the normal from the depth map; I find this simpler
   //than trying to recompute the exact normal within every possible object shader

   //pixel offset positions in screen space
   ivec2 dim = textureSize(depthmap,0);
   vec2 offset1 = vec2(0.0,1.0/float(dim.y));
   vec2 offset2 = vec2(1.0/float(dim.x),0.0);
   float depth1 = texture2D(depthmap, vTexCoords + offset1).r;
   float depth2 = texture2D(depthmap, vTexCoords + offset2).r;
   
   vec3 p1 = vec3(screen_position.xy+offset1, depth1 - depth);
   vec3 p2 = vec3(screen_position.xy+offset2, depth2 - depth);

   //convert to model space
   vec4 pos1 = projinv*vec4(p1,1);
   pos1 /= pos1.w;
   vec4 pos2 = projinv*vec4(p2,1);
   pos2 /= pos2.w;

   vec3 normal = normalize(cross(pos1.xyz-position, pos2.xyz-position)); //model normal, important we normalize in model space

   //pseudo random number

   float occlusion = 0.0;
   for(int c = 0; c < cycles; c++) {
   vec3 random = normalize(pseudorandomVec3(position+float(c)));
   for(int i=0; i < samples; i++) {

      vec3 ray = radius * reflect(sample_sphere[i],random);
      vec3 hemi_ray = position + sign(dot(ray,normal)) * ray; //model space
      vec4 hemi_screen = projectionMatrix*vec4(hemi_ray,1.0);
      hemi_screen /= hemi_screen.w;
      
      float occ_depth = texture2D(depthmap, clamp(hemi_screen.xy,0.0,1.0)).r;
      float difference = hemi_screen.z - occ_depth;
      
      occlusion += step(0.0, difference) * (1.0-smoothstep(0.0, area, difference));
   }
   }
   float ao = 1.0 - total_strength * occlusion * (1.0 / float(cycles*samples));
   gl_FragDepthEXT = clamp(ao+base,0.0,1.0);

}