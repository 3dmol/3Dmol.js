
uniform highp sampler3D data;
uniform highp sampler2D colormap;
uniform highp sampler2D depthmap;


uniform mat4 textmat;
uniform mat4 projinv;
uniform mat4 projectionMatrix;

uniform float step;
uniform float subsamples;
uniform float maxdepth;
uniform float transfermin;
uniform float transfermax;
in  vec4 mvPosition;
out vec4 color;
void main(void) {

   vec4 pos = mvPosition;
   bool seengood = false;
   float i = 0.0;
   color = vec4(1,1,1,0);
   float increment = 1.0/subsamples;
   float maxsteps = (maxdepth*subsamples/step);
//there's probably a better way to do this..
//calculate farthest possible point in model coordinates
   vec4 maxpos = vec4(pos.x,pos.y,pos.z-maxdepth,1.0);
// convert to projection
   maxpos = projectionMatrix*maxpos;
   vec4 startp = projectionMatrix*pos;
// homogonize
   maxpos /= maxpos.w;
   startp /= startp.w;
//take x,y from start and z from max
   maxpos = vec4(startp.x,startp.y,maxpos.z,1.0);
//convert back to model space
   maxpos = projinv*maxpos;
   maxpos /= maxpos.w;
   float incr = step/subsamples;
//get depth from depthmap
//startp is apparently [-1,1]
   vec2 tpos = startp.xy/2.0+0.5;
   float depth = texture(depthmap, tpos).r;
//compute vector between start and end
   vec4 direction = maxpos-pos;
   for( i = 0.0; i <= maxsteps; i++) {
      vec4 pt = (pos+(i/maxsteps)*direction);
      vec4 ppt = projectionMatrix*pt;
      float ptdepth = ppt.z/ppt.w;
      ptdepth = ((gl_DepthRange.diff * ptdepth) + gl_DepthRange.near + gl_DepthRange.far) / 2.0;
      if(ptdepth > depth) break;
      pt = textmat*pt;
//       pt /= pt.w;
      if(pt.x >= -0.01 && pt.y >= -0.01 && pt.z >= -0.01 && pt.x <= 1.01 && pt.y <= 1.01 && pt.z <= 1.01) {
         seengood = true;
      } else if(seengood) {
         break;
      }
      if( pt.x < -0.01 || pt.x > 1.01 || pt.y < -0.01 || pt.y > 1.01 || pt.z < -0.01 || pt.z > 1.01  ){
          color.a = 0.0;
          continue;
      }
      else {
         float val = texture(data, pt.zyx).r;
         if(isinf(val)) continue; //masked out
         float cval = (val-transfermin)/(transfermax-transfermin); //scale to texture 0-1 range
         vec4 val_color = texture(colormap, vec2(cval,0.5));
         color.rgb = color.rgb*color.a + (1.0-color.a)*val_color.a*val_color.rgb;
         color.a += (1.0 - color.a) * val_color.a; 
         if(color.a > 0.0) color.rgb /= color.a;
//          color = vec4(pt.x, pt.y, pt.z, 1.0);
      }
//       color = vec4(pt.x, pt.y, pt.z, 0.0)
    }
}

        