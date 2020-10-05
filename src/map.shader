
#shader vertex

#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;

uniform mat4 mvp;
uniform mat4 view_rotation; // the quaternion attached to what the mouse has done

out vec3 frag_pos;
out vec3 normal_transfer;
out vec4 colour_transfer;

void main() {

   gl_Position = mvp * vec4(position, 1.0);

   frag_pos = position;
   normal_transfer = normalize(normal);
   colour_transfer = colour;
}

#shader fragment

#version 330 core

in vec3 frag_pos;
in vec3 normal_transfer;
in vec4 colour_transfer;

layout(location = 0) out vec4 out_col;

struct LightSource {
   bool is_on;
   bool directional;
   vec3 position;
   vec3 direction_in_molecule_coordinates_space;
   vec4 ambient;
   vec4 diffuse;
   vec4 specular;
   vec4 halfVector;
   vec3 spotDirection;
   float spotExponent;
   float spotCutoff;
   float spotCosCutoff;
   float constantAttenuation;
   float linearAttenuation;
   float quadraticAttenuation;
};

struct Material {
   float shininess;
   float specular_strength;
   vec4 specular;
};
uniform Material material;

uniform LightSource light_sources[2];
uniform vec4 eye_position;
uniform vec4 background_colour;
uniform bool is_perspective_projection;
uniform bool do_depth_fog;
uniform bool do_diffuse_lighting; // really: all lighting vs ambient only (for demo, really)
uniform bool do_cosine_dependent_opacity;
uniform float map_opacity;
uniform bool do_fresnel;
uniform float fresnel_bias;
uniform float fresnel_scale;
uniform float fresnel_power;

float get_fog_amount(float depth_in) {

   if (! is_perspective_projection) {
      return depth_in * depth_in;
   } else {
      // needs tweaking
      float d = depth_in;
      float d4 = d * d * d * d;
      return d4;
   }

}


void main() {

   // get specular_strength from the material
   float specular_strength = material.specular_strength; // 1.5 is very shiny, shiny is good for transparent maps
   vec4 specular_light_colour = material.specular;

   // a light direction of 0,0,1 is good for fresnelly outlining (well, it used to be)

   float fog_amount = get_fog_amount(gl_FragCoord.z);

   for (int i=0; i<2; i++) {
      if (light_sources[i].is_on) {

         vec3 light_dir = normalize(light_sources[i].direction_in_molecule_coordinates_space);
         // light_dir = normalize(vec3(2,1,-5));
         float dp = dot(normal_transfer, light_dir);

         // we can't have specular lights where there is no diffuse light
         if (dp <= 0.0)
            specular_strength = 0.0;

         dp = max(dp, 0.0); // no negative dot products for diffuse for now, also, zero is avoided.

         float m = clamp(gl_FragCoord.z, 0.0f, 1.0f);

         vec3 eye_pos_3 =  eye_position.xyz;

         vec3 view_dir = eye_pos_3 - frag_pos;
         view_dir = normalize(view_dir);

         vec3 norm_2 = normal_transfer;
         norm_2 = normalize(norm_2);
         vec3 reflect_dir = reflect(-light_dir, norm_2);
         float dp_view_reflect = dot(view_dir, reflect_dir);
         dp_view_reflect = max(dp_view_reflect, 0.0);
         // when the exponent is low, the specular_strength needs to be reduced
         // a low exponent means lots of the map is specular (around the edges)
         float shininess = material.shininess;
         // shininess = 55;
         float spec = pow(dp_view_reflect, shininess);
         vec4 col_specular = specular_strength * spec * specular_light_colour;

         vec4 colour_local = colour_transfer;

         vec4 col_1 = colour_local;  // ambient
         float ambient_strength = 0.5; // was 0.2 for apoferritin
         vec4 col_2 = colour_local * dp;
         vec4 col_3 = col_1 * ambient_strength + 0.8 * col_2 + col_specular;

         if (! do_diffuse_lighting)
            col_3 = col_1;

         vec4 col_4 = col_3;

         // if (do_depth_fog)
         // col_4 = mix(col_3, background_colour, fog_amount);

         col_4 = col_3;

         col_4.a = map_opacity;

         if (false) {
            // float a = 1.0 + 2.0 * dp_view_reflect;
            vec3 eye_to_pos_uv = normalize(frag_pos - eye_pos_3);
            float dp_eye = dot(normalize(normal_transfer), eye_to_pos_uv);
            float a = dp_eye;
            a = clamp(dp, 0.0, 1.0);
            a = 0.7 * (1 - a);
            col_4.a += a * a * a;
         }

         out_col += 0.6 * col_4;
      }
   }

   if (do_fresnel) {
      // R can be calculated in the vertex shader
      vec3 eye_to_frag_pos_uv = normalize(frag_pos - eye_position.xyz);
      float dp_eye = dot(normalize(normal_transfer), eye_to_frag_pos_uv);
      // (I.N) should be 1.0 if we are looking staight on to (perpendicular to) the surface.
      float I_dot_N = dp_eye;
      if (I_dot_N <= 0.0) {
         float bias  = 0.01;
         float scale = 0.5;
         float power = 2.2;
         // I_dot_N = clamp(I_dot_N, 0.0, 1.0); // should not be needed.
         float R0 = fresnel_bias + fresnel_scale * pow((1.0 + I_dot_N), fresnel_power);
         float R = clamp(R0, 0.0, 1.0);
         float g = 0.8; // grey
         vec4 fresnel_colour = R * vec4(g,g,g,1.0);
         // fresnel_colour = mix(colour_transfer, fresnel_colour, 0.5);

         // out_col = mix(fresnel_colour, out_col, 0.00125 * R);
         out_col += fresnel_colour;
      }
   }

   out_col = mix(out_col, background_colour, fog_amount);
   
}
