
#shader vertex

#version 330 core
// model.shader

layout(location = 0) in vec3 model_rotation_matrix_0;
layout(location = 1) in vec3 model_rotation_matrix_1;
layout(location = 2) in vec3 model_rotation_matrix_2;
layout(location = 3) in vec3 model_translation; // where the atom is in molecular space
layout(location = 4) in vec3 position; // origin-based cylinder
layout(location = 5) in vec3 normal;   // ditto
layout(location = 6) in vec4 colour;

uniform mat4 mvp;
uniform mat4 view_rotation; // the quaternion attached to what the mouse has done

out vec3 frag_pos_transfer;
out vec3 normal_transfer;
out vec4 colour_transfer;

void main() {

   mat3 model_rotation_matrix = mat3(model_rotation_matrix_0, model_rotation_matrix_1, model_rotation_matrix_2);

   vec3 p2 = position * model_rotation_matrix;
   vec4 p3 = vec4(p2 + model_translation, 1.0);
   gl_Position = mvp * p3;

   vec4 n1 = vec4(normal * model_rotation_matrix, 1.0);

   frag_pos_transfer = p3.xyz;
   normal_transfer = n1.xyz;
   colour_transfer = colour;
}

#shader fragment

#version 330 core
// model.shader

in vec3 frag_pos_transfer;
in vec3 normal_transfer;
in vec4 colour_transfer;

layout(location = 0) out vec4 outputColor;

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
uniform LightSource light_sources[2];

struct Material {
   float shininess;
   float specular_strength;
   vec4 specular;
};
uniform Material material;

uniform vec3 eye_position;
uniform vec4 background_colour;
uniform bool is_perspective_projection;
uniform bool do_depth_fog;
uniform bool do_diffuse_lighting; // really: all lighting vs ambient only (for demo, really)


float get_fog_amount(float depth_in) {

   if (! is_perspective_projection) {
      return depth_in;
   } else {
      // needs tweaking
      float d = depth_in;
      float d4 = d * d * d * d;
      return d * d;
   }

}


void main() {

   vec4 bg_col = background_colour;
   vec4 sum_col = vec4(0,0,0,0);

   float specular_strength = material.specular_strength;
   float shininess = material.shininess;
   outputColor = vec4(0,0,0,0);

   // if (gl_FragCoord.z < 0.5) discard; // useful later maybe (in another shader)

   for (int i=0; i<2; i++) {
      if (light_sources[i].is_on) {

         // sum_col += vec4(0.14, 0, 0.14,1);

         vec3 light_dir = light_sources[i].direction_in_molecule_coordinates_space;
         float dp_raw = dot(normal_transfer, light_dir);
         // we can't have specular lights where there is no diffuse light
         float dp = dp_raw;
         if (dp <= 0.0)
            specular_strength = 0.0;

         // dp = clamp(dp, 0.0, 1.0); // no negative dot products for diffuse
         dp = abs(dp); // no black model interiors

         vec4 lsa = light_sources[i].ambient;
         vec4 lsd = light_sources[i].diffuse;
         vec4 ambient  = colour_transfer * lsa *  0.15;
         vec4 diffuse  = colour_transfer * lsd * dp * 0.9;

         // specular
         vec3 eye_pos = eye_position;
         vec3 norm_2 = normalize(normal_transfer); // not needed, I think
         vec3 view_dir = normalize(eye_pos - frag_pos_transfer);
         vec3 reflect_dir = reflect(-light_dir, norm_2);
         reflect_dir = normalize(reflect_dir); // belt and braces
         float dp_view_reflect = dot(view_dir, reflect_dir);
         dp_view_reflect = clamp(dp_view_reflect, 0.0, 1.0);
         float spec = pow(dp_view_reflect, shininess);
         vec4 specular = specular_strength * spec * light_sources[i].specular;

         sum_col += ambient; // + diffuse + specular;
         sum_col += diffuse;

         // vec3 light_to_eye = normalize(eye_pos - 100.0 * light_dir);
         // sum_col = vec4(0.5 * view_dir + vec3(0.5,0.5,0.5), 1.0);
         // sum_col = vec4(0.5 * light_to_eye + vec3(0.5,0.5,0.5), 1.0);
      }
   }

   float fog_amount = 0.0;
   if (do_depth_fog)
      fog_amount = get_fog_amount(gl_FragCoord.z);

   outputColor += mix(sum_col, bg_col, fog_amount);
   outputColor = colour_transfer;
   outputColor = sum_col;

   // If we are doing depth blur, we need to make the forground objects look more
   // like the background colour (*not* using alpha channel), something like this:
   //
   // uniform bool do_depth_fog
   // uniform float focus_blur_z_depth;
   // uniform float focus_blur_strength;
   //
   // if (do_depth_blur) (
   //    float d = 0.5 * (gl_FragCoord.z + 1); // converts to range 0 -> 1 (I think)
   //    float rc_z = focus_blur_z_depth;
   //    float dd = d - rc_z;
   //    float mf = abs(focus_blur_strength * dd);
   //    if (dd < 0.0) {
   //       outputColor = mix(outputColor, bg_col, mf);
   //    }
   // }
   //
   // Fingers crossed that will work well with the depth-of-fields shader pass.
   //
   // Also something similar for the map shader, I guess.

}


