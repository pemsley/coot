
#ifdef USE_PYTHON
#include <Python.h>
#endif // USE_PYTHON

#include <epoxy/gl.h>

#include "globjects.h"
#include "trackball.h"
#include "graphics-info.h"

#include "draw.hh"

// #define GRAPHICS_TESTING

#ifdef GRAPHICS_TESTING

// int programID_global = -1;
// int location_global = -1;
// GLuint VertexArrayID = -1;

// #define glGenVertexArrays glGenVertexArraysAPPLE
// #define glDeleteVertexArrays glDeleteVertexArraysAPPLE
// #define glBindVertexArray glBindVertexArrayAPPLE

#endif // GRAPHICS_TESTING



void
stereo_projection_setup_maybe(GtkWidget *widget, short int in_stereo_flag) {

   bool do_first  = false;
   bool do_second = false;

   if (graphics_info_t::display_mode_use_secondary_p()) {
      if (widget == graphics_info_t::glarea_2) {
	 do_second = true;
      } else {
	 do_first = true;
      }
   }

   if (in_stereo_flag == IN_STEREO_HARDWARE_STEREO) {
      if (graphics_info_t::which_eye == graphics_info_t::LEFT_EYE)
         do_second = true;
      if (graphics_info_t::which_eye == graphics_info_t::RIGHT_EYE)
         do_first = true;
   }

   if (do_first || do_second) {
      float skew_factor = 0.05 * graphics_info_t::hardware_stereo_angle_factor;
      float view_skew_matrix[16];

      // identity matrices first
      for(unsigned int ii=0; ii<16; ii++) view_skew_matrix[ii]   = 0.0;
      for(unsigned int ii=0; ii<4;  ii++) view_skew_matrix[ii*5] = 1.0;
      float trans_fac = 0.038;

      if (do_first) {
	 view_skew_matrix[8] = skew_factor; // 8 because this is the transpose
	 glMultMatrixf(view_skew_matrix);
	 glTranslatef(trans_fac, 0.0, 0.0);
      } else {
	 view_skew_matrix[8] = -skew_factor;
	 glMultMatrixf(view_skew_matrix);
	 glTranslatef(-trans_fac, 0.0, 0.0);
      }
   }
}


void draw_molecular_triangles(GtkWidget *widget) {

#ifdef USE_MOLECULES_TO_TRIANGLES

      // Martin's triangular molecules
      //
      // centre of the screen
      FCXXCoord pos(graphics_info_t::RotationCentre_x(),
		    graphics_info_t::RotationCentre_y(),
		    graphics_info_t::RotationCentre_z());
      // where is the eye?  That's what we want.
      // front plane is at z=0;
      coot::Cartesian tp_1_cart = unproject_xyz(widget->allocation.width/2,
						widget->allocation.height/2, 1);
      FCXXCoord tp_1(tp_1_cart.x(), tp_1_cart.y(), tp_1_cart.z());
      FCXXCoord diff = tp_1 - pos;
      FCXXCoord eye_pos = pos + diff * 5.0;
      // std::cout << "eye_pos: " << eye_pos << "\n";
      // coot::Cartesian eye_cart = pos + 20 * diff;
      // FCXXCoord eye_pos(eye_cart.x(), eye_cart.y(), eye_cart.z());
      if (graphics_info_t::mol_tri_scene_setup) {
	 if (graphics_info_t::mol_tri_renderer) {

	    //Can retrieve reference to the light if so preferred
	    // This doesn't move the lights
	    // FCXXCoord random_trans(50.0 * coot::util::random()/float(RAND_MAX),
	    // 		              50.0 * coot::util::random()/float(RAND_MAX),
	    //                        50.0 * coot::util::random()/float(RAND_MAX));
	    FCXXCoord light_pos = pos + diff * 10; //  + random_trans;
	    FCXXCoord neg_light_pos = pos + diff * 10; // - random_trans;

	    graphics_info_t::mol_tri_scene_setup->getLight(0)->setTranslation(light_pos);
	    graphics_info_t::mol_tri_scene_setup->getLight(1)->setTranslation(neg_light_pos);

	    for (int ii=graphics_info_t::n_molecules()-1; ii>=0; ii--) {
	       if (graphics_info_t::is_valid_model_molecule(ii)) {
		  if (graphics_info_t::molecules[ii].draw_it) {
		     if (graphics_info_t::molecules[ii].molrepinsts.size()) {
			// molrepinsts get added to mol_tri_scene_setup when then are made
			// turns on glLighting.
			graphics_info_t::mol_tri_scene_setup->renderWithRendererFromViewpoint(graphics_info_t::mol_tri_renderer, eye_pos);
		     }
		  }
	       }
	    }
	    glDisable(GL_LIGHTING);
	 }
      }

#endif // USE_MOLECULES_TO_TRIANGLES

}


void
display_density_level_maybe() {

   if (graphics_info_t::display_density_level_on_screen == 1) {
      if (graphics_info_t::display_density_level_this_image == 1) {

	 //	 std::cout << "DEBUG:: screen label: "
	 // << graphics_info_t::display_density_level_screen_string
	 // << std::endl;

	 bool is_bb = graphics_info_t::background_is_black_p();

	 GLfloat white[3] = {1.0, 1.0, 1.0};
	 GLfloat black[3] = {0.0, 0.0, 0.0};

	 if (is_bb)
	    glColor3fv(white);
	 else
	    glColor3fv(black);

	 glPushMatrix();
	 glLoadIdentity();

 	 glMatrixMode(GL_PROJECTION);
	 glPushMatrix();
 	 glLoadIdentity();

	 // Disable the fog so that the density level text will not
	 // change intensity in a zoom-dependent way:
	 glPushAttrib(GL_ENABLE_BIT);
	 glDisable(GL_FOG);

	 // glRasterPos3f();
	 graphics_info_t::printString_for_density_level(graphics_info_t::display_density_level_screen_string,
							0.0, 0.95, -0.9);

         glPopAttrib();
	 glPopMatrix();
	 glMatrixMode(GL_MODELVIEW);

	 glPopMatrix();

      }
   }
   graphics_info_t::display_density_level_this_image = 0;

   // glScalef (2.0, 2.0, 2.0);
   // glRasterPos3f(2.0, 28.0,0.0);

   // This does resonable things when changing the contour level (it
   // waits until each change has been made before returning (rather
   // than combining them into one change), but terrible things to
   // regular rotation (unusable jerkiness).
//    while (gtk_events_pending())
//       gtk_main_iteration();

}


#include "c-interface-generic-objects.h"

coot::Cartesian eye_position() {

   // eye position is the screen centre rotated by graphics_info_t::quat matrix
   // and translated by a length related to graphics_info_t::zoom

   coot::Cartesian rc(graphics_info_t::RotationCentre_x(),
		      graphics_info_t::RotationCentre_y(),
		      graphics_info_t::RotationCentre_z());

   float dist = 0.5 * graphics_info_t::zoom;

   GL_matrix glm;
   clipper::Coord_orth eye_dir(0,0,1);
   glm.from_quaternion(graphics_info_t::quat);
   clipper::Mat33<double> m = glm.to_clipper_mat();

   clipper::Coord_orth rot_dir(m * eye_dir);
   coot::Cartesian rot_dir_c(rot_dir.x(), rot_dir.y(), rot_dir.z());

   coot::Cartesian eye_position(rc + rot_dir_c * dist);

   return eye_position;
}

void
debug_eye_position(GtkWidget *widget) {

   coot::Cartesian rc(graphics_info_t::RotationCentre_x(),
		      graphics_info_t::RotationCentre_y(),
		      graphics_info_t::RotationCentre_z());

   coot::Cartesian ep = eye_position();

   coot::Cartesian pt((ep + rc) * 0.5);

   int go = generic_object_index("eye position");
   if (go == -1)
      go = new_generic_object_number("eye position");

   to_generic_object_add_point(go, "red", 4, pt.x(), pt.y(), pt.z());
   set_display_generic_object(go, 1);
}

shader_program_source
parse_shader(const std::string &file_name) {

   enum class ShaderType { NONE = -1, VERTEX = 0, FRAGMENT = 1 };

   ShaderType type = ShaderType::NONE;
   shader_program_source ss;
   std::ifstream f(file_name.c_str());
   if (f) {
      std::string line;
      while(std::getline(f, line)) {
	 if (line.find("#shader") != std::string::npos) {
	    if (line.find("vertex") != std::string::npos)
	       type = ShaderType::VERTEX;
	    if (line.find("fragment") != std::string::npos)
	       type = ShaderType::FRAGMENT;
	 } else {
	    if (type == ShaderType::VERTEX)
	       ss.VertexSource += line + "\n";
	    if (type == ShaderType::FRAGMENT)
	       ss.FragmentSource += line + "\n";
	 }
      }
   } else {
      std::cout << "Failed to open " << file_name  << std::endl;
   }
   return ss;
}

unsigned int compile_shader(const std::string &source, unsigned int type) {

#ifdef GRAPHICS_TESTING
   std::string type_s = "vertex";
   if (type == GL_FRAGMENT_SHADER)
      type_s = "fragment";
   unsigned int id = glCreateShader(type);
   const char *s = source.c_str();
   int l = source.size() + 1;
   glShaderSource(id,  1,  &s, &l);
   glCompileShader(id);

   int result;
   glGetShaderiv(id, GL_COMPILE_STATUS, &result);
   if (result == GL_FALSE) {
      int length;
      glGetShaderiv(id, GL_INFO_LOG_LENGTH, &length);
      char message[length+1];
      glGetShaderInfoLog(id, length, &length, message);
      std::cout << "Failed to compile " << type_s << " shader: " << message << std::endl;
   } else {
      std::cout << "glCompileShader() result was good for " << type_s << " shader " << std::endl;
   }
   return id;
#else
   std::cout << "compile_shader() return non-graphics testing 0 " << std::endl;
   return 0;
#endif
}

std::string file_to_string(const std::string &file_name) {

   std::ifstream f(file_name.c_str());
   if (f) {
      std::string s((std::istreambuf_iterator<char>(f)),
		    std::istreambuf_iterator<char>());
      return s;
   } else {
      return std::string("");
   }
}

unsigned int CreateShader(const std::string &vertex_shader, const std::string &fragment_shader) {

   unsigned int program  = glCreateProgram();
   unsigned int vs = compile_shader(vertex_shader, GL_VERTEX_SHADER);
   unsigned int fs = compile_shader(fragment_shader, GL_FRAGMENT_SHADER);

   glAttachShader(program, vs);
   glAttachShader(program, fs);
   glLinkProgram(program);
   glValidateProgram(program);

   glDeleteShader(vs);
   glDeleteShader(fs);

   return program;

}

void setup_for_single_triangle() {


   std::cout << "---------- setup for single triangle: " << glGetString(GL_VERSION) << std::endl;

   {
      float positions[12] = {
	 -0.5,  -0.5, -0.0,
   	 -0.5,   0.5, -0.0,
  	  0.5,   0.5, -0.0,
	  0.5,  -0.5, -0.0
      };

      unsigned int indices[8] { 0,1,1,2,2,3,3,0 };

      // GLuint VertexArrayID;
      glGenVertexArrays(1, &VertexArrayID);
      glBindVertexArray(VertexArrayID);

      GLuint vertexbuffer;
      glGenBuffers(1, &vertexbuffer);
      glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
      glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 12, &positions[0], GL_STATIC_DRAW);
      glEnableVertexAttribArray(0);
      glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

      unsigned int ibo;
      glGenBuffers(1, &ibo);
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * 8, &indices[0], GL_STATIC_DRAW);

      std::cout << "----------- parse and create shader " << std::endl;
      shader_program_source sps = parse_shader("Basic.shader");
      unsigned int programID = CreateShader(sps.VertexSource, sps.FragmentSource);
      programID_global = programID;
      std::cout << "----------- created shader program " << programID << std::endl;

      glBindAttribLocation(programID, 0, "position");

      // int ul = glGetUniformLocation(programID, "u_Color");
      // std::cout << "Got glGetUniformLocation for u_Color " << ul << std::endl;
      // location_global = ul;

   }
}

void draw_single_triangle() {

#ifdef GRAPHICS_TESTING

   glBindVertexArray(VertexArrayID);
   glUseProgram(programID_global);
   glDrawElements(GL_LINES, 8, GL_UNSIGNED_INT, nullptr);
   glBindVertexArray(0); // unbind
   glUseProgram(0);

#endif // GRAPHICS_TESTING
}

/* Basic.shader

#shader vertex

#version 120

uniform mat4 mygl_ModelViewMatrix;
uniform mat4 mygl_ProjectionMatrix;
uniform mat3 mygl_NormalMatrix;

attribute vec3 position;

void main() {

   gl_Position = gl_ModelViewProjectionMatrix * vec4(position, 1.0);

}

#shader fragment

#version 120

void main() {

  gl_FragColor = vec4(0.6, 0.1, 0.4, 1.0);

}
*/
