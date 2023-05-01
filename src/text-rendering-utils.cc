
#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <iostream>
#include <map>
#include <string>

#include <gtk/gtk.h>
#include <epoxy/gl.h>

// GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/ext.hpp>

// FreeType
#include <ft2build.h>
#include FT_FREETYPE_H

#ifdef THIS_IS_HMT
#else
#include "graphics-info.h"
#endif

#include "Shader.hh"

// Properties
// const GLuint WIDTH = 800, HEIGHT = 600;

#include "ft-character.hh"

// these are shared for HUD text and atom labels - I am not sure that that's a good idea.
GLuint VAO_for_text, VBO_for_text;


int setup_hud_text(int widget_width, int widget_height, Shader &shader, bool for_atom_label_flag) {

   // std::cout << ":::::::::::::::::::: setup hud_text()" << std::endl;

   GLenum err = glGetError();
   // std::cout << "RenderText start with err " << err << std::endl;

   // Set OpenGL options
   // glEnable(GL_CULL_FACE);
   // glEnable(GL_BLEND);

   // shader.init("hud-text.shader", Shader::Entity_t::HUD_TEXT);

   glm::mat4 projection = glm::ortho(0.0f, static_cast<GLfloat>(widget_width),
                                     0.0f, static_cast<GLfloat>(widget_height));
   shader.Use();
   if (for_atom_label_flag) {

      // GLuint projection_uniform_location = shader.atom_label_projection_uniform_location;
      GLuint projection_uniform_location = shader.mvp_uniform_location;
      glUniformMatrix4fv(projection_uniform_location, 1, GL_FALSE, glm::value_ptr(projection));

      err = glGetError();
      if (err) std::cout << "error in setup_hud_text() RenderText Aa " << err << std::endl;

   } else {
      GLuint projection_uniform_location = shader.hud_projection_uniform_location;
      glUniformMatrix4fv(projection_uniform_location, 1, GL_FALSE, glm::value_ptr(projection));
      err = glGetError(); if (err) std::cout << "RenderText Ab " << err << std::endl;
   }

   graphics_info_t g;
   g.load_freetype_font_textures();

   if (g.vera_font_loaded) {
      // Configure VAO/VBO for texture quads
      glGenVertexArrays(1, &VAO_for_text);
      glGenBuffers(1, &VBO_for_text);
      glBindVertexArray(VAO_for_text);
      glBindBuffer(GL_ARRAY_BUFFER, VBO_for_text);
      glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 6 * 4, NULL, GL_DYNAMIC_DRAW);
      glEnableVertexAttribArray(0);
      glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), 0);
      glBindBuffer(GL_ARRAY_BUFFER, 0);
      glBindVertexArray(0);
   }

    // debug_ft_characters();

   // std::cout << ":::::::::::::::::::: setup hud_text() done" << std::endl;
   return 0;  // return success status here?

}

