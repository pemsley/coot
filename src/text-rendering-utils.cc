
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

void debug_ft_characters() {

   std::map<GLchar, FT_character>::const_iterator it;
   std::map<GLchar, FT_character> &ft_characters = graphics_info_t::ft_characters;
   for (it=ft_characters.begin(); it!=ft_characters.end(); it++) {
      std::cout << "debug ft_characters " << it->first << " " << it->second.TextureID << std::endl;
   }

}

int setup_hud_text(int widget_width, int widget_height, Shader &shader, bool for_atom_label_flag) {

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

      GLuint projection_uniform_location = shader.atom_label_projection_uniform_location;
      glUniformMatrix4fv(projection_uniform_location, 1, GL_FALSE, glm::value_ptr(projection));

      err = glGetError(); if (err) std::cout << "RenderText Aa " << err << std::endl;
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

    return 0;  // return success status here?

}

void RenderText(Shader &shader, std::string text, GLfloat x, GLfloat y, GLfloat scale, glm::vec3 color)
{

    // Activate corresponding render state
    GLenum err = glGetError(); if (err) std::cout << "RenderText start err " << err << std::endl;

    std::map<GLchar, FT_character> &ft_characters = graphics_info_t::ft_characters;

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    shader.Use();
    err = glGetError(); if (err) std::cout << "RenderText A0 " << err << std::endl;
    glUniform3f(glGetUniformLocation(shader.get_program_id(), "textColour"), color.x, color.y, color.z);
    glActiveTexture(GL_TEXTURE0);
    err = glGetError(); if (err) std::cout << "RenderText A error " << err << std::endl;
    glBindVertexArray(VAO_for_text);
    err = glGetError(); if (err) std::cout << "RenderText B error " << err << std::endl;

    // Iterate through all characters
    std::string::const_iterator c;
    for (c = text.begin(); c != text.end(); c++) {
        err = glGetError(); if (err) std::cout << "RenderText loop start for " << *c << " " << err << std::endl;
        // const FT_character &ch = ft_characters[*c];
        std::map<GLchar, FT_character>::const_iterator it = ft_characters.find(*c);
        if (it == ft_characters.end()) {
           std::cout << "RenderText() Failed to lookup for " << *c << std::endl;
           continue;
        };
        const FT_character &ch = it->second;

        GLfloat xpos = x + ch.Bearing.x * scale;
        GLfloat ypos = y - (ch.Size.y - ch.Bearing.y) * scale;

        GLfloat w = ch.Size.x * scale;
        GLfloat h = ch.Size.y * scale;
        // Update VBO for each character, 2 triangles
        GLfloat vertices[6][4] = {
            { xpos,     ypos + h,   0.0, 0.0 },
            { xpos,     ypos,       0.0, 1.0 },
            { xpos + w, ypos,       1.0, 1.0 },

            { xpos,     ypos + h,   0.0, 0.0 },
            { xpos + w, ypos,       1.0, 1.0 },
            { xpos + w, ypos + h,   1.0, 0.0 }
        };
        // Render glyph texture over quad
        // std::cout << "about to bind to textureid " << ch.TextureID << std::endl;
        glBindTexture(GL_TEXTURE_2D, ch.TextureID);
        err = glGetError(); if (err) std::cout << "RenderText loop C glBindTexture() " << err << std::endl;
        // Update content of VBO memory
        glBindBuffer(GL_ARRAY_BUFFER, VBO_for_text);
        err = glGetError(); if (err) std::cout << "RenderText loop C glBindBuffer() " << err << std::endl;
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices); // Be sure to use glBufferSubData and not glBufferData

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        err = glGetError(); if (err) std::cout << "RenderText loop C glBindBuffer() again " << err << std::endl;
        // Render quad
        glDrawArrays(GL_TRIANGLES, 0, 6);
        err = glGetError(); if (err) std::cout << "RenderText loop C glDrawArrays() " << err << std::endl;
        // Now advance cursors for next glyph (note that advance is number of 1/64 pixels)
        x += (ch.Advance >> 6) * scale; // Bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
    }
    glBindVertexArray(0);
    err = glGetError(); if (err) std::cout << "error:: RenderText end D " << err << std::endl;
    glBindTexture(GL_TEXTURE_2D, 0);
    err = glGetError(); if (err) std::cout << "error:: RenderText end D " << err << std::endl;
}

void render_atom_label(Shader &shader, std::string text, glm::vec3 projected_point,
                       GLfloat scale, glm::vec3 color) {

   return;
   std::cout << "---------- render_atom_label() " << text << " " << glm::to_string(projected_point)
             << std::endl;

   std::map<GLchar, FT_character> &ft_characters = graphics_info_t::ft_characters;

   // need to pass widget width and height, I think.

   // Activate corresponding render state
   GLenum err = glGetError(); if (err) std::cout << "render_atom_label start err " << err << std::endl;

   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   shader.Use();
   err = glGetError(); if (err) std::cout << "error:: render_atom_label A0 " << err << std::endl;
   // GLuint loc = glGetUniformLocation(shader.get_program_id(), "textColour");
   // GLuint loc = shader.atom_label_textColour_uniform_location;

   glm::vec3 atom_colour(color.x, color.y, color.z);
   shader.set_vec3_for_uniform("textColour", atom_colour);
   err = glGetError(); if (err) std::cout << "error:: render_atom_label A1 " << err << std::endl;
   glActiveTexture(GL_TEXTURE0);
   err = glGetError(); if (err) std::cout << "error:: render_atom_label A3 " << err << std::endl;
   glBindVertexArray(VAO_for_text);
   err = glGetError(); if (err) std::cout << "error:: render_atom_label B " << err << std::endl;

   // projection = glm::ortho(0.0f, static_cast<GLfloat>(900),
   // 0.0f, static_cast<GLfloat>(900));
   // GLuint projection_uniform_location = shader.atom_label_projection_uniform_location;
   // glUniformMatrix4fv(projection_uniform_location, 1, GL_FALSE, glm::value_ptr(projection));
   //    err = glGetError(); if (err) std::cout << "RenderText Aa " << err << std::endl;

   float x = projected_point.x;
   float y = projected_point.y;

   // Iterate through all characters
   std::string::const_iterator c;
   for (c = text.begin(); c != text.end(); c++) {
      err = glGetError(); if (err) std::cout << "render_atom_label loop start for "
                                             << *c << " " << err << std::endl;
      // const FT_character &ch = ft_characters[*c];
      std::map<GLchar, FT_character>::const_iterator it = ft_characters.find(*c);
      if (it == ft_characters.end()) {
         std::cout << "render_atom_label(): Failed to lookup for " << *c << std::endl;
         continue;
      };
      const FT_character &ch = it->second;

      GLfloat xpos = x + ch.Bearing.x * scale;
      GLfloat ypos = y - (ch.Size.y - ch.Bearing.y) * scale;

      GLfloat w = ch.Size.x * scale;
      GLfloat h = ch.Size.y * scale;
      // Update VBO for each character, 2 triangles
      GLfloat vertices[6][4] = {
                                { xpos,     ypos + h,   0.0, 0.0 },
                                { xpos,     ypos,       0.0, 1.0 },
                                { xpos + w, ypos,       1.0, 1.0 },

                                { xpos,     ypos + h,   0.0, 0.0 },
                                { xpos + w, ypos,       1.0, 1.0 },
                                { xpos + w, ypos + h,   1.0, 0.0 }
      };
      // Render glyph texture over quad
      // std::cout << "about to bind to textureid " << ch.TextureID << std::endl;
      glBindTexture(GL_TEXTURE_2D, ch.TextureID);
      err = glGetError(); if (err) std::cout << "render_atom_label loop C glBindTexture() " << err << std::endl;
      // Update content of VBO memory
      glBindBuffer(GL_ARRAY_BUFFER, VBO_for_text);
      err = glGetError(); if (err) std::cout << "render_atom_label loop C glBindBuffer() " << err << std::endl;
      glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices); // Be sure to use glBufferSubData and not glBufferData

      glBindBuffer(GL_ARRAY_BUFFER, 0);
      err = glGetError(); if (err) std::cout << "render_atom_label loop C glBindBuffer() again " << err << std::endl;
      // Render quad
      glDrawArrays(GL_TRIANGLES, 0, 6);
      err = glGetError(); if (err) std::cout << "render_atom_label loop C glDrawArrays() " << err << std::endl;
      // Now advance cursors for next glyph (note that advance is number of 1/64 pixels)
      x += (ch.Advance >> 6) * scale; // Bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
   }
   glBindVertexArray(0);
   err = glGetError(); if (err) std::cout << "render_atom_label end D " << err << std::endl;
   glBindTexture(GL_TEXTURE_2D, 0);
   err = glGetError(); if (err) std::cout << "render_atom_label end D " << err << std::endl;
}

void draw_hud_text(int widget_width, int widget_height, Shader &shader) {
   // this puts the text top right.
   float x = widget_width  - float(widget_width)/180.0 * 130.0;
   float y = widget_height -  30.0;
   // std::cout << "x " << x << " y " << y << std::endl;
   x = 3;
   y = 3;
   float scale = 1.0f;
   RenderText(shader, "Welcome to Coot",  x, y, scale, glm::vec3(0.6, 0.7, 0.7f));
}
