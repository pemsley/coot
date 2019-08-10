
#include <iostream>
#include <string>

#include <glm/ext.hpp>
#include <glm/gtx/string_cast.hpp>  // to_string()

#include <gtk/gtk.h>
#include <epoxy/gl.h>

#include "graphics-info.h"

void RenderText(Shader &shader, std::string text, GLfloat x, GLfloat y, GLfloat scale, glm::vec3 color) {

   GLenum err = glGetError();
   std::cout << "RenderText start with err " << err << std::endl;

   shader.Use();
   err = glGetError(); std::cout << "RenderText A " << err << std::endl;
   glUniform3f(glGetUniformLocation(shader.get_program_id(), "textColour"), color.x, color.y, color.z);
   err = glGetError(); std::cout << "RenderText B " << err << std::endl;
   glActiveTexture(GL_TEXTURE0);
   err = glGetError(); std::cout << "RenderText C " << err << std::endl;
   glBindVertexArray(graphics_info_t::hud_text_vertexarray_id);
   err = glGetError(); std::cout << "RenderText D " << err << std::endl;

   // Iterate through all characters
   std::string::const_iterator c;
   for (c = text.begin(); c != text.end(); c++) {
        std::map<GLchar, FT_character>::const_iterator it = graphics_info_t::ft_characters.find(*c);
       if (it == graphics_info_t::ft_characters.end()) {
          std::cout << "ooooooooooooooooops failure to look up" << *c << std::endl;
          continue;
       } else {
          std::cout << "Found " << *c << std::endl;
       }
       FT_character ch = it->second; // cppy atm

        GLfloat xpos = x + ch.Bearing.x * scale;
        GLfloat ypos = y - (ch.Size.y - ch.Bearing.y) * scale;

        std::cout << "Debug ch: " << ch.Bearing.x << " " << ch.Bearing.y << " " << ch.Size.x
                  << " " << ch.Size.y << " " << ch.TextureID << std::endl;

        GLfloat w = ch.Size.x * scale;
        GLfloat h = ch.Size.y * scale;
        // Update VBO for each character
        GLfloat vertices[6][4] = {
            { xpos,     ypos + h,   0.0, 0.0 },
            { xpos,     ypos,       0.0, 1.0 },
            { xpos + w, ypos,       1.0, 1.0 },

            { xpos,     ypos + h,   0.0, 0.0 },
            { xpos + w, ypos,       1.0, 1.0 },
            { xpos + w, ypos + h,   1.0, 0.0 }
        };
        // Render glyph texture over quad
        glBindTexture(GL_TEXTURE_2D, ch.TextureID);
        err = glGetError(); std::cout << "   glBindTexture " << err << std::endl;
        // Update content of VBO memory
        glBindBuffer(GL_ARRAY_BUFFER, graphics_info_t::hud_text_array_buffer_id);
        err = glGetError(); std::cout << "   glBindBuffer " << err << std::endl;
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices); // Be sure to use glBufferSubData and not glBufferData
        err = glGetError(); std::cout << "   glBufferSubData " << err << std::endl;

        // glBindBuffer(GL_ARRAY_BUFFER, 0);
        // Render quad
        glDrawArrays(GL_TRIANGLES, 0, 6);
        // Now advance cursors for next glyph (note that advance is number of 1/64 pixels)
        x += (ch.Advance >> 6) * scale; // Bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
    }
    glBindVertexArray(0);
    glBindTexture(GL_TEXTURE_2D, 0);
}


void draw_hud_text() {

   GLenum err = glGetError();
   std::cout << "draw_hud_text() start " << err << std::endl;
   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   err = glGetError(); std::cout << "draw_hud_text() post-blend " << err << std::endl;

   glm::vec3 color(0.7, 0.7, 0.4);
   Shader &shader = graphics_info_t::shader_for_hud_text;
   shader.Use();
   std::string text = "Hello - this is Coot";
   GLfloat x = 0.1;
   GLfloat y = 0.1;
   GLfloat scale = 0.1;

   GLuint program_id = shader.get_program_id();
   glm::mat4 hud_projection = glm::ortho(0.0f, 900.0f, 0.0f, 900.0f);
   err = glGetError(); std::cout << "draw_hud_text() pre-glUniformMatrix4fv() " << err << std::endl;

   GLuint loc = glGetUniformLocation(program_id, "projection");
   err = glGetError(); std::cout << "draw_hud_text() post glGetUniformLocation() " << err << std::endl;

   std::cout << "loc:" << loc << std::endl;
   glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(hud_projection));
   err = glGetError(); std::cout << "draw_hud_text() glUniformMatrix4fv() " << err << std::endl;

   RenderText(shader, text, 50.0f, 50.0f, 0.6f, color);
}

