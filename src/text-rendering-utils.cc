
#include <iostream>
#include <map>
#include <string>

#include <gtk/gtk.h>
#include <epoxy/gl.h>

// GLEW
// #define GLEW_STATIC
// #include <GL/glew.h>

// GLFW
// #include <GLFW/glfw3.h>

// GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// FreeType
#include <ft2build.h>
#include FT_FREETYPE_H

#include "Shader.hh"

// Properties
// const GLuint WIDTH = 800, HEIGHT = 600;

/// Holds all state information relevant to a character as loaded using FreeType
struct FT_character {
    GLuint TextureID;   // ID handle of the glyph texture
    glm::ivec2 Size;    // Size of glyph
    glm::ivec2 Bearing;  // Offset from baseline to left/top of glyph
    GLuint Advance;    // Horizontal offset to advance to next glyph
};

std::map<GLchar, FT_character> ft_characters;
GLuint VAO, VBO;

void debug_ft_characters() {

   std::map<GLchar, FT_character>::const_iterator it;
   for (it=ft_characters.begin(); it!=ft_characters.end(); it++) {
      std::cout << "debug ft_characters " << it->first << " " << it->second.TextureID << std::endl;
   }

}

void RenderText(Shader &shader, std::string text, GLfloat x, GLfloat y, GLfloat scale, glm::vec3 color);

// The MAIN function, from here we start our application and run the Game loop
int setup_hud_text(int widget_width, int widget_height, Shader &shader) {

   GLenum err = glGetError();
   std::cout << "RenderText start with err " << err << std::endl;

    // Set OpenGL options
    // glEnable(GL_CULL_FACE);
    // glEnable(GL_BLEND);

    shader.init("shaders/hud-text.shader", Shader::Entity_t::HUD_TEXT);
    glm::mat4 projection = glm::ortho(0.0f, static_cast<GLfloat>(widget_width), 0.0f, static_cast<GLfloat>(widget_height));
    shader.Use();
    glUniformMatrix4fv(glGetUniformLocation(shader.get_program_id(), "projection"), 1, GL_FALSE, glm::value_ptr(projection));
    err = glGetError(); if (err) std::cout << "RenderText Aa " << err << std::endl;
    // FreeType
    FT_Library ft;
    // All functions return a value different than 0 whenever an error occurred
    if (FT_Init_FreeType(&ft))
        std::cout << "ERROR::FREETYPE: Could not init FreeType Library" << std::endl;

    // Load font as face
    std::string vera_font_file_name = "fonts/Vera.ttf";
    bool font_loaded = false;
    FT_Face face;
    if (FT_New_Face(ft, vera_font_file_name.c_str(), 0, &face)) {
        std::cout << "ERROR::FREETYPE: Failed to load font" << std::endl;
    } else {
       font_loaded = true;
       // Set size to load glyphs as
       FT_Set_Pixel_Sizes(face, 0, 48);

       // Disable byte-alignment restriction
       glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

       // Load first 128 characters of ASCII set
       for (GLubyte c = 0; c < 128; c++)
       {
           // Load character glyph
           if (FT_Load_Char(face, c, FT_LOAD_RENDER))
           {
               std::cout << "ERROR::FREETYTPE: Failed to load Glyph" << std::endl;
               continue;
           }
           // Generate texture
           GLuint texture;
           glGenTextures(1, &texture);
           glBindTexture(GL_TEXTURE_2D, texture);
           glTexImage2D(
               GL_TEXTURE_2D,
               0,
               GL_RED,
               face->glyph->bitmap.width,
               face->glyph->bitmap.rows,
               0,
               GL_RED,
               GL_UNSIGNED_BYTE,
               face->glyph->bitmap.buffer
           );
           // Set texture options
           glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
           glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
           glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
           glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
           // Now store character for later use
           GLuint face_glyph_advance_x = face->glyph->advance.x;
           // std::cout << "ID handle of the glyph texture " << texture << std::endl;
           FT_character character = {
               texture,
               glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
               glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
               face_glyph_advance_x
           };
           ft_characters.insert(std::pair<GLchar, FT_character>(c, character));
       }
       glBindTexture(GL_TEXTURE_2D, 0);
       // Destroy FreeType once we're finished
       FT_Done_Face(face);
       FT_Done_FreeType(ft);
    }

    if (font_loaded) {
       // Configure VAO/VBO for texture quads
       glGenVertexArrays(1, &VAO);
       glGenBuffers(1, &VBO);
       glBindVertexArray(VAO);
       glBindBuffer(GL_ARRAY_BUFFER, VBO);
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

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    shader.Use();
    err = glGetError(); if (err) std::cout << "RenderText A0 " << err << std::endl;
    glUniform3f(glGetUniformLocation(shader.get_program_id(), "textColour"), color.x, color.y, color.z);
    glActiveTexture(GL_TEXTURE0);
    err = glGetError(); if (err) std::cout << "RenderText A error " << err << std::endl;
    glBindVertexArray(VAO);
    err = glGetError(); if (err) std::cout << "RenderText B error " << err << std::endl;

    // Iterate through all characters
    std::string::const_iterator c;
    for (c = text.begin(); c != text.end(); c++) {
        err = glGetError(); if (err) std::cout << "RenderText loop start for " << *c << " " << err << std::endl;
        // const FT_character &ch = ft_characters[*c];
        std::map<GLchar, FT_character>::const_iterator it = ft_characters.find(*c);
        if (it == ft_characters.end()) {
           std::cout << "Failed to lookup for " << *c << std::endl;
           continue;
        };
        const FT_character &ch = it->second;

        GLfloat xpos = x + ch.Bearing.x * scale;
        GLfloat ypos = y - (ch.Size.y - ch.Bearing.y) * scale;

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
        // std::cout << "about to bind to textureid " << ch.TextureID << std::endl;
        glBindTexture(GL_TEXTURE_2D, ch.TextureID);
        err = glGetError(); if (err) std::cout << "RenderText loop C glBindTexture() " << err << std::endl;
        // Update content of VBO memory
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
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
    err = glGetError(); if (err) std::cout << "RenderText end D " << err << std::endl;
    glBindTexture(GL_TEXTURE_2D, 0);
    err = glGetError(); if (err) std::cout << "RenderText end D " << err << std::endl;
}

void draw_hud_text(int widget_width, int widget_height, Shader &shader) {
   // this puts the text top right.
   float x = widget_width  - float(widget_width)/180.0 * 130.0;
   float y = widget_height -  30.0;
   // std::cout << "x " << x << " y " << y << std::endl;
   RenderText(shader, "Welcome to Coot",  x, y, 0.5f, glm::vec3(0.6, 0.7, 0.7f));
}
