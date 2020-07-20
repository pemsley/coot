
#ifndef FT_CHARACTER_HH
#define FT_CHARACTER_HH

// let's try not to need OpenGL headers to compile whith this header

// FreeType
#include <ft2build.h>
#include FT_FREETYPE_H
#include "ft-character.hh"


struct FT_character {
   unsigned int TextureID;  // ID handle of the glyph texture
   glm::ivec2 Size;       // Size of glyph
   glm::ivec2 Bearing;    // Offset from baseline to left/top of glyph
   FT_Pos     Advance;    // Offset to advance to next glyph
};

#endif // FT_CHARACTER_HH
