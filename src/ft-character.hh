/*
 * src/ft-character.hh
 *
 * Copyright 2020 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#ifndef FT_CHARACTER_HH
#define FT_CHARACTER_HH

// let's try not to need OpenGL headers to compile whith this header
//  glm is OK
#include <glm/glm.hpp>

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
