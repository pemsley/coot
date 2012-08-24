/*
     pygl/freetype_font.h: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
*/


#ifndef _CCP4MG_FREETYPE_FONT_
#define _CCP4MG_FREETYPE_FONT_
#include <string>

#ifdef USE_FREETYPE
#include <ft2build.h>
#include FT_FREETYPE_H
#endif

class MGFontInfo;
MGFontInfo LoadFreeTypeFont(const std::string &foundry, const std::string &family,  const std::string &weight,  const std::string &slant,  const std::string &encoding, int size);
MGFontInfo LoadFreeTypeFont(const char *fontpath, int size);
MGFontInfo LoadFreeTypeFont(const char *name);

#ifdef USE_FREETYPE
MGFontInfo LoadFreeTypeFont(FT_Face face, const std::string &foundry, const std::string &family,
		          const std::string &weight,const  std::string &slant, const std::string &encoding,
			  int size);
#endif

int LoadAllFreeTypeFonts();
int FindNextFreeTypeFontSize(const MGFontInfo &finfo, int step);
#endif /*_CCP4MG_FREETYPE_FONT_  */
