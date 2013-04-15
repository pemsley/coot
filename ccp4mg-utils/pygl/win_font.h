/*
     pygl/win_font.h: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2005 University of York, CCLRC

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


#if defined(_WIN32) || defined(__WIN32__)

#ifndef _CCP4MG_WIN_FONT_
#define _CCP4MG_WIN_FONT_
#include <string>
class MGFontInfo;
MGFontInfo LoadWinFont(const std::string &foundry, const std::string &family,  const std::string &weight,  const std::string &slant,  const std::string &encoding, int size);
MGFontInfo LoadWinFont(const char *name);
int LoadAllWinFonts();
int FindNextWinFontSize(const MGFontInfo &finfo, int step);
#endif /*_CCP4MG_WIN_FONT_  */

#endif
