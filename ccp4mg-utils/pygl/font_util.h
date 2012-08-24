/*
     pygl/font_util.h: CCP4MG Molecular Graphics Program
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


#ifndef _CCP4MG_FONT_UTIL_
#define _CCP4MG_FONT_UTIL_
#include "font_info.h"
#include <ppmutil.h>
#include <string>

image_info MakeExampleStringImage(const MGFontInfo &finfo, const std::string &fg="#ffffff", const std::string &bg="#000000", const bool &scale_alpha=false);
image_info MakeStringImage(const MGFontInfo &finfo, const std::string &str, const std::string &fg="#ffffff", const std::string &bg="#000000", const bool &scale_alpha=false);
int WriteStringImage(const MGFontInfo &finfo, const std::string &str, const std::string &filename, const std::string &fg="#ffffff", const std::string &bg="#000000", const bool &scale_alpha=false);
image_info OverlayTextOnImage(const image_info &iinfo, const MGFontInfo &finfo, int x, int y, const std::string &str, const std::string &fg="#000000", const std::string &bg="#ffffff", const bool &scale_alpha=false);
#endif /* _CCP4MG_FONT_UTIL_ */
