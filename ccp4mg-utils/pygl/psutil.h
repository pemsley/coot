/*
     pygl/psutil.h: CCP4MG Molecular Graphics Program
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


#ifndef _CCP4MG_PSUTIL_
#define _CCP4MG_PSUTIL_

#include <fstream>
#include <string>

void PSTail(std::ofstream &fp);
void PSHeader(std::ofstream &fp, const std::string &filename, const std::string &papersize);
void PSHeader(std::ofstream &fp, const std::string &filename);
void zsortps(const std::string &filename, int fog=0);

#endif
