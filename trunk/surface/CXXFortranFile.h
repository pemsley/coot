/* 
 * 
 * Copyright 2004 by The University of Oxford
 * Author: Martin Noble, Jan Gruber
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */
/*
 *  CXXFortranFile.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Fri Jan 23 2004.
 *  
 *
 */

#ifndef CXXFortranFile_included
#define CXXFortranFile_included

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

using namespace std;

union fuichar {
  char uchars[4];
  int i;
  float f;
};


class CXXFortranFile{
	
private:
	ifstream inputStream;
	ofstream outputStream;
	int status;
	char mode[32];
	int isLikeSGI, isLikeI386;
	void init();
public:
	CXXFortranFile();
	CXXFortranFile(string filePath, const char *mode);
	~CXXFortranFile();
	int getFortranData(char *buffer, const int itemSize, const int nItems, const int FortranDataType);
	int putFortranData(char *buffer, const int itemSize, const int nItems, const int FortranDataType);
	int bad();

	void prepareSGIInt(int&);
	void prepareSGIFloat(float&);
	void prepareSGIShort(short&);

	
	static const int NoError = 0;
	static const int OpenError = 1;
	static const int ReadError = 2;
	static const int SplitRecord = 3;
	
	static const int FortranStringData = 0;
	static const int FortranCharData = 1;
	static const int FortranShortData = 2;
	static const int FortranIntData = 3;
	static const int FortranFloatData = 4;
};

#endif
