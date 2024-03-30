/*
 *  CXXFortranFile.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Fri Jan 23 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
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
	size_t getFortranData(char *buffer, const size_t itemSize, const size_t nItems, const int FortranDataType);
	int putFortranData(char *buffer, const size_t itemSize, const size_t nItems, const int FortranDataType);
	int bad();

	void prepareSGIInt(int);
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
