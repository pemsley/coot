/*
 *  CXXFortranFile.cpp
 *  CXXSurface
 *
 *  Created by Martin Noble on Fri Jan 23 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include "CXXFortranFile.h"
#include <string.h>

void CXXFortranFile::init(){
  fuichar testFuchar;
  testFuchar.uchars[0]= 0;
  testFuchar.uchars[1]= 1;
  testFuchar.uchars[2]= 2;
  testFuchar.uchars[3]= 3;

  //Examine values of members of this union to recognise native byte order etc
	
  if (testFuchar.i == 66051 && fabs (1.0 - (testFuchar.f/9.25572e-41))<0.01){
    isLikeSGI = 1;
    isLikeI386 = 0;
    cout << "This processors native number representations are SGI-like (hence will not be modified)\n";
  }

  else if (testFuchar.i == 50462976 && fabs (1.0 - (testFuchar.f/3.82047e-37))<0.01){
    isLikeSGI = 0;
    isLikeI386 = 1;
    cout << "This processors native number representations are Linux-like (hence will be modified)\n";
  }

  else {
    isLikeSGI = 1;
    isLikeI386 = 0;
    cout << "This processors native number representations are of unknown type (will pretend to be SGI, but not portable)\n";
    cout << "int " << testFuchar.i << endl;
    cout << "float " << testFuchar.f << endl;
  }
}


void CXXFortranFile::prepareSGIInt (int anInteger){
  if (isLikeSGI) return;
  else if (isLikeI386) {
    fuichar internal;
    internal.i = anInteger;
    unsigned char store = internal.uchars[3];
    internal.uchars[3] = internal.uchars[0];
    internal.uchars[0] = store;
    store = internal.uchars[2];
    internal.uchars[2] = internal.uchars[1];
    internal.uchars[1] = store;
  }
  return;
}

void CXXFortranFile::prepareSGIFloat (float &aFloat){
  if (isLikeSGI) return;
  else if (isLikeI386) {
    fuichar internal;
    internal.f = aFloat;
    unsigned char store = internal.uchars[3];
    internal.uchars[3] = internal.uchars[0];
    internal.uchars[0] = store;
    store = internal.uchars[2];
    internal.uchars[2] = internal.uchars[1];
    internal.uchars[1] = store;
    aFloat = internal.f;
  }
  return;
}

void CXXFortranFile::prepareSGIShort (short &anInteger){
  if (isLikeSGI) return;
  else if (isLikeI386) {
    fuichar internal;
    internal.i = anInteger;
    unsigned char store = internal.uchars[3];
    internal.uchars[3] = internal.uchars[0];
    internal.uchars[0] = store;
    store = internal.uchars[2];
    internal.uchars[2] = internal.uchars[1];
    internal.uchars[1] = store;
    anInteger = internal.i;
  }
  return;
}

CXXFortranFile::CXXFortranFile()
{
  init();
}

CXXFortranFile::~CXXFortranFile()
{
  if (!strcmp(mode, "r")){
    if (!(inputStream.bad())) {
      inputStream.close();
    }
  }
  else if (!strcmp(mode, "w")){
    if (!outputStream.bad()){
      outputStream.close();
    }
  }
}

CXXFortranFile::CXXFortranFile(string filePath, const char *modeRequested)
{
  init();
  strcpy(mode ,modeRequested);
  if (!strcmp(mode,"r")){
    inputStream.open(filePath.c_str());
    if (inputStream.bad()) {
      status=CXXFortranFile::OpenError;
    }
    else {
      status=CXXFortranFile::NoError;
    }
  }
  else if (!strcmp(mode,"w")){
    outputStream.open(filePath.c_str());
    if (outputStream.bad()) {
      status=CXXFortranFile::OpenError;
    }
    else {
      status=CXXFortranFile::NoError;
    }
  }
}	

size_t CXXFortranFile::getFortranData(char *buffer, const size_t itemSize, const size_t nItems, const int fortranDataType)
{

  

  size_t recordSize;
    size_t nBytesToRead, nToReadNow;
    int padding;
	
  nBytesToRead = nItems*itemSize;
	
  // Fortran record starts with an integer that is the record length
  inputStream.read ((char *)&recordSize, 4);
  prepareSGIInt(int(recordSize));
	
  //The record may be padded:  calculate how much larger a record is than the required data
  padding = int(recordSize-nBytesToRead);
  if (padding>=0){
    nToReadNow = nBytesToRead;
  }
  else {
    nToReadNow = recordSize;
    status = SplitRecord;
  }
	
  // Read the actual data
  inputStream.read (buffer, nToReadNow);
	
  // Skip over the padding 
  inputStream.seekg(padding,ios::cur);
	
  cout << recordSize << " " << nBytesToRead << endl;
  // Fortran record ends with an integer that says how long the record actually was
  inputStream.read ((char *)&recordSize, 4);
  prepareSGIInt(int(recordSize));
	
  //PostProcessing
  fuichar internal;
  switch (fortranDataType) {
  case FortranStringData:
    buffer[nBytesToRead]='\0';
    break;
  case FortranFloatData:
    for (int i=0; i<nItems; i++){
      for (int j=0; j<4; j++){
	internal.uchars[j] = buffer[4*i+j];
      }
      prepareSGIFloat(internal.f);
      for (int j=0; j<4; j++){
	buffer[4*i+j] = internal.uchars[j];
      }
    }
    break;
  case FortranIntData:
    for (int i=0; i<nItems; i++){
      for (int j=0; j<4; j++){
	internal.uchars[j] = buffer[4*i+j];
      }
      prepareSGIInt(internal.i);
      for (int j=0; j<4; j++){
	buffer[4*i+j] = internal.uchars[j];
      }
    }
    break;
  }
	
  status = inputStream.bad();
	
  return (nToReadNow);
}

int CXXFortranFile::putFortranData(char *buffer, const size_t itemSize, const size_t nItems, const int fortranDataType)
{
  size_t recordSize;
    size_t nBytesToWrite;
    size_t padding;
	
  nBytesToWrite = nItems*itemSize;
  recordSize = (nBytesToWrite >=80 ? nBytesToWrite: 80);
  recordSize = (!(recordSize%4) ? recordSize : recordSize+(4-recordSize%4));
	
  // Fortran record starts with an integer that is the record length
  {
    int rs = int(recordSize);
    prepareSGIInt(rs);
    outputStream.write ((char *)&rs, 4);
  }
  //The record may be padded:  calculate how much larger a record is than the provided data
  padding = recordSize-nBytesToWrite;
	
  //PreProcessing
  fuichar internal;
  switch (fortranDataType) {
  case FortranStringData:
    buffer[nBytesToWrite]='\0';
    break;
  case FortranFloatData:
    for (int i=0; i<nItems; i++){
      for (int j=0; j<4; j++){
	internal.uchars[j] = buffer[4*i+j];
      }
      prepareSGIFloat(internal.f);
      for (int j=0; j<4; j++){
	buffer[4*i+j] = internal.uchars[j];
      }
    }
    break;
  case FortranIntData:
    for (int i=0; i<nItems; i++){
      for (int j=0; j<4; j++){
	internal.uchars[j] = buffer[4*i+j];
      }
      prepareSGIInt(internal.i);
      for (int j=0; j<4; j++){
	buffer[4*i+j] = internal.uchars[j];
      }
    }
    break;
  }
	
  // write the actual data
  outputStream.write (buffer, nBytesToWrite);
	
  // Skip over the padding
  char *paddingBuffer = new char[padding];
  outputStream.write(paddingBuffer, padding);
  delete [] paddingBuffer;
	
  // Fortran record ends with an integer that says how long the record actually was
  outputStream.write ((char *)&recordSize, 4);
	
  status = outputStream.bad();
	
  return (0);
}

int CXXFortranFile::bad()
{
  return status;
}
