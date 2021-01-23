//   CCP4 Molecular Graphics Program
//
//   Copyright 2004 The University of York
//   Author: Stuart McNicholas and Liz Potterton
//
//   This program is free software and is distributed under the terms
//   and conditions of the CCP4 licence agreement as `Part 0' (Annex 2)
//   software, which is version 2.1 of the GNU Lesser General Public
//   Licence (LGPL) with the following additional clause:
//
//      `You may also combine or link a "work that uses the Library"
//      to produce a work containing portions of the Library, and
//      distribute that work under terms of your choice, provided that
//      you give prominent notice with each copy of the work that the
//      specified version of the Library is used in it, and that you
//      include or provide public access to the complete corresponding
//      machine-readable source code for the Library including whatever
//      changes were used in the work. (i.e. If you make changes to the
//      Library you must distribute those, but you do not need to
//      distribute source or object code to those portions of the work
//      not covered by this licence.)'
//
//   Note that this clause grants an additional right and does not
//   impose any additional restriction, and so does not affect
//   compatibility with the GNU General Public Licence (GPL). If you
//   wish to negotiate other terms, please contact the maintainer.
//   You can redistribute it and/or modify the program under the terms
//   of the GNU Lesser General Public License as published by the Free
//   Software Foundation; either version 2.1 of the License, or (at
//   your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//   Lesser General Public License for more details.
//
//   You should have received a copy of the CCP4 licence and/or GNU
//   Lesser General Public License along with this program; if not,
//   write to the CCP4 Secretary, Daresbury Laboratory, Warrington
//   WA4 4AD, UK. The GNU Lesser General Public can also be obtained
//   writing to the Free Software Foundation, Inc., 51 Franklin
//   Street, Fifth Floor, Boston, MA 02110-1301, USA


#ifdef USE_LIBGIF
#ifdef __cplusplus      //include the gif library  header as a C file.
extern "C"{
#endif
#include <gif_lib.h>
#ifdef __cplusplus
}
#endif
#endif

#ifdef USE_LIBPNG
#include <png.h>
#endif

#ifdef USE_LIBXPM
#include <X11/xpm.h>
#endif

#ifdef USE_LIBTIFF
#include <tiffio.h>
#endif

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <limits.h>
#include <iostream>
#include "ppmutil.h"
#include <math.h>
#include <setjmp.h>
#include <map>
#include <algorithm>

#if defined(linux) || defined (__APPLE_CC__) || defined (__MINGW32__)
#include <stdint.h>
#endif
#if defined(sgi) || defined(sun)
#include <sys/types.h>
#endif
#if defined(__osf__)
#include <inttypes.h>
#endif

#ifndef int16_t
typedef signed short    int16_t;
#endif
#ifndef uint16_t
typedef unsigned short  uint16_t;
#endif
#ifndef int32_t
typedef signed int      int32_t;
#endif
#ifndef uint32_t
typedef unsigned int    uint32_t;
#endif
#ifdef WIN32
#define strcasecmp _stricmp
#endif

#ifdef _USE_DL_
#include "ppmdl.h"
#if defined (linux)
#define LIBPNG_SOLIBRARY "libpng.so"
#define LIBGIF_SOLIBRARY "libungif.so"
#define LIBJPG_SOLIBRARY "libjpeg.so"
#define LIBTIFF_SOLIBRARY "libtiff.so"
#define LIBXPM_SOLIBRARY "libXpm.so"
#elif defined (__APPLE_CC__)
#define LIBPNG_SOLIBRARY "libpng.dylib"
#define LIBGIF_SOLIBRARY "libungif.dylib"
#define LIBJPG_SOLIBRARY "libjpeg.dylib"
#define LIBTIFF_SOLIBRARY "libtiff.dylib"
#define LIBXPM_SOLIBRARY "libXpm.dylib"
#elif  defined(_WIN32) || defined(__WIN32__)
#define LIBPNG_SOLIBRARY "libpng.dll"
#define LIBGIF_SOLIBRARY "libungif.dll"
#define LIBJPG_SOLIBRARY "jpeg-62.dll"
#define LIBTIFF_SOLIBRARY "libtiff.dll"
#define LIBXPM_SOLIBRARY "libXpm.dll" // Unlikely to exist.
#else
#define LIBPNG_SOLIBRARY "libpng.so"
#define LIBGIF_SOLIBRARY "libungif.so"
#define LIBJPG_SOLIBRARY "libjpeg.so"
#define LIBTIFF_SOLIBRARY "libtiff.so"
#define LIBXPM_SOLIBRARY "libXpm.so"
#endif
#endif

class ImageInfoFReadExc {};
class ImageInfoFWriteExc {};
class ImageInfoReadPPMExc {};
class ImageInfoWritePPMExc {};
class ImageInfoReadJPGExc {};
class ImageInfoWriteJPGExc {};
class ImageInfoReadBMPExc {};
class ImageInfoWriteBMPExc {};
class ImageInfoReadPNGExc {};
class ImageInfoWritePNGExc {};
class ImageInfoReadYUVExc {};
class ImageInfoWriteYUVExc {};
class ImageInfoReadRGBExc {};
class ImageInfoWriteRGBExc {};
class ImageInfoReadGIFExc {};
class ImageInfoWriteGIFExc {};
class ImageInfoReadXBMExc {};
class ImageInfoWriteXBMExc {};
class ImageInfoReadXPMExc {};
class ImageInfoWriteXPMExc {};

#ifdef _USE_DL_
static int have_png;
static int have_jpeg;
static int have_gif;
static int have_tiff;
static int have_xpm;
#endif

void make_little_endian(void *ptr, size_t size){
#if defined(__osf__) || defined(linux )
  return;
#else
  char bit1;
  char bit2;
  char bit3;
  char bit4;
  char bit5;
  char bit6;
  char bit7;
  char bit8;
  unsigned char* charptr = (unsigned char*)ptr;
  if(size==1)
    return;
  if(size==2){
    bit1 = *charptr;
    bit2 = *(charptr+1);
    *charptr = bit2;
    *(charptr+1) = bit1;
    return;
  }
  if(size==4){
    bit1 = *charptr;
    bit2 = *(charptr+1);
    bit3 = *(charptr+2);
    bit4 = *(charptr+3);
    *charptr = bit4;
    *(charptr+1) = bit3;
    *(charptr+2) = bit2;
    *(charptr+3) = bit1;
    return;
  }
  if(size==8){
    bit1 = *charptr;
    bit2 = *(charptr+1);
    bit3 = *(charptr+2);
    bit4 = *(charptr+3);
    bit5 = *(charptr+4);
    bit6 = *(charptr+5);
    bit7 = *(charptr+6);
    bit8 = *(charptr+7);
    *charptr = bit8;
    *(charptr+1) = bit7;
    *(charptr+2) = bit6;
    *(charptr+3) = bit5;
    *(charptr+4) = bit4;
    *(charptr+5) = bit3;
    *(charptr+6) = bit2;
    *(charptr+7) = bit1;

    return;
  }
#endif
}

size_t my_read(void *ptr, size_t size, size_t nmemb, FILE *stream){
  size_t nmemb_read;
  if((nmemb_read=fread(ptr,size,nmemb,stream))!=nmemb){
    printf("fread error\n");
    throw ImageInfoFReadExc();
  }
  return nmemb_read;
}

size_t my_write(const void *ptr, size_t size, size_t nmemb, FILE *stream){
  size_t nmemb_written;
  if((nmemb_written=fwrite(ptr,size,nmemb,stream))!=nmemb){
    printf("fwrite error\n");
    throw ImageInfoFWriteExc();
  }
  return nmemb_written;
}

char *get_suffix(const char *filename){

  char *suffix;
  suffix = new char[strlen(filename)];

  int len = strlen(filename);
  for(int i=len-1;i>=0;i--){
    if((filename[i] == '.')&&i<len-1){
      strncpy(suffix,filename+i+1,len-i-1);
      sprintf(suffix+len-i-1,"%c",'\0');
      break;
    }
  }

  return suffix;

}

void image_info::convert_colourspace(int colourspace_type_in){
  if(colourspace_type==IMAGEINFO_MONO)
    convert_greyscale();
  if(colourspace_type==IMAGEINFO_RGB)
    convert_rgb();
  if(colourspace_type==IMAGEINFO_RGBA)
    convert_rgba();
  if(colourspace_type==IMAGEINFO_MONOA)
    convert_greyscalea();
}

void image_info::set_colourspace_type(int colourspace_type_in){
  colourspace_type = colourspace_type_in;
  if(colourspace_type==IMAGEINFO_MONO)
    colourspace = IMAGEINFO_MONO_SIZE;
  if(colourspace_type==IMAGEINFO_RGB)
    colourspace = IMAGEINFO_RGB_SIZE;
  if(colourspace_type==IMAGEINFO_RGBA)
    colourspace = IMAGEINFO_RGBA_SIZE;
  if(colourspace_type==IMAGEINFO_MONOA)
    colourspace = IMAGEINFO_MONOA_SIZE;
}

int image_info::get_next_ppm_token(FILE *fp){
  char c;
  char buf[1024];
  int i;

  do{
    c=fgetc(fp);
    if(c=='#'){while((c=fgetc(fp))!='\n'){};}
  }while(isspace(c));
  ungetc(c,fp);

  i=0;
  do{
    c=fgetc(fp);
    sprintf(buf+i,"%c",c);i++;
    if(c=='#'){while((c=fgetc(fp))!='\n'){};}
  }while(isdigit(c));
  ungetc(c,fp);
  sprintf(buf+i,"%c",0);

  return atoi(buf);

}

image_info::~image_info(){

  if(pixels)
    delete [] pixels;

}

image_info::image_info(const std::vector<image_info> &patches, const std::vector<std::vector<int> > &pattern_in){

  /* Way of stitching together several images to create a larger one */
  /* Assumes patches are all the same size for the moment */
  /* Pattern descibes how they form a patchwork, eg: {{0,1},{2,3}} is
     0 1
     2 3
     and {{0,1,2,3}} is
     0 1 2 3
  */

  pixels = 0;
  std::vector<std::vector<int> > pattern = pattern_in;

  if(patches.size()==0){
    *this = image_info(patches[0]);
    return;
  }

  if(pattern.size()==0){
    /* Assume ordered square arrangement of input */
    if( (sqrt(double(patches.size()))*sqrt(double(patches.size())) - double(patches.size())) > 1e-5){
      std::cout << "No stitching pattern given and non-square number of input pictures:\nQuitting\n";
      return;
    }
    int ii = 0;
    unsigned dim = unsigned(floor(sqrt(double(patches.size()))));
    std::cout << "Assuming square of size " << dim << "\n";
    pattern.clear();
    for(unsigned int i = 0; i< dim; i++){
      pattern.push_back(std::vector<int>(dim));
      for(unsigned int j = 0; j< dim; j++){
	pattern[i][j] = ii;
	ii++;
      }
    }
  }

  width = 0;
  for(unsigned int i = 0; i< pattern[0].size(); i++){
    width += patches[pattern[0][i]].get_width();
  }

  height = 0;
  for(unsigned int j = 0; j< pattern.size(); j++){
    height += patches[pattern[j][0]].get_height();
  }

  colourspace = patches[0].get_colourspace();
  colourspace_type = patches[0].colourspace_type;

  pixels = new unsigned char[width*height*colourspace];

  for(int ii = 0; ii<width*height*colourspace;ii++)
    pixels[ii] = 0;

  for(unsigned int j = 0; j< pattern.size(); j++){
    int current_row_width = 0;
    for(unsigned int i = 0; i< pattern[j].size(); i++){
      for(int row = 0; row<patches[pattern[j][0]].get_height(); row++){
	for(int col = 0; col<patches[pattern[j][i]].get_width(); col++){

	  for(int ic=0;ic<colourspace;ic++)

	    pixels[width*(j*patches[pattern[j][0]].get_height()+row)*colourspace +
		  (current_row_width+col)*colourspace +ic ] 

	      = patches[pattern[j][i]].pixels[(row*patches[pattern[j][i]].get_width()+col)*colourspace +ic ];
	  
	}
      }
      current_row_width += patches[pattern[j][i]].get_width();
    }
  }
}

int image_info::write(const char *filename, int quality) const {
  char *suffix;
  suffix = get_suffix(filename);

  std::cout << "::write() suffix is " << suffix << std::endl;

  if(!pixels){
    printf("No pixel data in image_info object, will not write %s\n",filename);
    return 0;
  }

  if((!strcasecmp(suffix,"jpg"))||(!strcasecmp(suffix,"jpeg"))){
    try {
      writejpg(filename, quality);
    }
    catch (...) {
      printf("Error writing JPG file\n");
      return 0;
    }
  }else if(!strcasecmp(suffix,"png")){
    try {
       writepng(filename);
    }
    catch (...) {
      printf("Error writing PNG file\n");
      return 0;
    }
  }else if(!strcasecmp(suffix,"ppm")){
    try {
      writeppm(filename);
    }
    catch (...) {
      printf("Error writing PPM file\n");
      return 0;
    }
  }else if(!strcasecmp(suffix,"yuv")){
    try {
      writeyuv(filename);
    }
    catch (...) {
      printf("Error writing YUV file\n");
      return 0;
    }
  }else if(!strcasecmp(suffix,"gif")){
    try {
      writegif(filename);
    }
    catch (...) {
      printf("Error writing GIF file\n");
      return 0;
    }
  }else if(!strcasecmp(suffix,"wbmp")){
    try {
      writewbmp(filename);
    }
    catch (...) {
      printf("Error writing WBMP file\n");
      return 0;
    }
  }else if((!strcasecmp(suffix,"tif"))||(!strcasecmp(suffix,"tiff"))){
    try {
      writetif(filename);
    }
    catch (...) {
      printf("Error writing TIF file\n");
      return 0;
    }
  }else if(!strcasecmp(suffix,"xbm")){
    try {
      writexbm(filename);
    }
    catch (...) {
      printf("Error writing X11 bitmap file\n");
      return 0;
    }
  }else if(!strcasecmp(suffix,"xpm")){
    try {
      writexpm(filename);
    }
    catch (...) {
      printf("Error writing XPM file\n");
      return 0;
    }
  }else if(!strcasecmp(suffix,"bmp")){
    try {
      writebmp(filename);
    }
    catch (...) {
      printf("Error writing BMP file\n");
      return 0;
    }
  }else{
    printf("Cannot write %s\n",filename);
    printf("Can only write in png, ppm (raw), yuv, tiff, bmp, wbmp,\n"); // Need to add pbm/pgm/rgba
    printf("xpm*, xbm, gif* and jpeg at the moment.\n");
    printf("* = 256 colours max\n");
    return 0;
  }

  return 1;
}

void image_info::writeppm(const char *filename) const {
  FILE *fp;
  if ((fp = fopen(filename, "wb")) == NULL) {
    fprintf(stderr, "can't open %s in writeppm\n", filename);
    throw ImageInfoWritePPMExc();
  }

  fprintf(fp, "P6\n# CREATOR: Coot using CCP4's Write PPM util\n%d %d\n%d\n",
          width, height, UCHAR_MAX);

  if(colourspace_type!=IMAGEINFO_RGB){
    image_info tmp = image_info(*this);
    tmp.convert_rgb();
    tmp.write(filename);
    fclose(fp);
    return;
  }

  for(int i=0;i<height;i++){
    try {
      my_write(pixels + i*width*IMAGEINFO_RGB_SIZE, sizeof(unsigned char), IMAGEINFO_RGB_SIZE*width, fp);
    }
    catch (...) {
      throw ImageInfoWritePPMExc();
    }
  }
  fclose(fp);
}

void image_info::invert_colourmap(){

  image_info inverted_iinfo;

  inverted_iinfo.width = width;
  inverted_iinfo.height = height;
  inverted_iinfo.pixels = new unsigned char[width*height*colourspace];

  if(colourspace_type==IMAGEINFO_MONO||colourspace_type==IMAGEINFO_RGB||colourspace_type==IMAGEINFO_CMAP){
    for(int i=0;i<width*height*colourspace;i++)
     inverted_iinfo.pixels[i] = 255 - pixels[i];
  }else if(colourspace_type==IMAGEINFO_MONOA){
    for(int i=0;i<width*height*colourspace;i++)
      if(i%2!=1){
        inverted_iinfo.pixels[i] = 255 - pixels[i];
      }else{
        inverted_iinfo.pixels[i] = pixels[i];
      }
  }else if(colourspace_type==IMAGEINFO_RGBA){
    for(int i=0;i<width*height*colourspace;i++)
      if(i%4!=3){
        inverted_iinfo.pixels[i] = 255 - pixels[i];
      }else{
        inverted_iinfo.pixels[i] = pixels[i];
      }
  }

  memcpy(pixels,inverted_iinfo.pixels,width*height*colourspace);

}

image_info_yuv_t image_info::getyuv(bool subsample) const {

  image_info yuv_iinfo(width,height,pixels,colourspace_type);
  yuv_iinfo.colourspace_type = colourspace_type;

  yuv_iinfo.convert_yuv();
  yuv_iinfo.ScaleImage(((width+1)/2)*2, ((height+1)/2)*2);

  unsigned char *y = new unsigned char[yuv_iinfo.width*yuv_iinfo.height];
  unsigned char *u = new unsigned char[yuv_iinfo.width*yuv_iinfo.height];
  unsigned char *v = new unsigned char[yuv_iinfo.width*yuv_iinfo.height];

  int ii = 0;
  for(int i=0;i<yuv_iinfo.width*yuv_iinfo.height*yuv_iinfo.colourspace;i+=3){
    y[ii] = yuv_iinfo.pixels[i];
    u[ii] = yuv_iinfo.pixels[i+1];
    v[ii] = yuv_iinfo.pixels[i+2];
    ii++;
  }

  /* U,V subsample */
  if(subsample){
    unsigned char *usub = new unsigned char[(yuv_iinfo.width*yuv_iinfo.height)/4];
    unsigned char *vsub = new unsigned char[(yuv_iinfo.width*yuv_iinfo.height)/4];
    for(int i=0;i<yuv_iinfo.height;i+=2){
      for(int j=0;j<yuv_iinfo.width;j+=2){
        usub[i*yuv_iinfo.width/4+j/2] = (unsigned char)(float(u[i*yuv_iinfo.width+j] + u[i*yuv_iinfo.width+j+1] + u[(i+1)*yuv_iinfo.width+j] + u[(i+1)*yuv_iinfo.width+j+1])/4.0);
        vsub[i*yuv_iinfo.width/4+j/2] = (unsigned char)(float(v[i*yuv_iinfo.width+j] + v[i*yuv_iinfo.width+j+1] + v[(i+1)*yuv_iinfo.width+j] + v[(i+1)*yuv_iinfo.width+j+1])/4.0);
      }
    }
    u = usub;
    v = vsub;
  }
  image_info_yuv_t yuv;
  yuv.y = y;
  yuv.u = u;
  yuv.v = v;
  yuv.width = yuv_iinfo.width;
  yuv.height = yuv_iinfo.height;
  yuv._dum = 0; // compatibility (hopefully) with libfame.

  return yuv;
}

void image_info::writeyuv(const char *filename) const {

  image_info_yuv_t yuv = getyuv(true);
  FILE *fp;

  if ((fp = fopen(filename, "wb")) == NULL) {
    fprintf(stderr, "can't open %s in writeyuv\n", filename);
    throw ImageInfoWriteYUVExc();
  }

  for(unsigned i=0;i<yuv.height;i++){
    try {
      my_write(yuv.y + i*yuv.width, sizeof(unsigned char), yuv.width, fp);
      fflush(fp);
    }
    catch (...) {
      throw ImageInfoWriteYUVExc();
    }
  }
  for(unsigned i=0;i<yuv.height/2;i++){
    try {
      my_write(yuv.u + i*yuv.width/2, sizeof(unsigned char), yuv.width/2, fp);
      fflush(fp);
    }
    catch (...) {
      throw ImageInfoWriteYUVExc();
    }
  }
  for(unsigned i=0;i<yuv.height/2;i++){
    try {
      my_write(yuv.v + i*yuv.width/2, sizeof(unsigned char), yuv.width/2, fp);
      fflush(fp);
    }
    catch (...) {
      throw ImageInfoWriteYUVExc();
    }
  }
  delete [] yuv.y;
  delete [] yuv.u;
  delete [] yuv.v;

  fclose(fp);
}

void image_info::convert_yuv(){
  if(colourspace_type==IMAGEINFO_YUV){
    return;
  }

  convert_rgb();
  for(int i=0;i<width*height*colourspace;i+=3){
    double y = pixels[i]*0.299 + pixels[i+1]*0.587 + pixels[i+2]*0.114;
    double u = pixels[i]*-0.169 + pixels[i+1]*-0.332 + pixels[i+2]*0.5 +128;
    double v = pixels[i]*0.5 + pixels[i+1]*-0.419 + pixels[i+2]*-0.0813 +128;
    pixels[i]   = (unsigned char)(floor(y));
    pixels[i+1] = (unsigned char)(floor(u));
    pixels[i+2] = (unsigned char)(floor(v));
  }
  colourspace_type = IMAGEINFO_YUV;
  colourspace = IMAGEINFO_YUV_SIZE;
}

void image_info::convert_greyscale(){
  if(colourspace_type==IMAGEINFO_MONO)
    return;
  convert_rgb();

  image_info grey_iinfo;
  grey_iinfo.width = width;
  grey_iinfo.height = height;
  grey_iinfo.pixels = new unsigned char[width*height*IMAGEINFO_MONO_SIZE];

  int ii = 0;
  for(int i=0;i<width*height*colourspace;i+=3){
    double y = pixels[i]*0.299 + pixels[i+1]*0.587 + pixels[i+2]*0.114;
    grey_iinfo.pixels[ii]   = (unsigned char)(floor(y)+0.5);
    ii++;
  }

  colourspace = IMAGEINFO_MONO_SIZE;
  colourspace_type = IMAGEINFO_MONO;
  if(pixels)
  delete [] pixels;
  pixels = new unsigned char[width*height*IMAGEINFO_MONO_SIZE];
  memcpy(pixels,grey_iinfo.pixels,width*height*IMAGEINFO_MONO_SIZE);

}


void image_info::convert_rgba(){
  if(colourspace_type==IMAGEINFO_RGBA)
    return;

  image_info inverted_iinfo;
  inverted_iinfo.width = width;
  inverted_iinfo.height = height;

  inverted_iinfo.pixels = new unsigned char[width*height*IMAGEINFO_RGBA_SIZE];

  if(colourspace_type==IMAGEINFO_MONO){
    int ii=0;
    for(int i=0;i<width*height*IMAGEINFO_MONO_SIZE;i++){
      inverted_iinfo.pixels[ii++] = pixels[i];
      inverted_iinfo.pixels[ii++] = pixels[i];
      inverted_iinfo.pixels[ii++] = pixels[i];
      inverted_iinfo.pixels[ii++] = 255;
    }
  }else if(colourspace_type==IMAGEINFO_MONOA){
    int ii=0;
    for(int i=0;i<width*height*IMAGEINFO_MONOA_SIZE;i+=2){
      inverted_iinfo.pixels[ii++] = pixels[i];
      inverted_iinfo.pixels[ii++] = pixels[i];
      inverted_iinfo.pixels[ii++] = pixels[i];
      inverted_iinfo.pixels[ii++] = pixels[i+1];
    }
  }else if(colourspace_type==IMAGEINFO_RGB){
    int ii = 0;
    for(int i=0;i<width*height*IMAGEINFO_RGB_SIZE;i+=3){
      inverted_iinfo.pixels[ii++] = pixels[i];
      inverted_iinfo.pixels[ii++] = pixels[i+1];
      inverted_iinfo.pixels[ii++] = pixels[i+2];
      inverted_iinfo.pixels[ii++] = 255;
    }
  }else{
    convert_rgb();
    convert_rgba();
  }

  inverted_iinfo.colourspace = IMAGEINFO_RGBA_SIZE;

  colourspace = IMAGEINFO_RGBA_SIZE;
  colourspace_type = IMAGEINFO_RGBA;

  delete [] pixels;
  pixels = new unsigned char[width*height*IMAGEINFO_RGBA_SIZE];
  memcpy(pixels,inverted_iinfo.pixels,width*height*IMAGEINFO_RGBA_SIZE);

}

void image_info::convert_greyscalea(){
  if(colourspace_type==IMAGEINFO_MONOA)
    return;

  image_info inverted_iinfo;
  inverted_iinfo.width = width;
  inverted_iinfo.height = height;

  inverted_iinfo.pixels = new unsigned char[width*height*IMAGEINFO_MONOA_SIZE];

  if(colourspace_type==IMAGEINFO_MONO){
    int ii=0;
    for(int i=0;i<width*height*IMAGEINFO_MONO_SIZE;i++){
      inverted_iinfo.pixels[ii++] = pixels[i];
      inverted_iinfo.pixels[ii++] = 255;
    }
  }else if(colourspace_type==IMAGEINFO_RGB){
    int ii = 0;
    for(int i=0;i<width*height*IMAGEINFO_RGB_SIZE;i+=3){
      inverted_iinfo.pixels[ii++] = (unsigned char)(floor(double(pixels[i])*0.299 + double(pixels[i+1])*0.587 + double(pixels[i+2])*0.114)+0.5);
      inverted_iinfo.pixels[ii++] = 255;
    }
  }else if(colourspace_type==IMAGEINFO_RGBA){
    int ii = 0;
    for(int i=0;i<width*height*IMAGEINFO_RGBA_SIZE;i+=4){
      inverted_iinfo.pixels[ii++] = (unsigned char)(floor(double(pixels[i])*0.299 + double(pixels[i+1])*0.587 + double(pixels[i+2])*0.114)+0.5);
      inverted_iinfo.pixels[ii++] = pixels[i+4];
    }
  }else{
    convert_greyscale();
    convert_greyscalea();
  }

  inverted_iinfo.colourspace = IMAGEINFO_MONOA_SIZE;

  colourspace = IMAGEINFO_MONOA_SIZE;
  colourspace_type = IMAGEINFO_MONOA;

  delete [] pixels;
  pixels = new unsigned char[width*height*IMAGEINFO_MONOA_SIZE];
  memcpy(pixels,inverted_iinfo.pixels,width*height*IMAGEINFO_MONOA_SIZE);

}


void image_info::convert_rgb(){

  if(colourspace_type==IMAGEINFO_RGB){
    return;
  }

  image_info inverted_iinfo;
  inverted_iinfo.width = width;
  inverted_iinfo.height = height;

  inverted_iinfo.pixels = new unsigned char[width*height*IMAGEINFO_RGB_SIZE];

  if(colourspace_type==IMAGEINFO_MONO){
    int ii=0;
    for(int i=0;i<width*height*IMAGEINFO_MONO_SIZE;i++){
      inverted_iinfo.pixels[ii++] = pixels[i];
      inverted_iinfo.pixels[ii++] = pixels[i];
      inverted_iinfo.pixels[ii++] = pixels[i];
    }
  }

  if(colourspace_type==IMAGEINFO_MONOA){
    int ii=0;
    for(int i=0;i<width*height*IMAGEINFO_MONOA_SIZE;i+=2){
      inverted_iinfo.pixels[ii++] = pixels[i];
      inverted_iinfo.pixels[ii++] = pixels[i];
      inverted_iinfo.pixels[ii++] = pixels[i];
    }
  }

  if(colourspace_type==IMAGEINFO_RGBA){
    int ii = 0;
    for(int i=0;i<width*height*colourspace;i++){
      if(i%4!=3){
        inverted_iinfo.pixels[ii] = pixels[i];
        ii++;
      }
    }
  }

  inverted_iinfo.colourspace = IMAGEINFO_RGB_SIZE;
  colourspace = IMAGEINFO_RGB_SIZE;

  colourspace_type = IMAGEINFO_RGB;
  delete [] pixels;
  pixels = new unsigned char[width*height*IMAGEINFO_RGB_SIZE];
  memcpy(pixels,inverted_iinfo.pixels,width*height*IMAGEINFO_RGB_SIZE);

}

void image_info::invert(){

  image_info inverted_iinfo;

  inverted_iinfo.width = width;
  inverted_iinfo.height = height;
  inverted_iinfo.pixels = new unsigned char[width*height*IMAGEINFO_RGBA_SIZE];

  int iprime = 0;
  for(int i=height-1;i>=0;i--,iprime++)
    memcpy(inverted_iinfo.pixels+ i*width*colourspace, pixels + iprime*width*colourspace, sizeof(unsigned char)*colourspace*width);

  memcpy(pixels,inverted_iinfo.pixels,width*height*colourspace);

}

image_info::image_info(){
  pixels = 0;
  pixels = 0;
  colourspace_type = -1;
}

image_info const& image_info::operator=(image_info const &iinfo){

  if (this != &iinfo) {
  pixels = 0;
    copy(iinfo.get_width(),iinfo.get_height(),iinfo.get_pixels(),iinfo.get_colourspace_type());
  }
  return *this;
}

image_info::image_info(const image_info &iinfo){
  pixels = 0;
  copy(iinfo.get_width(),iinfo.get_height(),iinfo.get_pixels(),iinfo.get_colourspace_type());

}

image_info::image_info(int width_in, int height_in, unsigned char *pixels_in, int colourspace_type_in){
  pixels = 0;
  copy(width_in, height_in, pixels_in, colourspace_type_in);
}

void image_info::copy(int width_in, int height_in, unsigned char *pixels_in, int colourspace_type_in){

  width = width_in;
  height = height_in;
  colourspace_type = colourspace_type_in; 
  if(colourspace_type==IMAGEINFO_MONO)
    colourspace = IMAGEINFO_MONO_SIZE;
  if(colourspace_type==IMAGEINFO_RGB)
    colourspace = IMAGEINFO_RGB_SIZE;
  if(colourspace_type==IMAGEINFO_RGBA)
    colourspace = IMAGEINFO_RGBA_SIZE;
  if(colourspace_type==IMAGEINFO_MONOA)
    colourspace = IMAGEINFO_MONOA_SIZE;

  if(pixels)
    delete [] pixels;
  pixels = new unsigned char[width*height*colourspace];
  memcpy(pixels,pixels_in,width*height*colourspace);
  
}

void image_info::set_bitmap_data(int width_in, int height_in, unsigned char *bm){
    int bm_width = (width_in + 7) / 8;
    colourspace = IMAGEINFO_RGB_SIZE;
    colourspace_type = IMAGEINFO_RGB;
    width = 8*bm_width;
    height = height_in;
    pixels = new unsigned char[8*bm_width*height*IMAGEINFO_RGB_SIZE];
    unsigned int val, p;
    unsigned char *pixbuf = new unsigned char[1];
    for(int i=0;i<height;i++){
      for(int j=0;j<bm_width;j++){
	int index = (i*bm_width+j)*IMAGEINFO_RGB_SIZE*8;
	pixbuf[0] = bm[i*bm_width+j];
	val = pixbuf[0];
	p = ((val & 128)/128)*255;
	pixels[index]   = p;
	pixels[index+1] = p;
	pixels[index+2] = p;
	p = ((val & 64)/64)*255;
	pixels[index+3] = p;
	pixels[index+4] = p;
	pixels[index+5] = p;
	p = ((val & 32)/32)*255;
	pixels[index+6] = p;
	pixels[index+7] = p;
	pixels[index+8] = p;
	p = ((val & 16)/16)*255;
	pixels[index+9]  = p;
	pixels[index+10] = p;
	pixels[index+11] = p;
	p = ((val & 8)/8)*255;
	pixels[index+12] = p;
	pixels[index+13] = p;
	pixels[index+14] = p;
	p = ((val & 4)/4)*255;
	pixels[index+15] = p;
	pixels[index+16] = p;
	pixels[index+17] = p;
	p = ((val & 2)/2)*255;
	pixels[index+18] = p;
	pixels[index+19] = p;
	pixels[index+20] = p;
	p = ((val & 1))*255;
	pixels[index+21] = p;
	pixels[index+22] = p;
	pixels[index+23] = p;
      }
    }
    delete [] pixbuf;
}
 
image_info::image_info(const char *filename){
  pixels = 0;
  width = 0;
  height = 0;
  pixels = 0;
  colourspace = 0;
  colourspace_type = -1;

  read(filename); // And do what?
}

int image_info::read(const char *filename){

  FILE *fp;
  char *suffix;

  if ((fp = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "can't open %s in image_info::read(char *filename)\n", filename);
    return 0;
  }
  fclose(fp);

  suffix = get_suffix(filename);

  if(!strcasecmp(suffix,"png")){
    try { 
      readpng(filename);
    }
    catch (...) {
      printf("Error reading PNG file\n");
      return 0;
    }
  }
  else if((!strcasecmp(suffix,"jpg"))||(!strcasecmp(suffix,"jpeg"))){
    try {
      readjpg(filename);
    }
    catch (...) {
      printf("Error reading JPG file\n");
      return 0;
    }
  }
  else if(!strcasecmp(suffix,"rgba")){
    try {
      readrgba(filename);
    }
    catch (...) {
      printf("Error reading RGB file\n");
      return 0;
    }
    colourspace = IMAGEINFO_RGBA_SIZE;
    colourspace_type = IMAGEINFO_RGBA;
  }
  else if(!strcasecmp(suffix,"bmp")){
    try {
      readbmp(filename);
    }
    catch (...) {
      printf("Error reading BMP file\n");
      return 0;
    }
  }
  else if(!strcasecmp(suffix,"gif")){
    try {
      readgif(filename);
    }
    catch (...) {
      printf("Error reading GIF file\n");
      return 0;
    }
  }
  else if((!strcasecmp(suffix,"tif"))||!strcasecmp(suffix,"tiff")){
    try {
      readtif(filename);
    }
    catch (...) {
      printf("Error reading TIF file\n");
      return 0;
    }
  }
  else if((!strcasecmp(suffix,"ppm"))||(!strcasecmp(suffix,"pgm"))||(!strcasecmp(suffix,"pbm"))){
    try {
      readppm(filename);
    }
    catch (...) {
      printf("Error reading PPM file\n");
      return 0;
    }
  }else if(!strcasecmp(suffix,"xbm")){
    try {
      readxbm(filename);
    }
    catch (...) {
      printf("Error reading X11 bitmap file\n");
      return 0;
    }
  }else if(!strcasecmp(suffix,"xpm")){
    try {
      readxpm(filename);
    }
    catch (...) {
      printf("Error reading XPM file\n");
      return 0;
    }
  }else{
    printf("Can only read in pbm/pgm/ppm (raw or ascii),\n"); // Need to add wbmp
    printf("Windows BMP, raw RGBA, png, gif, tiff, xpm, xbm, \n");
    printf("and jpeg at the moment.\n");
    return 0;
  }
  return 1;
}

void image_info::writebmp(const char *filename) const {
  FILE * fp;
  if ((fp = fopen(filename, "wb")) == NULL) {
    fprintf(stderr, "can't open %s in writebmp\n", filename);
    throw ImageInfoWriteBMPExc();
    return;
  }
  int16_t type = 19778;
  try {
    my_write(&type, sizeof(int16_t), 1, fp);
  }
  catch (...) {
    throw ImageInfoWriteBMPExc();
  }
  int32_t size = 54+width*height*get_colourspace();
  try {
    my_write(&size, sizeof(int32_t), 1, fp);
  }
  catch (...) {
    throw ImageInfoWriteBMPExc();
  }
  int16_t res1=0,res2=0;
  try {
    my_write(&res1, sizeof(int16_t), 1, fp);
  }
  catch (...) {
    throw ImageInfoWriteBMPExc();
  }
  try {
    my_write(&res2, sizeof(int16_t), 1, fp);
  }
  catch (...) {
    throw ImageInfoWriteBMPExc();
  }
  int32_t offset;
  if(colourspace_type==IMAGEINFO_MONO)
    offset  = 1078;
  else
    offset  = 54;
  try {
    my_write(&offset, sizeof(int32_t), 1, fp);
  }
  catch (...) {
    throw ImageInfoWriteBMPExc();
  }
  /* End of header */
  int32_t header_size = 40;
  try {
    my_write(&header_size, sizeof(int32_t), 1, fp);
  }
  catch (...) {
    throw ImageInfoWriteBMPExc();
  }
  int32_t w = width;
  try {
    my_write(&w, sizeof(int32_t), 1, fp);
  }
  catch (...) {
    throw ImageInfoWriteBMPExc();
  }
  int32_t h = height;
  try {
    my_write(&h, sizeof(int32_t), 1, fp);
  }
  catch (...) {
    throw ImageInfoWriteBMPExc();
  }
  int16_t nplanes = 1;
  try {
    my_write(&nplanes, sizeof(int16_t), 1, fp);
  }
  catch (...) {
    throw ImageInfoWriteBMPExc();
  }
  int16_t bpp;
  if(colourspace_type == IMAGEINFO_RGB)
    bpp = 24;
  if(colourspace_type == IMAGEINFO_RGBA)
    bpp = 32;
  if(colourspace_type == IMAGEINFO_MONO)
    bpp = 8;
  try {
    my_write(&bpp, sizeof(int16_t), 1, fp);
  }
  catch (...) {
    throw ImageInfoWriteBMPExc();
  }
  int32_t comp = 0;
  try {
    my_write(&comp, sizeof(int32_t), 1, fp);
  }
  catch (...) {
    throw ImageInfoWriteBMPExc();
  }
  int32_t isize = width*height*colourspace;
  try {
    my_write(&isize, sizeof(int32_t), 1, fp);
  }
  catch (...) {
    throw ImageInfoWriteBMPExc();
  }
  int32_t xres = 600, yres = 600;
  try {
    my_write(&xres, sizeof(int32_t), 1, fp);
  }
  catch (...) {
    throw ImageInfoWriteBMPExc();
  }
  try {
    my_write(&yres, sizeof(int32_t), 1, fp);
  }
  catch (...) {
    throw ImageInfoWriteBMPExc();
  }
  int32_t ncol = 0, nimpcol = 0;
  try {
    my_write(&ncol, sizeof(int32_t), 1, fp);
  }
  catch (...) {
    throw ImageInfoWriteBMPExc();
  }
  try {
    my_write(&nimpcol, sizeof(int32_t), 1, fp);
  }
  catch (...) {
    throw ImageInfoWriteBMPExc();
  }
  /* End of info header */

  unsigned char *tmppix = 0;
  if(colourspace_type==IMAGEINFO_RGB){
    tmppix = new unsigned char[width*height*IMAGEINFO_RGB_SIZE];
    for (int i=0; i<width*height*IMAGEINFO_RGB_SIZE;i+=IMAGEINFO_RGB_SIZE){
      tmppix[i]   = pixels[i+2];
      tmppix[i+1] = pixels[i+1];
      tmppix[i+2] = pixels[i];
    }
  }
  if(colourspace_type==IMAGEINFO_RGBA){
    tmppix = new unsigned char[width*height*IMAGEINFO_RGBA_SIZE];
    for (int i=0; i<width*height*IMAGEINFO_RGBA_SIZE;i+=IMAGEINFO_RGBA_SIZE){
      tmppix[i]   = pixels[i+2];
      tmppix[i+1] = pixels[i+1];
      tmppix[i+2] = pixels[i];
      tmppix[i+3] = pixels[i+3];
    }
  }
  
  if(colourspace_type==IMAGEINFO_MONO){
    fseek(fp, 54, SEEK_SET);
    unsigned char junk[1024];
    for(int i=0; i<256;i++){
      junk[i*4] = i;
      junk[i*4+1] = i;
      junk[i*4+2] = i;
      junk[i*4+3] = 1;
    }
    my_write(junk, sizeof (unsigned char), 1024, fp);
  }

  int PhysicalWidth = ((width * (bpp / 8)) + 3) & (~3);
  int ScanPadding = PhysicalWidth - (width * ( bpp / 8 ));
  unsigned char *tmpbuf = new unsigned char[ScanPadding];
  for(int i=height-1;i>=0;i--){
    try {
      if(tmppix)
        my_write(tmppix+(i*width)*colourspace, sizeof (unsigned char), width*colourspace, fp);
      else
        my_write(pixels+(i*width)*colourspace, sizeof (unsigned char), width*colourspace, fp);
    }
    catch (...) {
      throw ImageInfoWriteBMPExc();
    }
    try {
      my_write(tmpbuf, sizeof (unsigned char), ScanPadding, fp);
    }
    catch (...) {
      throw ImageInfoWriteBMPExc();
    }
  }
  if(tmppix)
    delete [] tmppix;
   
  fclose(fp);
  return;
}

void image_info::readbmp(const char *filename)
{
    FILE *fp;
    int16_t type;
    int32_t offset;
    int16_t nplanes;
    int16_t bpp;
    unsigned char temp;

    fp = fopen(filename, "rb");

    try {
      my_read(&type, sizeof(int16_t), 1, fp);
      make_little_endian(&type,sizeof(int16_t));
    }
    catch (...) {
      throw ImageInfoReadBMPExc();
    }
    if(type!=19778){
      fprintf(stderr,"Error: not a BMP file\n");
        throw ImageInfoReadBMPExc();
      return;
    }
    fseek(fp, 8, SEEK_CUR);
    try {
      my_read(&offset, sizeof(int32_t), 1, fp);
      make_little_endian(&offset,sizeof(int32_t));
    }
    catch (...) {
      fprintf(stderr,"Error getting offset\n");
      throw ImageInfoReadBMPExc();
    }
    fseek(fp, 4, SEEK_CUR);
    try {
      my_read(&width, sizeof(int), 1, fp);
      make_little_endian(&width, sizeof(int));
    }
    catch (...) {
      fprintf(stderr,"Error getting width\n");
      throw ImageInfoReadBMPExc();
    }
    try {
      my_read(&height, sizeof(int), 1, fp);
      make_little_endian(&height, sizeof(int));
    }
    catch (...) {
      fprintf(stderr,"Error getting height\n");
      throw ImageInfoReadBMPExc();
    }
    try {
      my_read(&nplanes, sizeof(int16_t), 1, fp);
      make_little_endian(&nplanes, sizeof(int16_t));
    }
    catch (...) {
      fprintf(stderr,"Error getting nplanes\n");
      throw ImageInfoReadBMPExc();
    }
    if (nplanes != 1)    {
        printf("Error: number of Planes not 1!, %d\n",nplanes);
        throw ImageInfoReadBMPExc();
        return;
    }
    try {
      my_read(&bpp, sizeof(int16_t), 1, fp);
      make_little_endian(&bpp, sizeof(int16_t));
    }
    catch (...) {
      fprintf(stderr,"Error getting nplanes\n");
    }

    int32_t comp;
    try {
      my_read(&comp, sizeof(int32_t), 1, fp);
      make_little_endian(&comp, sizeof(int32_t));
    }
    catch (...) {
      fprintf(stderr,"Error getting comp\n");
    }
    int32_t bitmap_size;
    try {
      my_read(&bitmap_size, sizeof(int32_t), 1, fp);
      make_little_endian(&bitmap_size, sizeof(int32_t));
    }
    catch (...) {
      fprintf(stderr,"Error getting bitmap_size\n");
    }
    int32_t horiz_res;
    try {
      my_read(&horiz_res, sizeof(int32_t), 1, fp);
      make_little_endian(&horiz_res, sizeof(int32_t));
    }
    catch (...) {
      fprintf(stderr,"Error getting horiz_res\n");
    }
    int32_t vert_res;
    try {
      my_read(&vert_res, sizeof(int32_t), 1, fp);
      make_little_endian(&vert_res, sizeof(int32_t));
    }
    catch (...) {
      fprintf(stderr,"Error getting vert_res\n");
    }
    int32_t colours_used;
    try {
      my_read(&colours_used, sizeof(int32_t), 1, fp);
      make_little_endian(&colours_used, sizeof(int32_t));
    }
    catch (...) {
      fprintf(stderr,"Error getting colours_used\n");
    }

    int ncolours = 0;
    if(offset!=54){
      ncolours = (offset - 54)/4;
    }

    int32_t ncolheader;
    fseek(fp, offset-2*sizeof(int32_t), SEEK_SET);
    try {
      my_read(&ncolheader, sizeof(int32_t), 1, fp);
    }
    catch (...) {
      fprintf(stderr,"Error getting ncolheader\n");
    }

    if(bpp<0){
      /* broken file, assume mono */
      bpp = 8;
      ncolours = 0;
    }

    int PhysicalWidth = ((width * (bpp / 8)) + 3) & (~3);
    int ScanPadding = PhysicalWidth - (width * ( bpp / 8 ));
    unsigned char *tmpbuf = new unsigned char[ScanPadding];

    pixels  = new unsigned char[width*height*bpp/8];
    fseek(fp, offset, SEEK_SET);

    for(int i=height-1;i>=0;i--){
      try {
	my_read(pixels+(i*width)*bpp/8, sizeof (unsigned char), width*bpp/8, fp);
      }
      catch (...) {
        fprintf(stderr,"Error reading pixels, %d of %d\n",i,height);
	throw ImageInfoReadBMPExc();
      }
      try {
	my_read(tmpbuf, sizeof (unsigned char), ScanPadding, fp);
      }
      catch (...) {
        fprintf(stderr,"Error reading scanpadding\n");
	throw ImageInfoReadBMPExc();
      }
    }
    delete [] tmpbuf;

    if (bpp == 24){
      colourspace = IMAGEINFO_RGB_SIZE;
      colourspace_type = IMAGEINFO_RGB;
    }else if(bpp == 32){
      colourspace = IMAGEINFO_RGBA_SIZE;
      colourspace_type = IMAGEINFO_RGBA;
    }else if(bpp == 8&&ncolours==0){
      colourspace = IMAGEINFO_MONO_SIZE;
      colourspace_type = IMAGEINFO_MONO;
    }else if(bpp == 8&&ncolours>0&&ncolours<257){
      colourspace = IMAGEINFO_RGBA_SIZE;
      colourspace_type = IMAGEINFO_RGBA;
    }else{
        printf("Bits per Pixel not supported: (%d)\n",bpp);
        throw ImageInfoReadBMPExc();
        return;
    }

    if(ncolours>0&&bpp==8&&ncolours>0&&ncolours<257){
      fseek(fp, 54, SEEK_SET);
      unsigned char *cmap = new unsigned char[ncolours*4];
      for(int i=0;i<ncolours;i++){
        unsigned char r,g,b,a; 
        my_read(&r, sizeof(unsigned char), 1, fp);
        my_read(&g, sizeof(unsigned char), 1, fp);
        my_read(&b, sizeof(unsigned char), 1, fp);
        my_read(&a, sizeof(unsigned char), 1, fp);
        cmap[4*i] = b;
        cmap[4*i+1] = g;
        cmap[4*i+2] = r;
        cmap[4*i+3] = a;
      }
      
      unsigned char *cmapped_pixels  = new unsigned char[width*height];
      memcpy(cmapped_pixels,pixels,width*height);
      delete [] pixels;
      pixels  = new unsigned char[width*height*colourspace];
      for (int i=0; i<width*height;i++) {
        pixels[i*colourspace]   = cmap[colourspace*cmapped_pixels[i]];
        pixels[i*colourspace+1] = cmap[colourspace*cmapped_pixels[i]+1];
        pixels[i*colourspace+2] = cmap[colourspace*cmapped_pixels[i]+2];
        pixels[i*colourspace+3] = 255; // This is debatable ....

      }
      delete [] cmap;
    }

    fclose(fp);

    if(colourspace_type==IMAGEINFO_RGB){
      for (int i=0; i<width*height*IMAGEINFO_RGB_SIZE;i+=IMAGEINFO_RGB_SIZE)
	{
	  temp = pixels[i];
	  pixels[i] = pixels[i + 2];
	  pixels[i + 2] = temp;
	}
    }
    if(colourspace_type==IMAGEINFO_RGBA&&ncolours==0){
      for (int i=0; i<width*height*IMAGEINFO_RGBA_SIZE;i+=IMAGEINFO_RGBA_SIZE)
	{
	  temp = pixels[i];
	  pixels[i] = pixels[i + 2];
	  pixels[i + 2] = temp;
	}
    }
}

void image_info::readrgba(const char *filename){

  FILE *fp;
  unsigned long filesize;

  fp = fopen(filename, "rb");
  fseek(fp, 0, SEEK_END);
  filesize=ftell(fp);
  fseek(fp, 0, SEEK_SET);

  pixels = new unsigned char[filesize];

  try {
    my_read(pixels, sizeof(unsigned char), filesize, fp);
  }
  catch (...) {
    throw ImageInfoReadRGBExc();
  }

  fclose(fp);

  // Assume texture is square! how can we find out relative
  // dimensions? Are are they meaningless in rgba.

  width = (int)sqrt(double(filesize/4));
  height = width;

}

void image_info::readppm(const char *filename){
  FILE *fp;
  int c;
  int i;
  char buf[1024];
  int type;
  int maxval=1;

  fp = fopen(filename, "rb");

  i=0;
 
  while(i<2){
    c=fgetc(fp);
    if (c==EOF) {
      printf("Error reading ppm\n");
      throw ImageInfoReadPPMExc();
    }
    sprintf(buf+i,"%c",c);
    i++;
  }
  sprintf(buf+i,"%c",0);
  type = atoi(buf+1);

  width = get_next_ppm_token(fp);
  height = get_next_ppm_token(fp);
  if(type==2||type==5||type==3||type==6) maxval = get_next_ppm_token(fp);

  do{
    c=fgetc(fp);
    if (c==EOF) {
      printf("Error reading ppm\n");
      throw ImageInfoReadPPMExc();
    }
    if(c=='#'){
      while((c=fgetc(fp))!='\n'){
	if (c==EOF) {
	  printf("Error reading ppm\n");
	  throw ImageInfoReadPPMExc();
	}
      }
    }
  }while(isspace(c));
  fseek(fp,-1,SEEK_CUR);

  if(type==6){ // Colour Binary 
    colourspace = IMAGEINFO_RGB_SIZE;
    colourspace_type = IMAGEINFO_RGB;
    pixels = new unsigned char[width*height*IMAGEINFO_RGB_SIZE];
    try {
      my_read(pixels, sizeof(unsigned char), IMAGEINFO_RGB_SIZE*width*height, fp);
    }
    catch (...) {
      throw ImageInfoReadPPMExc();
    }
  }
  if(type==3){ // Colour Ascii 
    colourspace = IMAGEINFO_RGB_SIZE;
    colourspace_type = IMAGEINFO_RGB;
    pixels = new unsigned char[width*height*IMAGEINFO_RGB_SIZE];
    unsigned int r,g,b;
    for(int i=0;i<height;i++){
      for(int j=0;j<width;j++){
	c = fscanf(fp,"%d",&r); 
	if(c==0||c==EOF)
	  throw ImageInfoReadPPMExc();
	c = fscanf(fp,"%d",&g); 
	if(c==0||c==EOF)
	  throw ImageInfoReadPPMExc();
	c = fscanf(fp,"%d",&b);
	if(c==0||c==EOF)
	  throw ImageInfoReadPPMExc();
	pixels[i*width*IMAGEINFO_RGB_SIZE+j]   = r*255/maxval; 
	pixels[i*width*IMAGEINFO_RGB_SIZE+j+1] = g*255/maxval; 
	pixels[i*width*IMAGEINFO_RGB_SIZE+j+2] = b*255/maxval; 
      }
    }
  }
  if(type==5){ // Greyscale Binary
    pixels = new unsigned char[width*height*IMAGEINFO_MONO_SIZE];
    colourspace = IMAGEINFO_MONO_SIZE;
    colourspace_type = IMAGEINFO_MONO;
    for(int i=0;i<height;i++){
      for(int j=0;j<width;j++){
	try {
	  my_read(pixels + i*width*IMAGEINFO_MONO_SIZE+IMAGEINFO_MONO_SIZE*j, sizeof(unsigned char), 1, fp);
	}
	catch (...) {
	  throw ImageInfoReadPPMExc();
	}
	pixels[i*width*IMAGEINFO_MONO_SIZE+IMAGEINFO_MONO_SIZE*j] = pixels[i*width*IMAGEINFO_MONO_SIZE+IMAGEINFO_MONO_SIZE*j]*255/maxval;
      }
    }
  }
  if(type==2){ // Greyscale Ascii
    pixels = new unsigned char[width*height*IMAGEINFO_MONO_SIZE];
    colourspace = IMAGEINFO_MONO_SIZE;
    colourspace_type = IMAGEINFO_MONO;
    unsigned int val;
    for(int i=0;i<height;i++){
      for(int j=0;j<width;j++){
	c = fscanf(fp,"%d",&val);
	if(c==0||c==EOF)
	  throw ImageInfoReadPPMExc();
	pixels[i*width*IMAGEINFO_MONO_SIZE+IMAGEINFO_MONO_SIZE*j] = val*255/maxval; 
      }
    }
  }
  if(type==4){ // B/W Bitmap Binary
    pixels = new unsigned char[width*height*IMAGEINFO_MONO_SIZE];
    colourspace = IMAGEINFO_MONO_SIZE;
    colourspace_type = IMAGEINFO_MONO;
    unsigned int val, p;
    unsigned char pixbuf;
    for(int i=0;i<height;i++){
      for(int j=0;j<(width+7)/8;j++){
	try {
	  my_read(&pixbuf, sizeof(unsigned char), 1, fp);
	}
	catch (...) {
	  throw ImageInfoReadPPMExc();
	}
	val = pixbuf;
	p = ~((val & 128)*255/128);
        if(j*8<width)
	pixels[i*width+8*j]   = p;
	p = ~((val & 64)*255/64);
        if(j*8<width+1)
	pixels[i*width+8*j+1] = p;
	p = ~((val & 32)*255/32);
        if(j*8<width+2)
	pixels[i*width+8*j+2] = p;
	p = ~((val & 16)*255/16);
        if(j*8<width+3)
	pixels[i*width+8*j+3]  = p;
	p = ~((val & 8)*255/8);
        if(j*8<width+4)
	pixels[i*width+8*j+4] = p;
	p = ~((val & 4)*255/4);
        if(j*8<width+5)
	pixels[i*width+8*j+5] = p;
	p = ~((val & 2)*255/2);
        if(j*8<width+6)
	pixels[i*width+8*j+6] = p;
	p = ~((val & 1)*255);
        if(j*8<width+7)
	pixels[i*width+8*j+7] = p;
      }
    }
  }
  if(type==1){ // B/W Bitmap Ascii
    pixels = new unsigned char[width*height*IMAGEINFO_MONO_SIZE];
    colourspace = IMAGEINFO_MONO_SIZE;
    colourspace_type = IMAGEINFO_MONO;
    unsigned int val;
    for(int i=0;i<height;i++){
      for(int j=0;j<width;j++){
	c = fscanf(fp,"%d",&val);
	if(c==0||c==EOF)
	  throw ImageInfoReadPPMExc();
	pixels[i*width*IMAGEINFO_MONO_SIZE+IMAGEINFO_MONO_SIZE*j]   = val*255; 
      }
    }
  }

  fclose(fp);

}

void image_info::Rotate(){
  image_info inverted_iinfo;

  inverted_iinfo.width = height;
  inverted_iinfo.height = width;
  inverted_iinfo.pixels = new unsigned char[width*height*colourspace];

  int ii=0;
  for(int j=0;j<height;j++){
    for(int i=0;i<width;i++){
      int idx1 = (j * width + i ) * colourspace;
      int idx2 = (i * height + height - j -1 ) * colourspace;
      for(int k=0;k<colourspace;k++,ii++){
        inverted_iinfo.pixels[idx2+k] = pixels[idx1+k];
      }
    }
  }

  width  = inverted_iinfo.width;
  height = inverted_iinfo.height;
  delete [] pixels;
  pixels = new unsigned char[width*height*colourspace];
  memcpy(pixels,inverted_iinfo.pixels,width*height*colourspace);

}
void image_info::ScaleImage(int x, int y){

  if(x==get_width()&&y==get_height()) return;

  double current_ratio = (double)width / (double)height;
  if(x<0){
    x = int((double)y * current_ratio);
  }
  if(y<0){
    y = int((double)x / current_ratio);
  }
  image_info inverted_iinfo;

  inverted_iinfo.width = width;
  inverted_iinfo.height = height;
  inverted_iinfo.pixels = new unsigned char[x*y*colourspace];

  for(int j=0;j<y;j++){
    int hidx = j * height / y;
    for(int i=0;i<x;i++){
      int widx = i * width / x;
      for(int k=0;k<colourspace;k++){
        int idx = hidx*width*colourspace+widx*colourspace + k;
        inverted_iinfo.pixels[j*x*colourspace+i*colourspace+k] = pixels[idx];
      }
    }
  }

  width  = x;
  height = y;
  delete [] pixels;
  pixels = new unsigned char[width*height*colourspace];
  memcpy(pixels,inverted_iinfo.pixels,width*height*colourspace);

}

void image_info::readgif(const char *filename){

#ifdef USE_LIBGIF

#ifdef _USE_DL_
  if(!have_gif){
    if(init_gif(LIBGIF_SOLIBRARY)){
      have_gif = 1;
    }
  }
  if(!have_gif){
    printf("Gif library not found\n");
    printf("Please try reading from a supported format\n");
    return;
  }
#endif  

  GifFileType* giffile =  DGifOpenFileName(filename);

  width  = giffile->SWidth;
  height = giffile->SHeight;
  int err=GIF_OK;


  GifRecordType RecordType;
  ColorMapObject *ColorMap = giffile->SColorMap;
  GifByteType *Extension;
  int ExtCode;
  int got_gce = 0;
  int trans_index=0;
  unsigned char *gifpix;
  static int
            InterlacedOffset[] = { 0, 4, 2, 1 }, /* The way Interlaced image should. */
            InterlacedJumps[] = { 8, 8, 4, 2 };    /* be read - offsets and jumps... */

  do {
    if (DGifGetRecordType(giffile, &RecordType) == GIF_ERROR) {
      cout << "DGifGetRecordType error\n";
      throw ImageInfoReadGIFExc();
    }
    switch (RecordType) {
      case IMAGE_DESC_RECORD_TYPE:
        if (DGifGetImageDesc(giffile) == GIF_ERROR) {
          cout << "DGifGetImageDesc error\n";
	  throw ImageInfoReadGIFExc();
        }
	if (got_gce && trans_index >= 0){
           colourspace = IMAGEINFO_RGBA_SIZE;
	   colourspace_type = IMAGEINFO_RGBA;
	} else {
           colourspace = IMAGEINFO_RGB_SIZE;
	   colourspace_type = IMAGEINFO_RGB;
	}
        pixels = new unsigned char[width*height*colourspace];
        gifpix = new unsigned char[width];

        if (giffile->Image.Interlace){
	   for (int i = 0; i < 4; i++) {
             for (int irow = InterlacedOffset[i]; irow < height; irow += InterlacedJumps[i]) {
               err = DGifGetLine(giffile, gifpix, width);
               if(err!=GIF_OK){
		 cout << irow  << " " << i << " DGifGetLine error\n";
		 throw ImageInfoReadGIFExc();
	       }
               for(int j=0;j<width;j++){
	           pixels[(irow*width*colourspace)+(j*colourspace)]   = ColorMap->Colors[gifpix[j]].Red;
	           pixels[(irow*width*colourspace)+(j*colourspace)+1] = ColorMap->Colors[gifpix[j]].Green;
	           pixels[(irow*width*colourspace)+(j*colourspace)+2] = ColorMap->Colors[gifpix[j]].Blue;
	           if(colourspace_type==IMAGEINFO_RGBA)
		     if(trans_index == gifpix[j])
	               pixels[(irow*width*colourspace)+(j*colourspace)+3] = 0;
		     else
	               pixels[(irow*width*colourspace)+(j*colourspace)+3] = 255;
	       }
	     }
	   }
	}else{

          for(int i=0;i<height;i++){
             err = DGifGetLine(giffile, gifpix, width);
             if(err!=GIF_OK){
	       cout << i << " DGifGetLine error\n";
	       throw ImageInfoReadGIFExc();
	     }
             for(int j=0;j<width;j++){
	         pixels[(i*width*colourspace)+(j*colourspace)]   = ColorMap->Colors[gifpix[j]].Red;
	         pixels[(i*width*colourspace)+(j*colourspace)+1] = ColorMap->Colors[gifpix[j]].Green;
	         pixels[(i*width*colourspace)+(j*colourspace)+2] = ColorMap->Colors[gifpix[j]].Blue;
	         if(colourspace_type==IMAGEINFO_RGBA)
		   if(trans_index == gifpix[j])
	             pixels[(i*width*colourspace)+(j*colourspace)+3] = 0;
		   else
	             pixels[(i*width*colourspace)+(j*colourspace)+3] = 255;
	     }
          }
	}

	got_gce = 0;
	delete [] gifpix;
	break;
      case EXTENSION_RECORD_TYPE:
	if (DGifGetExtension(giffile, &ExtCode, &Extension) == GIF_ERROR) {
          cout << "DGifGetExtension error\n";
	  throw ImageInfoReadGIFExc();
        }
	if (ExtCode == 0xF9) {
        got_gce = 1;
        if (Extension[1] & 1)
          trans_index = Extension[4];
        else
          trans_index = -1;
	}
	while (Extension != NULL) {
	  if (DGifGetExtensionNext(giffile, &Extension) == GIF_ERROR) {
            cout << "DGifGetExtension error\n";
	    throw ImageInfoReadGIFExc();
          }
	}
	break;
      case TERMINATE_RECORD_TYPE:     
	break;
      default:
	break;
	{}
    }
  } while (RecordType != TERMINATE_RECORD_TYPE);

  DGifCloseFile(giffile);

#else
    printf("No gif support compiled into this program\n");
    printf("Please install lib(un)gif (if not already done) and\n");
    printf("rebuild adding -DUSE_LIBGIF to C_DEFINES\n");
    printf("and CXX_DEFINES. Or convert input file to a\n");
    printf("supported format.\n");
#endif /* USE_LIBGIF */

}

int MaxColourMapOverlap(unsigned char R, unsigned char G, unsigned char B, const std::vector<std::vector<unsigned char> > &colmap);

void image_info::writegif(const char *filename) const {

  image_info tmp = *this;
  tmp.convert_rgb();

#ifdef USE_LIBGIF

#ifdef _USE_DL_
  if(!have_gif){
    if(init_gif(LIBGIF_SOLIBRARY)){
      have_gif = 1;
    }
  }
  if(!have_gif){
    printf("Gif library not found\n");
    printf("Please try reading from a supported format\n");
    return;
  }
#endif  

  int err=GIF_OK;
  GifFileType* giffile =  EGifOpenFileName((char*)filename,FALSE);
  EGifSetGifVersion("87a");
  // All a bit tricky since we have to define a colourmap !! 
  std::vector<std::vector<unsigned char> > colmap = GetColourMap(256);
  ColorMapObject *ColorMap = MakeMapObject(colmap.size(),0);
  for(unsigned i=0;i<colmap.size();i++){
    ColorMap->Colors[i].Red = colmap[i][0];
    ColorMap->Colors[i].Green = colmap[i][1];
    ColorMap->Colors[i].Blue = colmap[i][2];
    //printf("%d %d %d\n",ColorMap->Colors[i].Red, ColorMap->Colors[i].Green,ColorMap->Colors[i].Blue); fflush(stdout);
  }

  err = EGifPutScreenDesc(giffile,width,height,ColorMap->BitsPerPixel,0,ColorMap);
  if(err!=GIF_OK)
    printf("Error: EGifPutScreenDesc\n");

  err = EGifPutImageDesc(giffile,0,0,width,height,0,0);
  if(err!=GIF_OK)
    printf("Error: EGifPutImageDesc\n");

  unsigned char *gifpix = new unsigned char[width];

  for(int i=0;i<height;i++){
    for(int j=0;j<width;j++){
      gifpix[j] = 0;
      gifpix[j] = MaxColourMapOverlap(tmp.pixels[i*width*tmp.colourspace+j*tmp.colourspace],tmp.pixels[i*width*tmp.colourspace+j*tmp.colourspace+1],tmp.pixels[i*width*tmp.colourspace+j*tmp.colourspace+2],colmap);
    }
    
    //printf("writing line %d\n",i); fflush(stdout);
    err = EGifPutLine(giffile, gifpix, width);
    if(err!=GIF_OK)
      printf("Error: EGifPutLine\n");
  }

  EGifCloseFile(giffile);
  delete [] gifpix;

#else
    printf("No gif support compiled into this program\n");
    printf("Please install lib(un)gif (if not already done) and\n");
    printf("rebuild adding -DUSE_LIBGIF to C_DEFINES\n");
    printf("and CXX_DEFINES. Or write output file in a\n");
    printf("supported format.\n");
#endif // USE_LIBGIF 

}


#ifdef USE_LIBJPEG

#ifdef __cplusplus      //include the jpeg library as a C file.
extern "C"{
#endif
#include<jpeglib.h>

struct my_error_mgr {
  struct jpeg_error_mgr pub;    /* "public" fields */

  jmp_buf setjmp_buffer;        /* for return to caller */
};

typedef struct my_error_mgr * my_error_ptr;

METHODDEF(void)
my_error_exit (j_common_ptr cinfo)
{
  /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
  my_error_ptr myerr = (my_error_ptr) cinfo->err;

  /* Always display the message. */
  /* We could postpone this until after returning, if we chose. */
  (*cinfo->err->output_message) (cinfo);

  /* Return control to the setjmp point */
  longjmp(myerr->setjmp_buffer, 1);
}
#ifdef __cplusplus
}
#endif
#endif

void image_info::writejpg(const char *filename, int quality) const {

#ifdef USE_LIBJPEG
#ifdef _USE_DL_
  if(!have_jpeg){
    if(init_jpeg(LIBJPG_SOLIBRARY)){ 
      have_jpeg = 1;
    }
  }
  if(!have_jpeg){
    printf("JPEG library not found\n");
    printf("Please try writing to a supported format\n");
    return;
  }
#endif  

  FILE * outfile;
  if ((outfile = fopen(filename, "wb")) == NULL) {
    fprintf(stderr, "can't open %s in write_jpeg\n", filename);
    throw ImageInfoWriteJPGExc();
    return;
  }

  struct jpeg_compress_struct cinfo;

  //struct jpeg_error_mgr jerr;
  //cinfo.err = jpeg_std_error(&jerr);

  struct my_error_mgr jerr;

  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = my_error_exit;
  /* Establish the setjmp return context for my_error_exit to use. */
  if (setjmp(jerr.setjmp_buffer)) {
    /* If we get here, the JPEG code has signaled an error.
     * We need to clean up the JPEG object, close the input file, and return.
     */
    jpeg_destroy_compress(&cinfo);
    fclose(outfile);
    throw ImageInfoWriteJPGExc();
  }

  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, outfile);
  cinfo.image_width  = width;
  cinfo.image_height = height;

  if(colourspace_type==IMAGEINFO_RGBA){
    image_info tmp(*this); 
    tmp.convert_rgb();
    tmp.write(filename);
    return;
  }
  if(colourspace_type==IMAGEINFO_MONOA){
    image_info tmp(*this); 
    tmp.convert_greyscale();
    tmp.write(filename);
    return;
  }
  if(colourspace_type==IMAGEINFO_RGB){
    cinfo.in_color_space = JCS_RGB;
  }else if(colourspace_type==IMAGEINFO_YUV){
    cinfo.in_color_space = JCS_YCbCr;
  }else if(colourspace_type==IMAGEINFO_MONO){
    cinfo.in_color_space = JCS_GRAYSCALE;
  }
  cinfo.input_components = colourspace;

  jpeg_set_defaults(&cinfo);

  jpeg_set_quality(&cinfo,quality,TRUE);

  jpeg_start_compress(&cinfo, TRUE);

  JSAMPROW row_pointer[1];	/* pointer to a single row */
  int row_stride;			/* physical row width in buffer */
  
  row_stride = cinfo.image_width * colourspace;	/* JSAMPLEs per row in image_buffer */

  while (cinfo.next_scanline < cinfo.image_height) {
    row_pointer[0] = & pixels[cinfo.next_scanline * row_stride];
    jpeg_write_scanlines(&cinfo, row_pointer, 1);
  }

  jpeg_finish_compress(&cinfo); 

  fclose(outfile);
  jpeg_destroy_compress(&cinfo);
#else
    printf("No jpeg support compiled into this program\n");
    printf("Please install libjpeg (if not already done) and\n");
    printf("rebuild adding -DUSE_LIBJPEG to C_DEFINES\n");
    printf("and CXX_DEFINES. Or write to a supported format.\n");
#endif 

}

void image_info::readtif(const char *filename){
#ifdef USE_LIBTIFF
#ifdef _USE_DL_
  if(!have_tiff){
    if(init_tiff(LIBTIFF_SOLIBRARY)){ 
      have_tiff = 1;
    }
  }
  if(!have_tiff){
    printf("tiff library not found\n");
    printf("Please try reading from a supported format\n");
    return;
  }
#endif  
  TIFF* tif = TIFFOpen(filename, "r");
  uint16 spp;
  TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
  if(spp!=1&&spp!=3&&spp!=4){
    uint32 tileWidth, tileLength;
    uint32 *bc, *bc2;
    bool istiled = true, isstripped = true;
    int success = TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileWidth);
    if(!success)  istiled = false;
    success = TIFFGetField(tif, TIFFTAG_TILELENGTH, &tileLength);
    if(!success)  istiled = false;
    success = TIFFGetField(tif, TIFFTAG_STRIPBYTECOUNTS, &bc);
    if(!success)  isstripped = false;
    uint32 stripsize = bc[0];
    success = TIFFGetField(tif,  TIFFTAG_TILEBYTECOUNTS, &bc2);
    if(!success)  isstripped = false;
    uint32 tilesize = bc2[0];
    if(istiled&&isstripped) {
      printf("Seems to be tiled and stripped in readtiff.\n");
      printf("Unlikley!!, quitting ...\n");
      return;
    }
    if(istiled) {
      pixels = new unsigned char[width*height*spp];
      colourspace = spp;
      if(spp==2){
        colourspace_type = IMAGEINFO_MONOA;
      }else{
        printf("Unsupported samples per pixel in readtiff: %d\n",spp);
        delete [] pixels;
        return;
      }
      /*
      printf("Samples per pixel size of %d tiled tiff files is currently not supported by CCP4MG\n",spp);
      printf("Please convert to another tiff format, eg rgb or rgb+alpha\n");
      printf("The following is developer info to help rectify the problem\n");
      printf("Scanline size:%d\n",TIFFScanlineSize(tif));
      printf("Number of tiles:%d\n",TIFFNumberOfTiles(tif));
      printf("Number of strips:%d\n",TIFFNumberOfStrips(tif));
      printf("Strip size:%d\n",TIFFStripSize(tif));
      printf("Strip size from TIFFGetField:%d\n",stripsize);
      printf("Tile size from TIFFGetField:%d\n",tilesize);
      printf("Tile width:%d\n",tileWidth);
      printf("Tile length:%d\n",tileLength);
      printf("Seems to be tiled\n");
      */
      uint32 imageWidth, imageLength;
      uint32 tileWidth, tileLength;
      uint32 x, y;
      tdata_t buf;
      uint16 config;

      TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &imageWidth);
      TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imageLength);
      TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileWidth);
      TIFFGetField(tif, TIFFTAG_TILELENGTH, &tileLength);
      TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &config);
      buf = _TIFFmalloc(TIFFTileSize(tif));
      if (config == PLANARCONFIG_CONTIG) {
        for (y = 0; y < imageLength; y += tileLength){
          for (x = 0; x < imageWidth; x += tileWidth){
            TIFFReadTile(tif, buf, x, y, 0, 0);
            for(unsigned jj=0;jj<tileLength;jj++){
              for(unsigned ii=0;ii<tileWidth*spp;ii++){
                pixels[(y+jj)*width*spp+x+ii] = ((unsigned char *)buf)[jj*tileWidth*spp+ii];
              }
            }
          }
        }
      } else if (config == PLANARCONFIG_SEPARATE) {
        uint16 s;
        for (s = 0; s < spp; s++){
          for (y = 0; y < imageLength; y += tileLength){
            for (x = 0; x < imageWidth; x += tileWidth){
              TIFFReadTile(tif, buf, x, y, 0, s);
              for(unsigned jj=0;jj<tileLength;jj++){
                for(unsigned ii=0;ii<tileWidth;ii++){
                  pixels[(y+jj)*width*spp+x+ii+s] = ((unsigned char *)buf)[jj*tileWidth+ii];
                }
              }
            }
          }
        }
      }
      _TIFFfree(buf);
    }
    if(isstripped) {
      pixels = new unsigned char[width*height*spp];
      colourspace = spp;
      if(spp==2){
        colourspace_type = IMAGEINFO_MONOA;
      }else{
        printf("Unsupported samples per pixel in readtiff: %d\n",spp);
        delete [] pixels;
        return;
      }
      uint32 imagelength;
      tdata_t buf;
      uint32 row;
      uint16 config;

      TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imagelength);
      TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &config);
      buf = _TIFFmalloc(TIFFScanlineSize(tif));
      if (config == PLANARCONFIG_CONTIG) {
        for (row = 0; row < imagelength; row++){
          TIFFReadScanline(tif, buf, row,0);
          /* Just possible this may ignore a colourmap, not sure */
          memcpy(pixels+row*TIFFScanlineSize(tif),(unsigned char *)buf,sizeof(unsigned char)*TIFFScanlineSize(tif));
        }
      } else if (config == PLANARCONFIG_SEPARATE) {
        uint16 s;
        for (s = 0; s < spp; s++){
          for (row = 0; row < imagelength; row++){
            TIFFReadScanline(tif, buf, row, s);
            for(int ii=0;ii<width;ii++){
              pixels[spp*ii+s+width*row*spp] = ((unsigned char *)buf)[ii];
            }
          }
        }
      }
      _TIFFfree(buf);
    }
  } else {
    uint32 *raster;
    raster = (uint32*) _TIFFmalloc(width*height*sizeof(uint32));
    if(!raster){
      printf("Error reading tif image\n");
      return;
    }
    TIFFReadRGBAImage(tif, width, height, raster, 0);
    pixels = new unsigned char[width*height*4];
    colourspace = 4;
    colourspace_type = IMAGEINFO_RGBA;
    int ii=0;
    for(int i=height-1;i>=0;i--){
      for(int j=0;j<width;j++){
         uint32 pix = raster[i*width+j];
         unsigned char a = (pix>>24)&255;
         unsigned char b = (pix>>16)&255;
         unsigned char g = (pix>>8)&255;
         unsigned char r = pix&255;
         pixels[ii++] = r;
         pixels[ii++] = g;
         pixels[ii++] = b;
         pixels[ii++] = a;
      }
    }
    _TIFFfree(raster);
  }
  TIFFClose(tif);
#else
    printf("No tiff support compiled into this program\n");
    printf("Please install libtiff (if not already done) and\n");
    printf("rebuild adding -DUSE_LIBTIFF to C_DEFINES\n");
    printf("and CXX_DEFINES. Or convert input file to a\n");
    printf("supported format.\n");
#endif  
}

void image_info::writetif(const char *filename) const {
#ifdef USE_LIBTIFF
#ifdef _USE_DL_
  if(!have_tiff){
    if(init_tiff(LIBTIFF_SOLIBRARY)){ 
      have_tiff = 1;
    }
  }
  if(!have_tiff){
    printf("tiff library not found\n");
    printf("Please try reading from a supported format\n");
    return;
  }
#endif  
  image_info tmp = *this;
  tmp.convert_rgba();
  TIFF* tif = TIFFOpen(filename, "w");
  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, tmp.width);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH, tmp.height);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE,8);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL,4);
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  TIFFSetField(tif, TIFFTAG_COMPRESSION,COMPRESSION_NONE);
  uint32 *raster = (uint32*) _TIFFmalloc(tmp.width*tmp.height*sizeof(uint32));

  int ii=0;
  for(int i=0;i<tmp.height;i++){
    for(int j=0;j<tmp.width;j++){
       unsigned char r,g, b,a;
       r = tmp.pixels[ii++];
       g = tmp.pixels[ii++];
       b = tmp.pixels[ii++];
       a = tmp.pixels[ii++];
       uint32 pix = 0;
       pix |= r;
       pix |= (g<<8);
       pix |= (b<<16);
       pix |= (a<<24);
       raster[i*tmp.width+j] = pix;
    }
  }
  for (int row = 0; row < tmp.height; row++)
    TIFFWriteScanline(tif, raster+row*tmp.width, row, 0);
  _TIFFfree(raster);
  TIFFClose(tif);
#else
    printf("No tiff support compiled into this program\n");
    printf("Please install libtiff (if not already done) and\n");
    printf("rebuild adding -DUSE_LIBTIFF to C_DEFINES\n");
    printf("and CXX_DEFINES. Or try writing to a\n");
    printf("supported format.\n");
#endif  
}

void image_info::readjpg(const char *filename){
#ifdef USE_LIBJPEG
#ifdef _USE_DL_
  if(!have_jpeg){
    if(init_jpeg(LIBJPG_SOLIBRARY)){ 
      have_jpeg = 1;
    }
  }
  if(!have_jpeg){
    printf("JPEG library not found\n");
    printf("Please try reading from a supported format\n");
    return;
  }
#endif  

  struct jpeg_decompress_struct cinfo;
  //struct jpeg_error_mgr jerr;
  //cinfo.err = jpeg_std_error(&jerr);

  FILE* infile;

  if ((infile = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "can't open %s in read_jpeg\n", filename);
    return;
  }

  struct my_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = my_error_exit;
  /* Establish the setjmp return context for my_error_exit to use. */
  if (setjmp(jerr.setjmp_buffer)) {
    /* If we get here, the JPEG code has signaled an error.
     * We need to clean up the JPEG object, close the input file, and return.
     */
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);
    throw ImageInfoReadJPGExc();
  }

  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, infile);
  jpeg_read_header(&cinfo, TRUE);

  jpeg_start_decompress(&cinfo);

  int row_stride = cinfo.output_width * cinfo.output_components;

  JSAMPARRAY buffer = (*cinfo.mem->alloc_sarray)
    ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);

  width = cinfo.image_width; height = cinfo.image_height;
  pixels = new unsigned char[width*height*cinfo.output_components];

  int i = 0;

  while (cinfo.output_scanline < cinfo.output_height) {
    jpeg_read_scanlines(&cinfo, buffer, 1);
    memcpy(pixels+row_stride*i, buffer[0], sizeof(unsigned char)*cinfo.output_components*width);
    i++;
  }

  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  fclose(infile);
  colourspace = cinfo.output_components;
  if(colourspace==IMAGEINFO_RGB_SIZE)
    colourspace_type = IMAGEINFO_RGB;
  else
    colourspace_type = IMAGEINFO_MONO;
#else
    printf("No jpeg support compiled into this program\n");
    printf("Please install libjpeg (if not already done) and\n");
    printf("rebuild adding -DUSE_LIBJPEG to C_DEFINES\n");
    printf("and CXX_DEFINES. Or convert input file to a\n");
    printf("supported format.\n");
#endif 
}

void image_info::readpng(const char *filename){

#ifdef USE_LIBPNG
#ifdef _USE_DL_
  if(!have_png){
    if(init_png(LIBPNG_SOLIBRARY)){
      have_png = 1;
    }
  }
  if(!have_png){
    printf("PNG library not found\n");
    printf("Please try reading from a supported format\n");
    return;
  }
#endif  
  FILE * infile;

  if ((infile = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "can't open %s in readpng\n", filename);
    return;
  }
  
  png_byte header[8];
  try {
    my_read(header, sizeof(png_byte), 8, infile);
  }
  catch (...) {
    throw ImageInfoReadPNGExc();
  }
  
  if(png_sig_cmp(header, 0, 8)) {
    fprintf(stderr, "%s is not a PNG in readpng\n", filename);
    fclose(infile);
    throw ImageInfoReadPNGExc();
  }
  
  png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
  if(!png_ptr){
    fprintf(stderr, "Cannot create png_ptr in readpng\n");
    fclose(infile);
    throw ImageInfoReadPNGExc();
  }

  png_infop info_ptr = png_create_info_struct(png_ptr);
  if(!info_ptr){
    png_destroy_read_struct(&png_ptr,0, 0);
    fprintf(stderr, "Cannot create png_info in readpng\n");
    fclose(infile);
    throw ImageInfoReadPNGExc();
  }

  png_infop end_info = png_create_info_struct(png_ptr);
  if(!end_info){
    png_destroy_read_struct(&png_ptr, &info_ptr, 0);
    fprintf(stderr, "Cannot create end_info in readpng\n");
    fclose(infile);
    throw ImageInfoReadPNGExc();
  }

  if(setjmp(png_jmpbuf(png_ptr))){
    png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
    fprintf(stderr, "setjmp error in readpng\n");
    fclose(infile);
    throw ImageInfoReadPNGExc();
  }

  png_init_io(png_ptr, infile);
  
  png_set_sig_bytes(png_ptr,8);

  png_uint_32 png_width, png_height;
  int bit_depth, color_type, interlace_type, compression_type, filter_method;

  png_read_info(png_ptr, info_ptr);

  png_get_IHDR(png_ptr, info_ptr, &png_width, &png_height,
       &bit_depth, &color_type, &interlace_type,
       &compression_type, &filter_method);

  width = png_width;
  height = png_height;

  int channels = png_get_channels(png_ptr, info_ptr);
  if(channels==1){
    if(color_type==PNG_COLOR_TYPE_GRAY){
      colourspace = IMAGEINFO_MONO_SIZE;
      colourspace_type = IMAGEINFO_MONO;
    }else if(color_type==PNG_COLOR_TYPE_PALETTE){
      colourspace = IMAGEINFO_RGBA_SIZE;
      colourspace_type = IMAGEINFO_RGBA;
    }else{
       fprintf(stderr, "Unsupported number of channels in readpng\n");
       fclose(infile);
       throw ImageInfoReadPNGExc();
    }
  }else if(channels==3){
    colourspace = IMAGEINFO_RGB_SIZE;
    colourspace_type = IMAGEINFO_RGB;
  }else if(channels==4){
    colourspace = IMAGEINFO_RGBA_SIZE;
    colourspace_type = IMAGEINFO_RGBA;
  }else if(channels==2){
    colourspace = IMAGEINFO_MONOA_SIZE;
    colourspace_type = IMAGEINFO_MONOA;
    //png_set_strip_alpha(png_ptr);
    //fprintf(stderr, "Warning Don't know how to deal properly with GRAY_ALPHA yet, stripping alpha channel\n");
  }else{
    fprintf(stderr, "Unsupported number of channels in readpng\n");
    fclose(infile);
    throw ImageInfoReadPNGExc();
  }

  if (bit_depth == 16)
    png_set_strip_16(png_ptr);

  if (bit_depth < 8)
    png_set_packing(png_ptr);

  png_bytep *row_pointers = (png_bytep *)png_malloc(png_ptr, height*sizeof(png_bytep));
  for (int i=0; i<height; i++)
    row_pointers[i]=(png_bytep)png_malloc(png_ptr, colourspace*width*sizeof(unsigned char));

  png_read_image(png_ptr, row_pointers);

  pixels = new unsigned char[width*height*colourspace];

#if PNG_LIBPNG_VER_MAJOR >= 1 && PNG_LIBPNG_VER_MINOR >=4
  png_colorp palette = 0;
  int num_trans = 0;
  png_bytep trans_alpha = 0;
#else
#define palette png_ptr->palette
#define num_trans png_ptr->num_trans
#define trans_alpha png_ptr->trans
#endif

  if (! palette) {
     std::cout << "Null palatte" << std::endl;
     return;
  }
  for (int i=0; i<height; i++){
    if(color_type==PNG_COLOR_TYPE_PALETTE){
      for (int j=0; j<width; j++){
        pixels[i*width*colourspace+j*colourspace  ] = palette[row_pointers[i][j]].red;
        pixels[i*width*colourspace+j*colourspace+1] = palette[row_pointers[i][j]].green;
        pixels[i*width*colourspace+j*colourspace+2] = palette[row_pointers[i][j]].blue;

        if(row_pointers[i][j]<num_trans){
         pixels[i*width*colourspace+j*colourspace+3] = trans_alpha[row_pointers[i][j]];
        }else
         pixels[i*width*colourspace+j*colourspace+3] = 255;
      }
    }else{
      memcpy(pixels+i*width*colourspace,row_pointers[i],sizeof(unsigned char)*colourspace*width);
    }
  }

  png_read_end(png_ptr, end_info);
  png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);

  fclose(infile);
#else
    printf("No png support compiled into this program\n");
    printf("Please install libpng (if not already done) and\n");
    printf("rebuild adding -DUSE_LIBPNG to C_DEFINES\n");
    printf("and CXX_DEFINES. Or read from a supported format.\n");
#endif 
}

typedef struct {
   unsigned char TypeField;
   unsigned char FixHeaderField;
   unsigned char width;
   unsigned char height;
} WAPBMPHEADER;

void image_info::writewbmp(const char *filename) const {
  WAPBMPHEADER wbh;
  wbh.TypeField = 0;
  wbh.FixHeaderField = 0;

  FILE *fp;
  fp = fopen(filename,"wb");

  if(!fp) {
    printf("Error opening output file in image_info::writewbmp\n");
    return;
  }
  image_info tmp = *this;
  tmp.ScaleImage(64,64); // Every image scaled to 64*64 at moment.
  tmp.convert_greyscale();
  wbh.width = tmp.get_width();
  wbh.height = tmp.get_height();
  size_t nmemb;
  if((nmemb=fwrite( &wbh, sizeof(WAPBMPHEADER),1,fp)) != 1) { 
    printf("Couldn't write WBMP-Header %ld %lu in image_info::writewbmp\n",long(nmemb),(long)sizeof(WAPBMPHEADER));
    fclose(fp);
    return;
  }
  unsigned char *bw_pixels = new unsigned char[height*width/8];

  for(int i=0;i<wbh.height;i++){
    for(int j=0;j<wbh.width/8;j++){
      bw_pixels[i*wbh.width/8+j] = 0;
      for(int k=0;k<8;k++){
	if(pixels[i*wbh.width + j*8 + 7-k] > 128)
	  bw_pixels[i*wbh.width/8+j] += 1<<k;
      }
    }
  }
  if( fwrite(bw_pixels, 1,wbh.height*wbh.width/8,fp) != size_t(wbh.height*wbh.width/8)) {
    printf("Couldn't write WBMP pixels in image_info::writewbmp.\n");
    fclose(fp);
    return;
  }
  fclose(fp);

}

void image_info::writepng(const char *filename) const {
#ifdef USE_LIBPNG
#ifdef _USE_DL_
  if(!have_png){
    if(init_png(LIBPNG_SOLIBRARY)){
      have_png = 1;
    }
  }
  if(!have_png){
    printf("PNG library not found\n");
    printf("Please try writing to a supported format\n");
    return;
  }
#endif
  FILE * outfile;

  if ((outfile = fopen(filename, "wb")) == NULL) {
    fprintf(stderr, "can't open %s in writepng\n", filename);
    throw ImageInfoWritePNGExc();
  }

  png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, 0,0,0);
  if(!png_ptr){
    fprintf(stderr, "Error creating png_structp in writepng\n");
    fclose(outfile);
    throw ImageInfoWritePNGExc();
  }
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if(!info_ptr){
    png_destroy_write_struct(&png_ptr,0);
    fclose(outfile);
    fprintf(stderr, "Error creating png_infop in writepng\n");
    throw ImageInfoWritePNGExc();
  }
 
  if(setjmp(png_jmpbuf(png_ptr))){
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(outfile);
    fprintf(stderr, "Error calling setjmp in writepng\n");
    throw ImageInfoWritePNGExc();
   }

  png_set_filter(png_ptr, 0, PNG_FILTER_NONE);

  png_init_io(png_ptr, outfile);

  if(colourspace_type==IMAGEINFO_RGB)
    png_set_IHDR(png_ptr, info_ptr, width, height,
		 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
		 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  else if(colourspace_type==IMAGEINFO_RGBA)
    png_set_IHDR(png_ptr, info_ptr, width, height,
		 8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE,
		 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  else if(colourspace_type==IMAGEINFO_MONOA)
    png_set_IHDR(png_ptr, info_ptr, width, height,
		 8, PNG_COLOR_TYPE_GRAY_ALPHA, PNG_INTERLACE_NONE,
		 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  else if(colourspace_type==IMAGEINFO_MONO)
    png_set_IHDR(png_ptr, info_ptr, width, height,
		 8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
		 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

  png_bytep *row_pointers = (png_bytep *)png_malloc(png_ptr, height*sizeof(png_bytep));

  for (int i=0; i<height; i++){
    row_pointers[i]=(png_bytep)png_malloc(png_ptr, colourspace*width*sizeof(unsigned char));
    row_pointers[i] = pixels + i*width*colourspace;
  }

  png_set_rows(png_ptr, info_ptr, row_pointers);

  if (setjmp(png_jmpbuf(png_ptr)))
    throw ImageInfoWritePNGExc();

  png_write_png(png_ptr, info_ptr, 0, 0);

  png_write_end(png_ptr, info_ptr);
  png_destroy_write_struct(&png_ptr, &info_ptr);

  fclose(outfile);
#else
    printf("No png support compiled into this program\n");
    printf("Please install libpng (if not already done) and\n");
    printf("rebuild adding -DUSE_LIBPNG to C_DEFINES\n");
    printf("and CXX_DEFINES. Or write to a supported format.\n");
    throw ImageInfoWritePNGExc();
#endif 
}

std::vector<std::string> image_info::GetSupportedReadFormats(){
  std::vector<std::string> formats;
  formats.push_back(std::string("ppm"));
  formats.push_back(std::string("pgm"));
  formats.push_back(std::string("pbm"));
  formats.push_back(std::string("bmp"));
  formats.push_back(std::string("rgba"));
  formats.push_back(std::string("PPM"));
  formats.push_back(std::string("PGM"));
  formats.push_back(std::string("PBM"));
  formats.push_back(std::string("BMP"));
  formats.push_back(std::string("RGBA"));
  formats.push_back(std::string("xbm"));
  formats.push_back(std::string("XBM"));
  formats.push_back(std::string("xpm"));
  formats.push_back(std::string("XPM"));
#ifdef _USE_DL_
  if(!have_png){
    if(init_png(LIBPNG_SOLIBRARY)){
      have_png = 1;
    }
  }
  if(!have_jpeg){
    if(init_jpeg(LIBJPG_SOLIBRARY)){
      have_jpeg = 1;
    }
  }
  if(!have_gif){
    if(init_gif(LIBGIF_SOLIBRARY)){
      have_gif = 1;
    }
  }
  if(!have_tiff){
    if(init_tiff(LIBTIFF_SOLIBRARY)){
      have_tiff = 1;
    }
  }

  if(have_png){
    formats.push_back(std::string("png"));
    formats.push_back(std::string("PNG"));
  }
  if(have_jpeg){
    formats.push_back(std::string("jpg"));
    formats.push_back(std::string("jpeg"));
    formats.push_back(std::string("JPG"));
    formats.push_back(std::string("JPEG"));
  }
  if(have_gif){
    formats.push_back(std::string("gif"));
    formats.push_back(std::string("GIF"));
  }
  if(have_tiff){
    formats.push_back(std::string("tif"));
    formats.push_back(std::string("tiff"));
    formats.push_back(std::string("TIF"));
    formats.push_back(std::string("TIFF"));
  }
#else
#ifdef USE_LIBJPEG
  formats.push_back(std::string("jpg"));
  formats.push_back(std::string("jpeg"));
  formats.push_back(std::string("JPG"));
  formats.push_back(std::string("JPEG"));
#endif
#ifdef USE_LIBPNG
  formats.push_back(std::string("png"));
  formats.push_back(std::string("PNG"));
#endif
#ifdef USE_LIBGIF
  formats.push_back(std::string("gif"));
  formats.push_back(std::string("GIF"));
#endif
#ifdef USE_LIBTIFF
  formats.push_back(std::string("tif"));
  formats.push_back(std::string("tiff"));
  formats.push_back(std::string("TIF"));
  formats.push_back(std::string("TIFF"));
#endif
#endif
  return formats;
}

std::vector<std::string> image_info::GetSupportedWriteFormats(){
  std::vector<std::string> formats;
  formats.push_back(std::string("ppm"));
  formats.push_back(std::string("PPM"));
  formats.push_back(std::string("xbm"));
  formats.push_back(std::string("XBM"));
  formats.push_back(std::string("bmp"));
  formats.push_back(std::string("BMP"));
  formats.push_back(std::string("wbmp"));
  formats.push_back(std::string("WBMP"));
  formats.push_back(std::string("xpm"));
  formats.push_back(std::string("XPM"));
#ifdef _USE_DL_
  if(!have_png){
    if(init_png(LIBPNG_SOLIBRARY)){
      have_png = 1;
    }
  }
  if(!have_jpeg){
    if(init_jpeg(LIBJPG_SOLIBRARY)){
      have_jpeg = 1;
    }
  }
  if(!have_gif){
    if(init_gif(LIBGIF_SOLIBRARY)){
      have_gif = 1;
    }
  }
  if(!have_tiff){
    if(init_tiff(LIBTIFF_SOLIBRARY)){
      have_tiff = 1;
    }
  }
  if(have_png){
    formats.push_back(std::string("png"));
    formats.push_back(std::string("PNG"));
  }
  if(have_jpeg){
    formats.push_back(std::string("jpg"));
    formats.push_back(std::string("jpeg"));
    formats.push_back(std::string("JPG"));
    formats.push_back(std::string("JPEG"));
  }
  if(have_gif){
    formats.push_back(std::string("gif"));
    formats.push_back(std::string("GIF"));
  }
  if(have_tiff){
    formats.push_back(std::string("tif"));
    formats.push_back(std::string("tiff"));
    formats.push_back(std::string("TIF"));
    formats.push_back(std::string("TIFF"));
  }
#else
#ifdef USE_LIBJPEG
  formats.push_back(std::string("jpg"));
  formats.push_back(std::string("jpeg"));
  formats.push_back(std::string("JPG"));
  formats.push_back(std::string("JPEG"));
#endif
#ifdef USE_LIBPNG
  formats.push_back(std::string("png"));
  formats.push_back(std::string("PNG"));
#endif
#ifdef USE_LIBGIF
  formats.push_back(std::string("gif"));
  formats.push_back(std::string("GIF"));
#endif
#ifdef USE_LIBTIFF
  formats.push_back(std::string("tif"));
  formats.push_back(std::string("tiff"));
  formats.push_back(std::string("TIF"));
  formats.push_back(std::string("TIFF"));
#endif
#endif
  return formats;
}

image_info image_info::GenerateMask(void){
  image_info mask = *this;
  if(get_colourspace_type()!=IMAGEINFO_RGBA){
    return mask;
  }
  for(int i=0;i<get_height();i++){
    for(int j=3;j<get_width()*IMAGEINFO_RGBA_SIZE;j+=4){
      if(j%4==3){
	int idx = i * get_width()*IMAGEINFO_RGBA_SIZE + j;
	if(pixels[idx]==0){
	  mask.pixels[idx-3] = 255;
	  mask.pixels[idx-2] = 255;
	  mask.pixels[idx-1] = 255;
	  mask.pixels[idx]   = 255;
	  pixels[idx-3] = 0;
	  pixels[idx-2] = 0;
	  pixels[idx-1] = 0;
	  pixels[idx]   = 0;
	}else{
	  mask.pixels[idx-3] = 0;
	  mask.pixels[idx-2] = 0;
	  mask.pixels[idx-1] = 0;
	  mask.pixels[idx]   = 0;
	}
      }
    }
  }
  mask.convert_rgb();
  return mask;
}

image_info image_info::GenerateMask(unsigned char R, unsigned char G, unsigned char B){

  image_info mask = *this;
  mask.convert_rgb();
  
  for(int i=0;i<get_height();i++){
    for(int j=0;j<get_width()*IMAGEINFO_RGB_SIZE;j+=3){
      int idx = i * get_width()*IMAGEINFO_RGB_SIZE+ j;
      if(mask.pixels[idx]==R&&mask.pixels[idx+1]==G&&mask.pixels[idx+2]==B){
	mask.pixels[idx]   = 255;
	mask.pixels[idx+1] = 255;
	mask.pixels[idx+2] = 255;
      }else{
	mask.pixels[idx]   = 0;
	mask.pixels[idx+1] = 0;
	mask.pixels[idx+2] = 0;
      }
    }
  }
  return mask;
}

std::vector<std::vector<unsigned char> > image_info::GetColourMap(int ncol) const {

  image_info tmp = *this;
  tmp.convert_rgb();

  std::vector<std::vector<unsigned char> > cmap;

  /*

  std::vector<std::vector<std::vector<int> > > hist(256);
  for(int i=0;i<256;i++){
    hist[i] = std::vector<std::vector<int> >(256);
    for(int j=0;j<256;j++){
      hist[i][j] = std::vector<int>(256);
    }
  }

  for(int i=0;i<height;i++){
    for(int j=0;j<width;j++){
      unsigned char r = pixels[(i*width+j)*colourspace];
      unsigned char g = pixels[(i*width+j)*colourspace+1];
      unsigned char b = pixels[(i*width+j)*colourspace+2];
      hist[r][g][b]++;
    }
  }

  int nhist = 0;
  for(int i=0;i<256;i++){
    for(int j=0;j<256;j++){
      for(int k=0;k<256;k++){
	if(hist[i][j][k]>0){
	  cmap.push_back(std::vector<unsigned char>());
	  cmap[nhist].push_back(i);
	  cmap[nhist].push_back(j);
	  cmap[nhist++].push_back(k);
	  //std::cout << i << " " << j << " " << k << ": " << hist[i][j][k] << "\n";
	}
      }
    }
  }
  
  std::cout << "nhist: " << nhist << "\n";
  std::cout << nhist/ncol << "\n";

  if(nhist<=ncol)
    return cmap;

  if((nhist/ncol)>1){
    std::vector<std::vector<unsigned char> > cmap_new;
    int bucket_size = nhist/ncol;
    for(int i=0;i<ncol;i++){
      cmap_new.push_back(std::vector<unsigned char>());
      int start = i*bucket_size;
      int end = (i+1)*bucket_size;
      int R = 0;
      int G = 0;
      int B = 0;
      unsigned int ncols = 0;
      for(int j=start;j<end;j++){
	R += cmap[j][0];// *hist[cmap[j][0]][cmap[j][1]][cmap[j][2]];
	G += cmap[j][1];// *hist[cmap[j][0]][cmap[j][1]][cmap[j][2]];
	B += cmap[j][2];// *hist[cmap[j][0]][cmap[j][1]][cmap[j][2]];
	//std::cout << R << " " << G << " " << B << "\n";
	//ncols += hist[cmap[j][0]][cmap[j][1]][cmap[j][2]];
	ncols ++;
      }
      cmap_new.back().push_back(R/ncols);
      cmap_new.back().push_back(G/ncols);
      cmap_new.back().push_back(B/ncols);
    }
    std::cout << "cmap_new.size() " << cmap_new.size() << "\n";
    for(unsigned i=0;i<cmap_new.size();i++){
      std::cout << (int)cmap_new[i][0] << " " << (int)cmap_new[i][1] << " " << (int)cmap_new[i][2] << "\n";
    }
    return cmap_new;
  }

  return cmap;

  
  */

  /* First subdivide picture into blocks. The number of blocks
   * is at most the number of requested colours.
   */ 
  int xdiv = 0;
  int ydiv = 0;
  int xlength = width;
  int ylength = height;
  int ncol_cur = 1;

  while ((ncol_cur<ncol) && (width/((1<<xdiv)+1)>0) && (height/((1<<ydiv)+1)>0)) {
    if(xlength>ylength){
      xdiv++;
      xlength = width / (1<<xdiv);
    }else{
      ydiv++;
      ylength = height / (1<<ydiv);
    }
    ncol_cur *= 2;
  }

  //printf("ncol obtained:%d\n",ncol_cur);
  //printf("width:%d, height:%d\n",width,height);
  printf("xdiv:%d, xlength:%d, ydiv:%d, ylength:%d\n",xdiv,xlength,ydiv,ylength);

  //printf("\n");

  std::vector<std::vector<unsigned char> > average_cols;
  /* Then calculate average colour of the pixels in each block. */
  for(int i=0;i<(1<<xdiv);i++){
    int xstart = i*xlength;
    int xend   = (i+1)*xlength-1;
    /* Make sure we don't miss pixels that are beyond our nice
     * power of 2 columns.
     */
    if(i==(1<<xdiv)-1)
      xend = width - 1;
    for(int j=0;j<(1<<ydiv);j++){
      int ystart = j*ylength;
      int yend   = (j+1)*ylength-1;
      /* Make sure we don't miss pixels that are beyond our nice
       * power of 2 rows.
       */
      if(j==(1<<ydiv)-1)
        yend = height - 1;
      //printf("(%d,%d) -- (%d,%d)\n",xstart,ystart,xend,yend);
      unsigned int R=0,G=0,B=0;
      int npixel = 0;
      for(int jj=ystart;jj<=yend;jj++){
        for(int ii=xstart;ii<=xend;ii++){
          unsigned char r = tmp.pixels[jj*width*IMAGEINFO_RGB_SIZE+ii*IMAGEINFO_RGB_SIZE];
          unsigned char g = tmp.pixels[jj*width*IMAGEINFO_RGB_SIZE+ii*IMAGEINFO_RGB_SIZE+1];
          unsigned char b = tmp.pixels[jj*width*IMAGEINFO_RGB_SIZE+ii*IMAGEINFO_RGB_SIZE+2];
          R += (unsigned int)r;
          G += (unsigned int)g;
          B += (unsigned int)b;
          npixel++;
        }
      }
      if (npixel > 0) { 
	 average_cols.push_back(std::vector<unsigned char>());
	 //printf("Average R:%d G:%d B:%d\n",R/npixel,G/npixel,B/npixel);
	 average_cols.back().push_back((unsigned char)(R/npixel));
	 average_cols.back().push_back((unsigned char)(G/npixel));
	 average_cols.back().push_back((unsigned char)(B/npixel));
      }
    }
  }
  /* Of course, the colourmap may be over sized as we
   * haven't done any checking for redundancy, ie. two
   * entries may be the same.
   */
  printf("Colourmap size %ld\n",long(average_cols.size()));

  return average_cols;

}

int MaxColourMapOverlap(unsigned char R, unsigned char G, unsigned char B, const std::vector<std::vector<unsigned char> > &colmap){

  int max_i = 0;

  int max_overlap = 765;
  for(unsigned i=0;i<colmap.size();i++){
    int overlap_r = abs(colmap[i][0] - R);
    int overlap_g = abs(colmap[i][1] - G);
    int overlap_b = abs(colmap[i][2] - B);
    int overlap = overlap_r + overlap_g + overlap_b;
    if(overlap==0) return i;
    if(overlap<max_overlap){
      max_overlap = overlap;
      max_i = i;
    }
  }
  //printf("%d\n",max_overlap);
  return max_i;

}

void image_info::readxbm(const char *filename){
  FILE *fp;
  fp = fopen(filename, "rb");
  if(!fp){
    throw ImageInfoReadXBMExc();
    return;
  }
#define LINESIZE 1024
  char line[LINESIZE];
  char junk[LINESIZE];
  bool have_height=false;
  bool have_width=false;
  bool have_bits=false;
  int ival;
  while((!have_width||!have_height)&&(fgets(line,LINESIZE,fp)!=0)){
    if(sscanf(line,"#define %s %d",junk,&ival)==2){
      std::string suff(junk);
      if(suff.find("_width")!=std::string::npos){
        if(suff.substr(suff.size()-6)==std::string("_width")){
          have_width = true;
          width = ival;
        }
      }
      if(suff.find("_height")!=std::string::npos){
        if(suff.substr(suff.size()-7)==std::string("_height")){
          have_height = true;
          height = ival;
        }
      }
    }
  }
  bool shorts = false;
  while((!have_bits)&&(fgets(line,LINESIZE,fp)!=0)){
    if(sscanf(line,"static char %s[] = {",junk)==1||sscanf(line,"static unsigned char %s[] = {",junk)==1||sscanf(line,"static short %s[] = {",junk)==1)
      have_bits = true;
    if(sscanf(line,"static short %s[] = {",junk)==1)
      shorts = true;
  }

  int padded_width = (width + 7)/8;
  if((width%16)>0&&(width%16)<9&&shorts)
    padded_width++;

  pixels = new unsigned char[width*height];
  colourspace = 1;
  colourspace_type = IMAGEINFO_MONO;

  int c;
  unsigned int val;
  unsigned int p;
  for(int i=0;i<height;i++){
    for(int j=0;j<padded_width;j++){
      c = fscanf(fp,"%s",junk); 
      if(c==0||c==EOF)
        throw ImageInfoReadPPMExc();
      std::string hex(junk);
      if(shorts){
        if(hex.find(',')!=std::string::npos)
          hex = hex.substr(0,hex.length()-2);
        sscanf(hex.c_str(),"%d",&val);
      }else{
        hex = hex.substr(2,2);
        sscanf(hex.c_str(),"%x",&val);
      }
      //std::cout << i << " " << j << " " << junk << " " << val << "\n"; std::cout.flush();
      p = ~((val & 1)*255);
      if(j*8<width)
        pixels[i*width+j*8]   = p;
      p = ~((val & 2)*255/2);
      if(j*8+1<width)
        pixels[i*width+j*8+1] = p;
      p = ~((val & 4)*255/4);
      if(j*8+2<width)
        pixels[i*width+j*8+2] = p;
      p = ~((val & 8)*255/8);
      if(j*8+3<width)
        pixels[i*width+j*8+3]  = p;
      p = ~((val & 16)*255/16);
      if(j*8+4<width)
        pixels[i*width+j*8+4] = p;
      p = ~((val & 32)*255/32);
      if(j*8+5<width)
        pixels[i*width+j*8+5] = p;
      p = ~((val & 64)*255/64);
      if(j*8+6<width)
        pixels[i*width+j*8+6] = p;
      p = ~((val & 128)*255/128);
      if(j*8+7<width)
        pixels[i*width+j*8+7] = p;
    }
  }

  fclose(fp);

}

void image_info::Dither(int mode){
  double alpha = 0.4375;
  double beta = 0.0625;
  double gamma = 0.3125;
  double delta = 0.1875;
  for(int i=0;i<height;i++){
    for(int j=0;j<width;j++){
       double p = floor(pixels[i*width+j]/255.0 + 0.5);
       double e = pixels[i*width+j]/255.0  - p;
       if(i%2==0){
       if(j<width-1)
         pixels[i*width+j+1] += (unsigned char)(alpha*e*255);
       if(i<height-1&&j<width-1)
         pixels[(i+1)*width+j+1] += (unsigned char)(beta*e*255);
       if(i<height-1)
         pixels[(i+1)*width+j] += (unsigned char)(gamma*e*255);
       if(i<height-1&&j>0)
         pixels[(i+1)*width+j-1] += (unsigned char)(delta*e*255);
       }else{
       if(j>0)
         pixels[i*width+j-1] += (unsigned char)(alpha*e*255);
       if((i<(height-1) && height-1>0) && (j>0))
         pixels[(i+1)*width+j-1] += (unsigned char)(beta*e*255);
       if(i<height-1)
         pixels[(i+1)*width+j] += (unsigned char)(gamma*e*255);
       if(i<height-1&&j<width-1)
         pixels[(i+1)*width+j+1] += (unsigned char)(delta*e*255);
       }
    }
  }
}

void image_info::writexbm(const char *filename) const {
  /* Are these valid XBMs, do they contain the padding bytes? */
  image_info tmp = *this;
  tmp.convert_greyscale();
  tmp.Dither();
  FILE * outfile;
  if ((outfile = fopen(filename, "wb")) == NULL) {
    fprintf(stderr, "can't open %s in writexbm\n", filename);
    throw ImageInfoWriteXBMExc();
  }


  unsigned char *bw_pixels = new unsigned char[tmp.height*(tmp.width+7)/8];
  for(int i=0;i<height;i++){
    for(int j=0;j<width/8;j++){
      bw_pixels[i*(width+7)/8+j] = 0;
      for(int k=0;k<8;k++){
        if(tmp.pixels[i*width + j*8 + k] > 128)
          bw_pixels[i*(width+7)/8+j] += 1<<k;
      }
    }
    int leftover = tmp.width%8;
    for(int k=0;k<leftover;k++){
      if(tmp.pixels[i*width + (width/8)*8 + k] > 128)
        bw_pixels[i*(width+7)/8+ width/8] += 1<<k;
    }
  }
  size_t n_chars = strlen(filename)+1;
  char *prefix = new char[n_chars];
  strncpy(prefix,filename, n_chars); 
  prefix[strlen(filename)-4] = '\0';
  fprintf(outfile,"#define %s_width %d\n",prefix,tmp.width);
  fprintf(outfile,"#define %s_height %d\n",prefix,tmp.height);
  fprintf(outfile,"static char %s_bits[] ={\n",prefix);
  delete [] prefix;
  for(int i=0;i<height;i++){
    for(int j=0;j<(width+7)/8;j++){
      fprintf(outfile,"0x%x, ",bw_pixels[i*(width+7)/8+j]);
      if((i*((width+7)/8)+j)%12==0) fprintf(outfile,"\n");
    }
  }
  
  fprintf(outfile,"};\n");
  fclose(outfile);
  delete [] bw_pixels;
}

void image_info::readxpm(const char *filename){
#ifdef USE_LIBXPM
#ifdef _USE_DL_
  if(!have_xpm){
    if(init_xpm(LIBXPM_SOLIBRARY)){
      have_xpm = 1;
    }
  }
  if(!have_xpm){
    printf("xpm library not found\n");
    printf("Please try reading from a supported format\n");
    return;
  }
#endif
  XpmImage image; 
  XpmInfo info; 
  int ret = XpmReadFileToXpmImage((char*)filename,&image,&info);
  if(ret!=XpmSuccess){
    if(ret==XpmOpenFailed)
      printf("XPM file %s does not exist\n",filename);
    if(ret==XpmFileInvalid)
      printf("%s is not an xpm file\n",filename);
    if(ret==XpmNoMemory)
      printf("XpmNoMemory\n");
    throw ImageInfoReadXPMExc();
  }
  width = image.width;
  height = image.height;

  bool transparent=false;
  for(unsigned i=0;i<image.ncolors;i++){
    std::string rgb(image.colorTable[i].c_color);
    if(rgb==std::string("None"))
      transparent = true;
  }
  if(transparent){
    colourspace_type = IMAGEINFO_RGBA;
    colourspace = IMAGEINFO_RGBA_SIZE;
  }else{
    colourspace_type = IMAGEINFO_RGB;
    colourspace = IMAGEINFO_RGB_SIZE;
  }

  pixels = new unsigned char[width*height*colourspace];
  unsigned int r,g,b;
  for(int i=0;i<height;i++){
    for(int j=0;j<width;j++){
      std::string rgb(image.colorTable[image.data[i*width+j]].c_color);
      if(rgb==std::string("None")){
        pixels[i*width*colourspace+j*colourspace] = 0;
        pixels[i*width*colourspace+j*colourspace+1] = 0;
        pixels[i*width*colourspace+j*colourspace+2] = 0;
        pixels[i*width*colourspace+j*colourspace+3] = 0;
      }else{
        sscanf(rgb.substr(1,2).c_str(),"%x",&r);
        sscanf(rgb.substr(3,2).c_str(),"%x",&g);
        sscanf(rgb.substr(5,2).c_str(),"%x",&b);
        pixels[i*width*colourspace+j*colourspace] = r;
        pixels[i*width*colourspace+j*colourspace+1] = g;
        pixels[i*width*colourspace+j*colourspace+2] = b;
      }
    }
  }
#else
  printf("No xpm support compiled into this program\n");
  printf("Please install libxpm (if not already done) and\n");
  printf("rebuild adding -DUSE_LIBXPM to C_DEFINES\n");
  printf("and CXX_DEFINES. Or convert input file to a\n");
  printf("supported format.\n");
  throw ImageInfoReadXPMExc();
#endif // USE_LIBXPM
}

void image_info::writexpm(const char *filename) const {
#ifdef USE_LIBXPM
#ifdef _USE_DL_
  if(!have_xpm){
    if(init_xpm(LIBXPM_SOLIBRARY)){
      have_xpm = 1;
    }
  }
  if(!have_xpm){
    printf("xpm library not found\n");
    printf("Please try reading from a supported format\n");
    return;
  }
  XpmImage image; 
  XpmInfo info; 

  image_info tmp = *this;
  tmp.convert_rgb();

  image.width = tmp.width;
  image.height = tmp.height;
  /* Make colour map maximum possible required size. */
  image.colorTable = new XpmColor[tmp.width*tmp.height];
  image.data = new unsigned[tmp.width*tmp.height];
  image.cpp = 2;
  std::map<std::string,int> cmap;

  int icolour = 1;
  unsigned int r,g,b;
  for(int i=0;i<height;i++){
    for(int j=0;j<width;j++){
      r = pixels[i*width*colourspace+j*colourspace];
      g = pixels[i*width*colourspace+j*colourspace+1];
      b = pixels[i*width*colourspace+j*colourspace+2];
      char col[8];
      snprintf(col,8,"#%02x%02x%02x",r,g,b);
      col[7] = '\0';
      std::string mapentry(col);
      if(cmap[mapentry]==0){
        cmap[mapentry]=icolour;
        char colstr[3];
        snprintf(colstr,3,"%02x",icolour-1);
        colstr[2] = '\0';
        image.colorTable[cmap[mapentry]-1].c_color = (char *)(mapentry.c_str());
        image.colorTable[cmap[mapentry]-1].string = new char[3];
        strncat(image.colorTable[cmap[mapentry]-1].string,colstr,2);
        icolour++;
      }
      image.data[i*width+j] = cmap[mapentry]-1;
    }
  }
  image.ncolors = cmap.size();

  info.valuemask = XpmColorTable|XpmComments;
  info.hints_cmt = "Hints";
  info.colors_cmt = "Colours";
  info.pixels_cmt = "Pixels";
  info.x_hotspot = 0;
  info.y_hotspot = 0;
  info.nextensions = 0;
  info.extensions = 0;

  int ret = XpmWriteFileFromXpmImage((char*)filename,&image,&info);
  if(ret!=XpmSuccess){
    if(ret==XpmOpenFailed)
      printf("Cannot open XPM file %s for writing\n",filename);
    if(ret==XpmNoMemory)
      printf("XpmNoMemory\n");
    throw ImageInfoReadXPMExc();
  }
  delete [] image.colorTable;
  delete [] image.data;
#endif
#else
  printf("No xpm support compiled into this program\n");
  printf("Please install libxpm (if not already done) and\n");
  printf("rebuild adding -DUSE_LIBXPM to C_DEFINES\n");
  printf("and CXX_DEFINES. Or write output to a\n");
  printf("supported format.\n");
  throw ImageInfoWriteXPMExc();
#endif // USE_LIBGIF
}

int image_info::Overlay(const image_info &overlay_in, int x, int y){
  image_info overlay = overlay_in;

  // Overlay mono
  if(overlay.get_colourspace_type()==IMAGEINFO_MONO){
    if(get_colourspace_type()==IMAGEINFO_MONOA)
      overlay.convert_greyscalea();
    else if(get_colourspace_type()==IMAGEINFO_RGB)
      overlay.convert_rgb();
    else if(get_colourspace_type()==IMAGEINFO_RGBA)
      overlay.convert_rgba();
    else if(get_colourspace_type()!=IMAGEINFO_MONO)
      convert_greyscale();
  }
  // Overlay monoa
  if(overlay.get_colourspace_type()==IMAGEINFO_MONOA){
    if(get_colourspace_type()==IMAGEINFO_MONO)
      convert_greyscalea();
    else if(get_colourspace_type()==IMAGEINFO_RGB){
      overlay.convert_rgba();
      convert_rgba();
    }else if(get_colourspace_type()==IMAGEINFO_RGBA)
      overlay.convert_rgba();
    else if(get_colourspace_type()!=IMAGEINFO_MONOA)
      convert_greyscalea();
  }
  // Overlay rgb
  if(overlay.get_colourspace_type()==IMAGEINFO_RGB){
    if(get_colourspace_type()==IMAGEINFO_MONO)
      convert_rgb();
    else if(get_colourspace_type()==IMAGEINFO_MONOA){
      convert_rgba();
      overlay.convert_rgba();
    }else if(get_colourspace_type()==IMAGEINFO_RGBA)
      overlay.convert_rgba();
    else if(get_colourspace_type()!=IMAGEINFO_RGB)
      convert_rgb();
  }
  // Overlay rgba
  if(overlay.get_colourspace_type()==IMAGEINFO_RGBA)
     convert_rgba();

  unsigned char *overlay_pix = overlay.get_pixels();
  for(int i=y;i<y+overlay.get_height()&&i<height;i++){
    for(int j=x;j<x+overlay.get_width()&&j<width;j++){
      if(colourspace_type==IMAGEINFO_MONO||colourspace_type==IMAGEINFO_RGB){
        for(int k=0;k<colourspace;k++){
          pixels[i*width*colourspace+j*colourspace+k] = overlay_pix[(i-y)*overlay.get_width()*colourspace+(j-x)*colourspace+k];
        }
      }else if(colourspace_type==IMAGEINFO_MONOA){
        float alpha = float(overlay_pix[(i-y)*overlay.get_width()*colourspace+(j-x)*colourspace+1])/255.0;
        float pix = (alpha*(overlay_pix[(i-y)*overlay.get_width()*colourspace+(j-x)*colourspace])+(1.0f-alpha)*(pixels[i*width*colourspace+j*colourspace]));
        if(pix>255.0) pix = 255.0;
        pixels[i*width*colourspace+j*colourspace]   = (unsigned char)pix;
      }else if(colourspace_type==IMAGEINFO_RGBA){
        float alpha = float(overlay_pix[(i-y)*overlay.get_width()*colourspace+(j-x)*colourspace+3])/255.0;
        float pix = (alpha*(overlay_pix[(i-y)*overlay.get_width()*colourspace+(j-x)*colourspace])+(1.0f-alpha)*(pixels[i*width*colourspace+j*colourspace]));
        if(pix>255.0) pix = 255.0;
        pixels[i*width*colourspace+j*colourspace]   = (unsigned char)pix;
        pix = (alpha*(overlay_pix[(i-y)*overlay.get_width()*colourspace+(j-x)*colourspace+1])+(1.0f-alpha)*(pixels[i*width*colourspace+j*colourspace+1]));
        if(pix>255.0) pix = 255.0;
        pixels[i*width*colourspace+j*colourspace+1] = (unsigned char)pix;
        pix = (alpha*(overlay_pix[(i-y)*overlay.get_width()*colourspace+(j-x)*colourspace+2])+(1.0f-alpha)*(pixels[i*width*colourspace+j*colourspace+2]));
        if(pix>255.0) pix = 255.0;
        pixels[i*width*colourspace+j*colourspace+2] = (unsigned char)pix;
      }
    }
  }
  return 0;
}

class sort_float_int {
  public:
    int operator()(const std::pair<float,int>&p1, const std::pair<float,int>&p2) const {
      return p1.first < p2.first;
    }
};

void image_info::Kuwahara(void){
  unsigned char *new_pixels = new unsigned char[(width)*(height)*colourspace];

  float mult_mono[] = {1.0};
  float mult_rgb[] = {0.299f,0.587f,0.114f};
  float *mult;

  if(colourspace==4)
     mult = mult_rgb;
  if(colourspace==1)
     mult = mult_mono;

  std::vector<std::pair<float,int> > area_map(4);
  std::vector<float> avge_pixel(colourspace);

  for(int ii=0;ii<height;ii++){
    for(int jj=0;jj<width;jj++){
      float a1_p1 = 0.0f;
      float a1_p2 = 0.0f;
      float a1_p3 = 0.0f;
      float a1_p4 = 0.0f;
      float a1_p5 = 0.0f;
      float a1_p6 = 0.0f;
      for(int kk=0;kk<colourspace&&kk<3;kk++){
        if(ii>1&&jj>1) a1_p1 += mult[kk]*(float)pixels[(ii-2)*(width)*colourspace + (jj-2)*colourspace + kk];
        if(ii>1&&jj>0) a1_p2 += mult[kk]*(float)pixels[(ii-2)*(width)*colourspace + (jj-1)*colourspace + kk];
        if(ii>1) a1_p3 += mult[kk]*(float)pixels[(ii-2)*(width)*colourspace + (jj)*colourspace + kk];
        if(ii>0&&jj>1) a1_p4 += mult[kk]*(float)pixels[(ii-1)*(width)*colourspace + (jj-2)*colourspace + kk];
        if(ii>0&&jj>0) a1_p5 += mult[kk]*(float)pixels[(ii-1)*(width)*colourspace + (jj-1)*colourspace + kk];
        if(ii>0) a1_p6 += mult[kk]*(float)pixels[(ii-1)*(width)*colourspace + (jj)*colourspace + kk];
      }
      float a1_mean = (a1_p1 + a1_p2 + a1_p3 + a1_p4 + a1_p5 + a1_p6)/6.0f;
      float a1_var = (a1_p1-a1_mean) * (a1_p1-a1_mean)
                 + (a1_p2-a1_mean) * (a1_p2-a1_mean)
                 + (a1_p3-a1_mean) * (a1_p3-a1_mean)
                 + (a1_p4-a1_mean) * (a1_p4-a1_mean)
                 + (a1_p5-a1_mean) * (a1_p5-a1_mean)
                 + (a1_p6-a1_mean) * (a1_p6-a1_mean);

      float a2_p1 = 0.0f;
      float a2_p2 = 0.0f;
      float a2_p3 = 0.0f;
      float a2_p4 = 0.0f;
      float a2_p5 = 0.0f;
      float a2_p6 = 0.0f;
      for(int kk=0;kk<colourspace&&kk<3;kk++){
        if(ii>1&&jj<width-1) a2_p1 += mult[kk]*(float)pixels[(ii-2)*(width)*colourspace + (jj+1)*colourspace + kk];
        if(ii>1&&jj<width-2) a2_p2 += mult[kk]*(float)pixels[(ii-2)*(width)*colourspace + (jj+2)*colourspace + kk];
        if(ii>0&&jj<width-1) a2_p3 += mult[kk]*(float)pixels[(ii-1)*(width)*colourspace + (jj+1)*colourspace + kk];
        if(ii>0&&jj<width-2) a2_p4 += mult[kk]*(float)pixels[(ii-1)*(width)*colourspace + (jj+2)*colourspace + kk];
        if(jj<width-1) a2_p5 += mult[kk]*(float)pixels[(ii)*(width)*colourspace + (jj+1)*colourspace + kk];
        if(jj<width-2) a2_p6 += mult[kk]*(float)pixels[(ii)*(width)*colourspace + (jj+2)*colourspace + kk];
      }
      float a2_mean = (a2_p1 + a2_p2 + a2_p3 + a2_p4 + a2_p5 + a2_p6)/6.0f;
      float a2_var = (a2_p1-a2_mean) * (a2_p1-a2_mean)
                 + (a2_p2-a2_mean) * (a2_p2-a2_mean)
                 + (a2_p3-a2_mean) * (a2_p3-a2_mean)
                 + (a2_p4-a2_mean) * (a2_p4-a2_mean)
                 + (a2_p5-a2_mean) * (a2_p5-a2_mean)
                 + (a2_p6-a2_mean) * (a2_p6-a2_mean);

      float a3_p1 = 0.0f;
      float a3_p2 = 0.0f;
      float a3_p3 = 0.0f;
      float a3_p4 = 0.0f;
      float a3_p5 = 0.0f;
      float a3_p6 = 0.0f;
      for(int kk=0;kk<colourspace&&kk<3;kk++){
        if(ii<height-1) a3_p1 += mult[kk]*(float)pixels[(ii+1)*(width)*colourspace + (jj)*colourspace + kk];
        if(ii<height-1&&jj<width-1) a3_p2 += mult[kk]*(float)pixels[(ii+1)*(width)*colourspace + (jj+1)*colourspace + kk];
        if(ii<height-1&&jj<width-2) a3_p3 += mult[kk]*(float)pixels[(ii+1)*(width)*colourspace + (jj+2)*colourspace + kk];
        if(ii<height-2) a3_p4 += mult[kk]*(float)pixels[(ii+2)*(width)*colourspace + (jj)*colourspace + kk];
        if(ii<height-2&&jj<width-1) a3_p5 += mult[kk]*(float)pixels[(ii+2)*(width)*colourspace + (jj+1)*colourspace + kk];
        if(ii<height-2&&jj<width-2) a3_p6 += mult[kk]*(float)pixels[(ii+2)*(width)*colourspace + (jj+2)*colourspace + kk];
      }
      float a3_mean = (a3_p1 + a3_p2 + a3_p3 + a3_p4 + a3_p5 + a3_p6)/6.0f;
      float a3_var = (a3_p1-a3_mean) * (a3_p1-a3_mean)
                 + (a3_p2-a3_mean) * (a3_p2-a3_mean)
                 + (a3_p3-a3_mean) * (a3_p3-a3_mean)
                 + (a3_p4-a3_mean) * (a3_p4-a3_mean)
                 + (a3_p5-a3_mean) * (a3_p5-a3_mean)
                 + (a3_p6-a3_mean) * (a3_p6-a3_mean);

      float a4_p1 = 0.0f;
      float a4_p2 = 0.0f;
      float a4_p3 = 0.0f;
      float a4_p4 = 0.0f;
      float a4_p5 = 0.0f;
      float a4_p6 = 0.0f;
      for(int kk=0;kk<colourspace&&kk<3;kk++){
        if(jj>1) a4_p1 += mult[kk]*(float)pixels[(ii)*(width)*colourspace + (jj-2)*colourspace + kk];
        if(jj>0) a4_p2 += mult[kk]*(float)pixels[(ii)*(width)*colourspace + (jj-1)*colourspace + kk];
        if(ii<height-1&&jj>0) a4_p3 += mult[kk]*(float)pixels[(ii+1)*(width)*colourspace + (jj-2)*colourspace + kk];
        if(ii<height-1&&jj>1) a4_p4 += mult[kk]*(float)pixels[(ii+1)*(width)*colourspace + (jj-1)*colourspace + kk];
        if(ii<height-2&&jj>0) a4_p5 += mult[kk]*(float)pixels[(ii+2)*(width)*colourspace + (jj-2)*colourspace + kk];
        if(ii<height-2&&jj>1) a4_p6 += mult[kk]*(float)pixels[(ii+2)*(width)*colourspace + (jj-1)*colourspace + kk];
      }
      float a4_mean = (a4_p1 + a4_p2 + a4_p3 + a4_p4 + a4_p5 + a4_p6)/6.0f;
      float a4_var = (a4_p1-a4_mean) * (a4_p1-a4_mean)
                 + (a4_p2-a4_mean) * (a4_p2-a4_mean)
                 + (a4_p3-a4_mean) * (a4_p3-a4_mean)
                 + (a4_p4-a4_mean) * (a4_p4-a4_mean)
                 + (a4_p5-a4_mean) * (a4_p5-a4_mean)
                 + (a4_p6-a4_mean) * (a4_p6-a4_mean);
      area_map[0] = std::pair<float,int>(a1_var,0);
      area_map[1] = std::pair<float,int>(a2_var,1);
      area_map[2] = std::pair<float,int>(a3_var,2);
      area_map[3] = std::pair<float,int>(a4_var,3);
      std::sort(area_map.begin(),area_map.end(),sort_float_int());
      
      if(area_map[0].second==0){
        for(int kk=0;kk<colourspace&&kk<3;kk++){
          if(ii>1&&jj>1) avge_pixel[kk] = (float)pixels[(ii-2)*(width)*colourspace + (jj-2)*colourspace + kk];
          if(ii>1&&jj>0) avge_pixel[kk] += (float)pixels[(ii-2)*(width)*colourspace + (jj-1)*colourspace + kk];
          if(ii>1)       avge_pixel[kk] += (float)pixels[(ii-2)*(width)*colourspace + (jj)*colourspace + kk];
          if(ii>0&&jj>1) avge_pixel[kk] += (float)pixels[(ii-1)*(width)*colourspace + (jj-2)*colourspace + kk];
          if(ii>0&&jj>0) avge_pixel[kk] += (float)pixels[(ii-1)*(width)*colourspace + (jj-1)*colourspace + kk];
          if(ii>0)       avge_pixel[kk] += (float)pixels[(ii-1)*(width)*colourspace + (jj)*colourspace + kk];
        }
      }else if(area_map[0].second==1){
        for(int kk=0;kk<colourspace&&kk<3;kk++){
          if(ii>1&&jj<width-1) avge_pixel[kk] = (float)pixels[(ii-2)*(width)*colourspace + (jj+1)*colourspace + kk];
          if(ii>1&&jj<width-2) avge_pixel[kk] += (float)pixels[(ii-2)*(width)*colourspace + (jj+2)*colourspace + kk];
          if(ii>0&&jj<width-1) avge_pixel[kk] += (float)pixels[(ii-1)*(width)*colourspace + (jj+1)*colourspace + kk];
          if(ii>0&&jj<width-2) avge_pixel[kk] += (float)pixels[(ii-1)*(width)*colourspace + (jj+2)*colourspace + kk];
          if(jj<width-1)       avge_pixel[kk] += (float)pixels[(ii)*(width)*colourspace + (jj+1)*colourspace + kk];
          if(jj<width-2)       avge_pixel[kk] += (float)pixels[(ii)*(width)*colourspace + (jj+2)*colourspace + kk];
        }
      }else if(area_map[0].second==2){
        for(int kk=0;kk<colourspace&&kk<3;kk++){
          if(ii<height-1)             avge_pixel[kk] = (float)pixels[(ii+1)*(width)*colourspace + (jj)*colourspace + kk];
          if(ii<height-1&&jj<width-1) avge_pixel[kk] += (float)pixels[(ii+1)*(width)*colourspace + (jj+1)*colourspace + kk];
          if(ii<height-1&&jj<width-2) avge_pixel[kk] += (float)pixels[(ii+1)*(width)*colourspace + (jj+2)*colourspace + kk];
          if(ii<height-2)             avge_pixel[kk] += (float)pixels[(ii+2)*(width)*colourspace + (jj)*colourspace + kk];
          if(ii<height-2&&jj<width-1) avge_pixel[kk] += (float)pixels[(ii+2)*(width)*colourspace + (jj+1)*colourspace + kk];
          if(ii<height-2&&jj<width-2) avge_pixel[kk] += (float)pixels[(ii+2)*(width)*colourspace + (jj+2)*colourspace + kk];
        }
      }else if(area_map[0].second==3){
        for(int kk=0;kk<colourspace&&kk<3;kk++){
          if(jj>1)              avge_pixel[kk] = (float)pixels[(ii)*(width)*colourspace + (jj-2)*colourspace + kk];
          if(jj>0)              avge_pixel[kk] += (float)pixels[(ii)*(width)*colourspace + (jj-1)*colourspace + kk];
          if(ii<height-1&&jj>0) avge_pixel[kk] += (float)pixels[(ii+1)*(width)*colourspace + (jj-2)*colourspace + kk];
          if(ii<height-1&&jj>1) avge_pixel[kk] += (float)pixels[(ii+1)*(width)*colourspace + (jj-1)*colourspace + kk];
          if(ii<height-2&&jj>0) avge_pixel[kk] += (float)pixels[(ii+2)*(width)*colourspace + (jj-2)*colourspace + kk];
          if(ii<height-2&&jj>1) avge_pixel[kk] += (float)pixels[(ii+2)*(width)*colourspace + (jj-1)*colourspace + kk];
        }
      }

      for(int kk=0;kk<colourspace&&kk<3;kk++){
        new_pixels[ii*(width)*colourspace + jj*colourspace + kk] = (unsigned char)(avge_pixel[kk]/6.0f);
      }
      if(colourspace==4)
        new_pixels[ii*(width)*colourspace + jj*colourspace + 3] = pixels[ii*(width)*colourspace + jj*colourspace + 3];
    }
  }

  pixels = new_pixels;

}
