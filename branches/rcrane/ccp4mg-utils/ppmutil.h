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


#ifndef __CCP4_PPM__
#define __CCP4_PPM__

typedef struct _image_info_yuv_t_ {
  unsigned int width, height, _dum;
  unsigned char *y;
  unsigned char *u;
  unsigned char *v;
} image_info_yuv_t ;

#include <stdio.h>
#include <vector>
#include <string>

enum { IMAGEINFO_MONO, IMAGEINFO_MONOA, IMAGEINFO_RGB, IMAGEINFO_RGBA, IMAGEINFO_YUV, IMAGEINFO_CMAP }; 
#define IMAGEINFO_MONO_SIZE 1
#define IMAGEINFO_MONOA_SIZE 2
#define IMAGEINFO_RGB_SIZE 3
#define IMAGEINFO_RGBA_SIZE 4
#define IMAGEINFO_YUV_SIZE 3
#define IMAGEINFO_CMAP_SIZE 1

class image_info {
  int get_next_ppm_token(FILE *fp);
  int width;
  int height;
  int colourspace;
  int colourspace_type;
  unsigned char *pixels;
  std::vector<std::vector<unsigned char> > colour_map;
  void writeppm(const char *filename) const ;
  void writeyuv(const char *filename) const ;
  void readppm(const char *filename);
  void readbmp(const char *filename);
  void writebmp(const char *filename) const ;
  void writejpg(const char *filename, int quality=100) const ;
  void readjpg(const char *filename);
  void writepng(const char *filename) const ;
  void readpng(const char *filename);
  void readtif(const char *filename);
  void writetif(const char *filename) const ;
  void readgif(const char *filename);
  void writegif(const char *filename) const ;
  void readrgba(const char *filename);// RAW RGBA
  void writewbmp(const char *filename) const ;
  void readxbm(const char *filename);
  void writexbm(const char *filename) const ;
  void readxpm(const char *filename);
  void writexpm(const char *filename) const ;
 public:
  image_info();
  image_info(const image_info &iiinfo);
  image_info(const std::vector<image_info> &patches, const std::vector<std::vector<int> > &pattern);
  image_info(const char *filename);
  image_info(int width_in, int height_in, unsigned char *pixels_in, int colourspace_type_in=IMAGEINFO_RGB);
  void set_bitmap_data(int width_in, int height_in, unsigned char *bm);
  ~image_info();
  void invert(void);
  void invert_colourmap(void);
  void convert_rgb(void);
  void convert_rgba(void);
  void convert_yuv(void);
  void convert_greyscale(void);
  void convert_greyscalea(void);
  void convert_colourspace(int colourspace_type_in);
  int read(const char *filename);
  int write(const char *filename, int quality=100) const;
  int get_width(void) const {return width;};
  int get_height(void) const {return height;};
  int get_colourspace(void) const {return colourspace;};
  int get_colourspace_type(void) const {return colourspace_type;};
  unsigned char* get_pixels(void) const {return pixels;};
  void set_width(int width_in) {width = width_in;};
  void set_height(int height_in) {height = height_in;};
  void set_pixels(unsigned char* pixels_in) {pixels = pixels_in;};
  void set_colourspace(int colourspace_in) {colourspace = colourspace_in;};
  void set_colourspace_type(int colourspace_type_in);
  void ScaleImage(int x, int y);
  void Rotate();
  void Dither(int mode=0);
  image_info const& operator=(image_info const &b);
  void copy(int width_in, int height_in, unsigned char *pixels_in, int colourspace_type_in);
  static std::vector<std::string> GetSupportedReadFormats();
  static std::vector<std::string> GetSupportedWriteFormats();
  image_info GenerateMask(void);
  image_info GenerateMask(unsigned char R, unsigned char G, unsigned char B);
  std::vector<std::vector<unsigned char> > GetColourMap(int ncol) const;
  int Overlay(const image_info &overlay, int x, int y);
  void Kuwahara();
  image_info_yuv_t getyuv(bool subsample=true) const;
};

#endif
