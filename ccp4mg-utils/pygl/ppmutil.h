/*
     pygl/ppmutil.h: CCP4MG Molecular Graphics Program
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

enum { IMAGEINFO_MONO, IMAGEINFO_MONOA, IMAGEINFO_RGB, IMAGEINFO_RGBA, IMAGEINFO_YUV, IMAGEINFO_CMAP, IMAGEINFO_CMYK, IMAGEINFO_ARGB }; 
#define IMAGEINFO_MONO_SIZE 1
#define IMAGEINFO_MONOA_SIZE 2
#define IMAGEINFO_RGB_SIZE 3
#define IMAGEINFO_RGBA_SIZE 4
#define IMAGEINFO_ARGB_SIZE 4
#define IMAGEINFO_YUV_SIZE 3
#define IMAGEINFO_CMAP_SIZE 1
#define IMAGEINFO_CMYK_SIZE 4

class image_info {
 protected:
  int get_next_ppm_token(FILE *fp);
  int width;
  int height;
  int colourspace;
  int colourspace_type;
  unsigned char *pixels;
  unsigned char *cmapped_pixels;
  unsigned char *cmap;
  int ncmapped_colours;
  std::vector<std::vector<unsigned char> > colour_map;
 public:
  image_info();
  image_info(const image_info &iiinfo);
  image_info(const std::vector<image_info> &patches, const std::vector<std::vector<int> > &pattern);
  image_info(const char *filename);
  image_info(int width_in, int height_in, unsigned char *pixels_in, int colourspace_type_in=IMAGEINFO_RGB);
  void set_bitmap_data(int width_in, int height_in, unsigned char *bm);
  virtual ~image_info();
  virtual void invert(void);
  virtual void invert_colourmap(void);
  virtual void convert_rgb(void);
  virtual void convert_rgba(void);
  virtual void convert_cmyk(void); /* Not to be used publically really. But we'll keep it here for future.*/
  virtual void convert_yuv(void);
  virtual void convert_greyscale(void);
  virtual void convert_greyscalea(void);
  virtual void convert_colourspace(int colourspace_type_in);
  virtual int read(const char *filename);
  virtual int write(const char *filename, int flags=100) const;
  virtual int get_width(void) const {return width;};
  virtual int get_height(void) const {return height;};
  virtual int get_colourspace(void) const {return colourspace;};
  virtual int get_colourspace_type(void) const {return colourspace_type;};
  virtual unsigned char* get_pixels(void) const {return pixels;};
  virtual unsigned char* get_cmapped_pixels(void) const {return cmapped_pixels;};
  virtual unsigned char* get_cmap(void) const {return cmap;};
  virtual int get_ncolours(void) const {return ncmapped_colours;};
  virtual void set_width(int width_in) {width = width_in;};
  virtual void set_height(int height_in) {height = height_in;};
  virtual void set_pixels(unsigned char* pixels_in) {pixels = pixels_in;};
  virtual void set_colourspace(int colourspace_in) {colourspace = colourspace_in;};
  virtual void set_colourspace_type(int colourspace_type_in);
  virtual void ScaleImage(int x, int y);
  virtual void Rotate();
  virtual void Dither(int mode=0);
  image_info const& operator=(image_info const &b);
  virtual void copy(int width_in, int height_in, unsigned char *pixels_in, int colourspace_type_in, unsigned char *cmapped_pixels_in, unsigned char *cmap_in, int ncolours_in);
  static std::vector<std::string> GetSupportedReadFormats();
  static std::vector<std::string> GetSupportedWriteFormats();
  image_info GenerateMask(void);
  image_info GenerateMask(unsigned char R, unsigned char G, unsigned char B);
  virtual std::vector<std::vector<unsigned char> > GetColourMap(int ncol) const;
  virtual int Overlay(const image_info &overlay, int x, int y);
  void Kuwahara();
  image_info_yuv_t getyuv(bool subsample=true) const;
  virtual void writeppm(int fildes) const ;
  virtual void writeyuv(int fildes) const ;
  virtual void writebmp(int fildes) const ;
  virtual void writejpg(int fildes, int quality=100) const ;
  virtual void writepng(int fildes) const ;
  virtual void writewbmp(int fildes) const ;
  virtual void writexbm(int fildes, const char *name=0) const ;
  virtual void writetif(int fildes, const char *filename, int flags=100) const ;
  virtual void writeppm(FILE *fp) const ;
  virtual void writeyuv(FILE *fp) const ;
  virtual void writebmp(FILE *fp) const ;
  virtual void writejpg(FILE *fp, int quality=100) const ;
  virtual void writepng(FILE *fp) const ;
  virtual void writewbmp(FILE *fp) const ;
  virtual void writexbm(FILE *fp, const char *name=0) const ;
  virtual void writeppm(const char *filename) const ;
  virtual void writeyuv(const char *filename) const ;
  virtual void writebmp(const char *filename) const ;
  virtual void writejpg(const char *filename, int quality=100) const ;
  virtual void writepng(const char *filename) const ;
  virtual void writetif(const char *filename, int flags=100) const ;
  virtual void writegif(const char *filename) const ;
  virtual void writewbmp(const char *filename) const ;
  virtual void writexbm(const char *filename) const ;
  virtual void writexpm(const char *filename) const ;
  virtual void readppm(const char *filename);
  virtual void readbmp(const char *filename);
  virtual void readjpg(const char *filename);
  virtual void readpng(const char *filename);
  virtual void readtif(const char *filename);
  virtual void readgif(const char *filename);
  virtual void readrgba(const char *filename);// RAW RGBA
  virtual void readxbm(const char *filename);
  virtual void readxpm(const char *filename);
};

#endif
