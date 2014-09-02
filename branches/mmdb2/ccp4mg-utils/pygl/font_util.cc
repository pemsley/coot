/*
     pygl/font_util.cc: CCP4MG Molecular Graphics Program
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


#include <stdio.h>
#include <iostream>
#include "font_util.h"
#include <ppmutil.h>
#include <string>
#include <vector>
#include <math.h>
#include <string.h>

image_info OverlayTextOnImage(const image_info &iinfo_in, const MGFontInfo &finfo, int x, int y, const std::string &str, const std::string &fg, const std::string &bg, const bool &dum){
  image_info iinfo = iinfo_in;
  if(!finfo.isValid()){
    std::cout << "Font not valid when attempting to overlay text on image\n";
    return iinfo;
  }
  image_info text = MakeStringImage(finfo,str,fg,bg,false);
  unsigned char *pixels = text.get_pixels();
  for(int i=0;i<text.get_height();i++){
    for(int j=0;j<text.get_width();j++){
       unsigned char r = pixels[(i*text.get_width()+j)*IMAGEINFO_RGBA_SIZE];
       unsigned char g = pixels[(i*text.get_width()+j)*IMAGEINFO_RGBA_SIZE+1];
       unsigned char b = pixels[(i*text.get_width()+j)*IMAGEINFO_RGBA_SIZE+2];
       unsigned char a = pixels[(i*text.get_width()+j)*IMAGEINFO_RGBA_SIZE+3];
       if(a>0){
         r = (unsigned char)(float(r) * 255.0/float(a));
         g = (unsigned char)(float(g) * 255.0/float(a));
         b = (unsigned char)(float(b) * 255.0/float(a));
         pixels[(i*text.get_width()+j)*IMAGEINFO_RGBA_SIZE] = r;
         pixels[(i*text.get_width()+j)*IMAGEINFO_RGBA_SIZE+1] = g;
         pixels[(i*text.get_width()+j)*IMAGEINFO_RGBA_SIZE+2] = b;
       }
    }
  }
  text.set_pixels(pixels);
  text.invert();
  iinfo.Overlay(text,x,y);
  return iinfo;
}

int WriteStringImage(const MGFontInfo &finfo, const std::string &str, const std::string &filename, const std::string &fg, const std::string &bg, const bool &scale_alpha){
  image_info im = MakeStringImage(finfo,str,fg,bg,scale_alpha);
  im.invert();
  return im.write(filename.c_str());
}

image_info MakeExampleStringImage(const MGFontInfo &finfo, const std::string &fg, const std::string &bg, const bool &scale_alpha){
  return MakeStringImage(finfo,finfo.Family(),fg,bg,scale_alpha);
}

image_info MakeStringImage(const MGFontInfo &finfo, const std::string &str, const std::string &fg, const std::string &bg, const bool &scale_alpha){
  image_info iinfo;
  if(!finfo.isValid()){
    std::cout << "Font not valid when attempting to MakeStringImage\n";
    return iinfo;
  }

  std::vector<image_info> letters;
  for(unsigned i=0;i<str.length();i++){
    letters.push_back(image_info(finfo.FullBitMapWidth(str[i]),finfo.Height(str[i]),finfo.FullBitMap(str[i]),IMAGEINFO_MONO));
  }
  int have_kerning = 0;
  if(finfo.HaveKerning()){
    have_kerning = 1;
  }
  int kern = 0;

  int width = 0;

#ifdef USE_FREETYPE
  FT_Vector kerning = {0,0};
#endif
  int this_advx = int(ceil(double(finfo.AdvanceX(str[0]))));
  int this_width = int(ceil(double(finfo.Width(str[0]))));
  if(this_advx>this_width)
    width += this_advx;
  else
    width += this_width;
  
  for(unsigned int i = 1; i< letters.size()-1; i++){
    width += int(ceil(finfo.AdvanceX(str[i])));
  }

  this_advx = int(ceil(double(finfo.AdvanceX(str[str.length()-1]))));
  this_width = int(ceil(double(finfo.Width(str[str.length()-1]))));
  if(letters.size()>1){
    if(this_advx>this_width)
      width += this_advx;
    else
      width += this_width;
  }
   width += int(ceil(finfo.LBearing(str[0])))+int(ceil(finfo.LBearing(str[str.length()-1])));

  width += 2*finfo.Width('W');
  width =  (width+7) & ~7;

  int max_height = 0;
  int h;

  int max_descent = 0;
  int d = 0;
  for(unsigned int i = 0; i< letters.size(); i++){
    if((d=int(ceil(finfo.Descent(str[i]))))>max_descent)
      max_descent = d;
  }
  for(unsigned int i = 0; i< letters.size(); i++){
    if((h=letters[i].get_height()+max_descent)>max_height)
      max_height = h;
  }
    
  int height = max_height;
  unsigned char *pixels = new unsigned char[width*height*2];
  memset(pixels,0,2*width*height*sizeof(unsigned char));
  int max_address=2*width*height;

  int current_advance = -int(ceil(finfo.LBearing(str[0])));;
  for(unsigned i=0;i<letters.size();i++){
    unsigned char *letter_pixels = letters[i].get_pixels();
    int bearing = 0;
    bearing = int(ceil(finfo.LBearing(str[i])));
#ifdef USE_FREETYPE
    if(i>0&&have_kerning){
      kerning = finfo.Kerning(str[i-1],str[i]);
      //std::cout << "kerning of " << str[i-1] << ", " << str[i] << ": " << kerning.x/64.0 << "\n";
      kern = int(ceil(kerning.x/64.0));
    }
#endif

    for(int row = 0; row<letters[i].get_height(); row++){
      for(int col = 0; col<letters[i].get_width(); col++){
        if((width*(row-int(ceil(finfo.Descent(str[i])))+max_descent) + current_advance + col + bearing + kern )*2>max_address){
          printf("Programming error in MakeStringImage %d is bigger than %d\n",(width*(row-int(ceil(finfo.Descent(str[i])))+max_descent) + current_advance + col + bearing + kern )*2,max_address);
        } else {
	pixels[(width*(row-int(ceil(finfo.Descent(str[i])))+max_descent) + current_advance + col + bearing + kern )*2]
	  += letter_pixels[letters[i].get_width()*row+col];
	pixels[(width*(row-int(ceil(finfo.Descent(str[i])))+max_descent) + current_advance + col + bearing + kern )*2+1]
	  =  pixels[(width*(row-int(ceil(finfo.Descent(str[i])))+max_descent) + current_advance + col + bearing + kern )*2];
        }
      }
    }
    current_advance += int(ceil(finfo.AdvanceX(str[i])
#ifdef USE_FREETYPE
				+kerning.x/64.0
#endif
				));
    //std::cout << "width of " << str[i] << ": " << finfo.Width(str[i]) << "\n";
    //std::cout << "advance of " << str[i] << ": " << finfo.AdvanceX(str[i]) << "\n";
    //std::cout << "lbearing of " << str[i] << ": " << finfo.LBearing(str[i]) << "\n";
  }
  //std::cout << "fg: " << fg << "\n";
  //std::cout << "bg: " << bg << "\n";
  
  unsigned int fg_r, fg_g, fg_b;
  unsigned int bg_r, bg_g, bg_b;

  std::string fg_red=fg.substr(1,2);
  std::string fg_green=fg.substr(3,2);
  std::string fg_blue=fg.substr(5,2);
  std::string bg_red=bg.substr(1,2);
  std::string bg_green=bg.substr(3,2);
  std::string bg_blue=bg.substr(5,2);
  
  sscanf(fg_red.c_str(),"%x",&fg_r);
  sscanf(fg_green.c_str(),"%x",&fg_g);
  sscanf(fg_blue.c_str(),"%x",&fg_b);
  sscanf(bg_red.c_str(),"%x",&bg_r);
  sscanf(bg_green.c_str(),"%x",&bg_g);
  sscanf(bg_blue.c_str(),"%x",&bg_b);

  unsigned char *rgbapixels = new unsigned char[width*height*IMAGEINFO_RGBA_SIZE];
  /* 
   * This scaling by alpha is what you want to do when target has no alpha channel,
   * but is incorrect if it does have one.
   */
  for(int i=0;i<height;i++){
    for(int j=0; j<width;j++){
       if(scale_alpha){
         rgbapixels[(i*width+j)*IMAGEINFO_RGBA_SIZE]   = (unsigned char)(int(pixels[(i*width+j)*IMAGEINFO_MONOA_SIZE])*(fg_r/256.0)  + int(255-pixels[(i*width+j)*IMAGEINFO_MONOA_SIZE])*(bg_r/256.0));
         rgbapixels[(i*width+j)*IMAGEINFO_RGBA_SIZE+1] = (unsigned char)(int(pixels[(i*width+j)*IMAGEINFO_MONOA_SIZE])*(fg_g/256.0)  + int(255-pixels[(i*width+j)*IMAGEINFO_MONOA_SIZE])*(bg_g/256.0));
         rgbapixels[(i*width+j)*IMAGEINFO_RGBA_SIZE+2] = (unsigned char)(int(pixels[(i*width+j)*IMAGEINFO_MONOA_SIZE])*(fg_b/256.0)  + int(255-pixels[(i*width+j)*IMAGEINFO_MONOA_SIZE])*(bg_b/256.0));
         rgbapixels[(i*width+j)*IMAGEINFO_RGBA_SIZE+3] = 255;
       } else {
         rgbapixels[(i*width+j)*IMAGEINFO_RGBA_SIZE]   = pixels[(i*width+j)*IMAGEINFO_MONOA_SIZE]*(unsigned char)(fg_r/255.0);
         rgbapixels[(i*width+j)*IMAGEINFO_RGBA_SIZE+1] = pixels[(i*width+j)*IMAGEINFO_MONOA_SIZE]*(unsigned char)(fg_g/255.0);
         rgbapixels[(i*width+j)*IMAGEINFO_RGBA_SIZE+2] = pixels[(i*width+j)*IMAGEINFO_MONOA_SIZE]*(unsigned char)(fg_b/255.0);
         rgbapixels[(i*width+j)*IMAGEINFO_RGBA_SIZE+3] = pixels[(i*width+j)*IMAGEINFO_MONOA_SIZE+1];
       }
    }
  }
  delete [] pixels;
  iinfo = image_info(width,height,rgbapixels,IMAGEINFO_RGBA);
  return iinfo;
}
