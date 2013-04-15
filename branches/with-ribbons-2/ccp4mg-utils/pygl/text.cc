/*
     pygl/text.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009-2010 University of York
     Copyright (C) 2012 STFC

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

#ifdef _WIN32
#include <windows.h>
#endif

#include "cprimitive.h"
#include "help.h"
#ifdef __APPLE_CC__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <vector>

#include <iostream>
#include <sstream>
#include <string>
#include <list>

#include "ppmutil.h"
#include "texture.h"
#include "font_info.h"

#include "rgbreps.h"

#include "mgutil.h"

#ifdef USE_FREETYPE
#include <ft2build.h>
#include FT_FREETYPE_H
#endif

typedef struct StringPos {
   std::string str;
   double x,y;
   std::vector<double> col;
} StringPos;

void BuildTextTextures(MGFontInfo &finfo){
  std::cout << "Building Textures for " << finfo.Family() << "\n"; std::cout.flush();
  glEnable(GL_TEXTURE_2D);
  std::vector <unsigned int> my_texture_id;
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

  for(int i=finfo.base;i<=finfo.last;i++){
    int width  = finfo.Width(i);
    int height = finfo.Height(i);
    unsigned char *pix = finfo.BitMap(i);
    if(!pix)
      std::cout << "Bitmap error!\n"; std::cout.flush();
    image_info tiinfo(width,height,finfo.BitMap(i),IMAGEINFO_MONO);
    if(width<64)
      width = 64;
    if(width>64&&width<128)
      width = 128;
    if(width>128&&width<256)
      width = 256;
    if(width>256)
      width = 256;
    
    if(height<64)
      height= 64;
    if(height>64&&height<128)
      height= 128;
    if(height>128&&height<256)
      height= 256;
    if(height>256)
      height= 256;

    tiinfo.ScaleImage(width,height);

    GLuint texName;
    glGenTextures(1, &texName);
    glBindTexture(GL_TEXTURE_2D, texName);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE_ALPHA, width, height, 0, GL_COLOR_INDEX, GL_UNSIGNED_BYTE, tiinfo.get_pixels());
    my_texture_id.push_back(texName);
  }
  finfo.SetTextureID(my_texture_id);
}

void SimpleText::BitMapFont(int count, int left, const MGFontInfo &finfo, int underline, int strikethrough){

  if(finfo.isValid()){

    if(underline||strikethrough){
      int ul_height = int(ceil(fabs(finfo.UnderLineWidth())));
#if defined(_WIN32) || defined(__WIN32__)
      // This suggests I have done something wrong elsewhere.
      int char_width = int(ceil(finfo.AdvanceX(count)-finfo.LBearing(count)));
#else
      int char_width = int(ceil(finfo.AdvanceX(count)+finfo.LBearing(count)));
#endif
      GLubyte *ulbm = new GLubyte[ul_height*char_width];
      for(int ii=0;ii<ul_height*char_width;ii++)
	ulbm[ii] = 0xff;
      if(finfo.format==FONT_BITMAP_MONO)
        glBitmap (0,0,0,0,0,-baseline+yskip*float(fn_size),0);
      else
        glBitmap (0,0,0,0,0,-baseline+yskip*finfo.BaseLineSkip(),0);
      if(underline){
	glBitmap (0,0,0,0,0, -finfo.UnderLinePosition(),0);
	glDrawPixels(char_width,ul_height,GL_COLOR_INDEX ,GL_UNSIGNED_BYTE,ulbm);
	glBitmap (0,0,0,0,0, finfo.UnderLinePosition(), 0);
      }
      if(strikethrough){
	glBitmap (0,0,0,0,0, -finfo.StrikeThroughPosition(),0);
	glDrawPixels(char_width,ul_height,GL_COLOR_INDEX ,GL_UNSIGNED_BYTE,ulbm);
	glBitmap (0,0,0,0,0, finfo.StrikeThroughPosition(), 0);
      }
      if(finfo.format==FONT_BITMAP_MONO)
        glBitmap (0,0,0,0,0,baseline+yskip*float(fn_size),0);
      else
        glBitmap (0,0,0,0,0,baseline+yskip*finfo.BaseLineSkip(),0);
      delete [] ulbm;
    }


    if(finfo.format==FONT_BITMAP_MONO){
      glDisable(GL_BLEND);
    glBitmap (finfo.Width(count),      // width 
              finfo.Height(count),     // height
             -finfo.LBearing(count),   // xorig
              finfo.Descent(count) - 1 + 
	            yskip*float(fn_size)+baseline,// yorig
              finfo.AdvanceX(count),         // xmove
              finfo.AdvanceY(count),         // ymove
              finfo.BitMap(count));         // bitmap
    } else {

      glEnable(GL_BLEND);
      GLfloat kern=0;
#ifdef USE_FREETYPE
#endif
      
      glBitmap (0,0,0,0,
		finfo.LBearing(count),
		-finfo.Descent(count)-baseline-yskip*finfo.BaseLineSkip(),
		0);

      glDrawPixels(finfo.Width(count), finfo.Height(count), GL_COLOR_INDEX ,GL_UNSIGNED_BYTE, finfo.BitMap(count));

      glBitmap (0,0,0,0,0,
		finfo.Descent(count)+baseline+yskip*finfo.BaseLineSkip(),
		0);

      /* And then advance raster position */
      glBitmap (0,0,0,0,
		finfo.AdvanceX(count)-finfo.LBearing(count)-kern,
		finfo.AdvanceY(count),
		0);
    }

  }

}

void SimpleText::BackspaceBitMapFont(int count, int left, const MGFontInfo &finfo){

  if(finfo.isValid()){

    if(finfo.format==FONT_BITMAP_MONO){
    glBitmap (0,                                         // width 
              0,                                         // height
	     -finfo.LBearing(count),   // xorig
	      finfo.Descent(count) - 1 +
	            yskip*float(fn_size)+baseline,// yorig
	      -finfo.AdvanceX(count),        // xmove
	      -finfo.AdvanceY(count),        // ymove
	      0);                                        // bitmap
    }else{
      GLfloat kern=0;
#ifdef USE_FREETYPE
#endif
          glBitmap (0,                                     // width 
		0,                                     // height
		0,   // xorig
		0,   // yorig
		-finfo.AdvanceX(count)-kern, // xmove
		-finfo.AdvanceY(count),  // ymove
		0);    
    } 
  }

}

void SimpleText::BlankBitMapFont(int count, int left, const MGFontInfo &finfo){

  if(finfo.isValid()){
    if(finfo.format==FONT_BITMAP_MONO){
    glBitmap (0,                                         // width 
              0,                                         // height
	     -finfo.LBearing(count),   // xorig
	      finfo.Descent(count) - 1 +
	            yskip*float(fn_size)+baseline,// yorig
	      finfo.AdvanceX(count),         // xmove
	      finfo.AdvanceY(count),         // ymove
	      0);                                        // bitmap
    }else{
      GLfloat kern=0;
#ifdef USE_FREETYPE
#endif
          glBitmap (0,                                     // width 
		0,                                     // height
		0,   // xorig
		0,   // yorig
		finfo.AdvanceX(count)-kern, // xmove
		finfo.AdvanceY(count),  // ymove
		0);    
    } 
  }

}

int SimpleText::LoadFont(){
  return 0;
  //int id = FontCache::LoadFont("",GetFontFamily(),GetFontWeight(),GetFontSlant(),"",atoi(GetFontSize().c_str()),1,0);
  return id;
}

int SimpleText::GetFontSize() const {
    return fn_size;
  int orig_size = fn_size;
  int mult = SimpleBillBoard::GetMagnification();
  const GLubyte *renderer = glGetString(GL_RENDERER);
  if(mult==1||renderer==0)
    return fn_size;
  int curr_size = orig_size*mult;
  return curr_size;
}

std::string SimpleText::GetFontFamily() const {return family;}
std::string SimpleText::GetFontWeight() const {return weight;}
std::string SimpleText::GetFontSlant() const {return slant;}
std::string SimpleText::GetText(void) const{ return text; }
 
void SimpleText::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_line();
  else
    set_draw_colour_line_override(col);
}

static int text_id = 0; //Ug

Text::Text(const Cartesian &vertex_in, const std::string &text_in, const Cartesian &origin_in, double alpha_in, const std::string &family_in, const int size_in, const std::string &weight_in, const std::string &slant_in) : SimpleText(vertex_in,text_in,origin_in,alpha_in,family_in,size_in,weight_in,slant_in){
}

BillBoardText::BillBoardText(const Cartesian &vertex_in, const std::string &text_in, const Cartesian &origin_in, double alpha_in, const std::string &family_in, const int size_in, const std::string &weight_in, const std::string &slant_in) : SimpleText(vertex_in,text_in,origin_in,alpha_in,family_in,size_in,weight_in,slant_in){
}

SimpleText::SimpleText(const Cartesian &vertex_in, const GLuint tex_id, const GLuint tex_id_b, const int width, const int height, const Cartesian &origin_in, const std::string &text_in, const int size_in,  double alpha_in){
  vertices.push_back(vertex_in);
  slant = "";
  weight = "normal";
  underline = false;
  text = text_in;
  origin = origin_in;
  texture_id = tex_id;
  fn_size = size_in;
  texture_id_b = tex_id_b;
  texture.set_width(width);
  texture.set_height(height);
  id = text_id++; // Need to get rid of this crap
  multicoloured = true;
  text_height=width;
  text_width=height;
  centered = false;
}
Text::Text(const Cartesian &vertex_in, const GLuint tex_id, const GLuint tex_id_b, const int width, const int height, const Cartesian &origin_in, const std::string &text_in, const int size_in, double alpha_in) : SimpleText(vertex_in,tex_id,tex_id_b,width,height,origin_in,text_in,size_in,alpha_in){
}
BillBoardText::BillBoardText(const Cartesian &vertex_in, const GLuint tex_id, const GLuint tex_id_b, const int width, const int height, const Cartesian &origin_in, const std::string &text_in, const int size_in, double alpha_in): SimpleText(vertex_in,tex_id,tex_id_b,width,height,origin_in,text_in,size_in,alpha_in){
}


SimpleText::SimpleText(const Cartesian &vertex_in, const std::string &text_in, const Cartesian &origin_in, double alpha_in, const std::string &family_in, const int size_in, const std::string &weight_in, const std::string &slant_in){
  vertices.push_back(vertex_in);
  text = text_in;
  origin = origin_in;
  alpha = alpha_in;
  family = family_in;
  fn_size = size_in;
  //fn_descent = 0;
  weight = weight_in;
  underline = false;
  slant = slant_in;
  yskip = 0.0;
  id = text_id++; // Need to get rid of this crap
  colour[0] = colour[1] = colour[2] = 0.0;
  text_height=0;
  text_width=0;
  //renderStringToPixmap(); 
  texture_id = 0;
  //initialize();
  multicoloured = false;
  centered = false;
}

void SimpleText::initialize(){

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glGenTextures(1,&texture_id);
  glBindTexture( GL_TEXTURE_2D, texture_id );
  glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,texture.get_width(),texture.get_height(),0,GL_RGBA,GL_UNSIGNED_BYTE,texture.get_pixels());
  glGenTextures(1,&texture_id_b);
  glBindTexture( GL_TEXTURE_2D, texture_id_b );
  glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,texture.get_width(),texture.get_height(),0,GL_RGBA,GL_UNSIGNED_BYTE,texture.get_pixels());

}

int SimpleText::GetID(void) const{
  return id;
}

void SimpleText::SetText(const std::string &text_in){
  text = text_in;
}

void SimpleText::SetFontName(const std::string &family_in, const int size_in, const std::string &weight_in, const std::string &slant_in){
  family = family_in;
  weight = weight_in;
  slant = slant_in;
  fn_size = size_in;
  //renderStringToPixmap(); 
  texture_id = 0;
  //initialize();
}

Text::~Text(){
}

void SimpleText::SetColour(float r, float g, float b, float a) {

  colour[0] = r;
  colour[1] = g;
  colour[2] = b;
  colour[3] = a;
  return;
  glDisable(GL_LIGHTING);
  glColor4f(r,g,b,a);

  GLfloat valsr[256];
  GLfloat valsg[256];
  GLfloat valsb[256];
  GLfloat valsa[256];
  for(int i=0;i<256;i++){
    valsr[i] = float(i)*r/256.0;
    valsg[i] = float(i)*g/256.0;
    valsb[i] = float(i)*b/256.0;
    valsa[i] = float(i)*a/256.0;
    if(valsr[i]>1.0) valsr[i]=1.0;
    if(valsg[i]>1.0) valsg[i]=1.0;
    if(valsb[i]>1.0) valsb[i]=1.0;
  }
  glPixelMapfv(GL_PIXEL_MAP_I_TO_R,256,valsr);
  glPixelMapfv(GL_PIXEL_MAP_I_TO_G,256,valsg);
  glPixelMapfv(GL_PIXEL_MAP_I_TO_B,256,valsb);
  glPixelMapfv(GL_PIXEL_MAP_I_TO_A,256,valsa);
  glPixelTransferi(GL_MAP_COLOR,GL_TRUE);
}

void SimpleText::SetDefaultColour(void) {

  GLfloat params[4];
  glGetFloatv(GL_COLOR_CLEAR_VALUE,params);

  double y = params[0]*0.299 + params[1]*0.587 + params[2]*0.114;

  if(y>0.5)
    SetColour(0.0,0.0,0.0,1.0);
  else
    SetColour(1.0,1.0,1.0,1.0);

}

void Text::SetRasterPosition(double x, double y, double z) const {

  glRasterPos3d(x,y,z);

  GLboolean valid=GL_TRUE;
  glGetBooleanv(GL_CURRENT_RASTER_POSITION_VALID,&valid);
  if(!valid){
    GLdouble rx,ry,rz;
    GLdouble modelMatrix[16];
    GLdouble projMatrix[16];
    GLint viewport[4];
    glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
    glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);
    glGetIntegerv(GL_VIEWPORT,viewport);
    gluUnProject((viewport[2]+viewport[0])/2,(viewport[3]+viewport[1])/2,0.5,modelMatrix,projMatrix,viewport,&rx,&ry,&rz);
    GLdouble winx,winy,winz;
    gluProject(x,y,z,modelMatrix,projMatrix,viewport,&winx,&winy,&winz);
    if(winz>=1.0||winz<=0.0)
      return;
    if(winx>viewport[0]&&winx<viewport[2]&&winy>viewport[1]&&winy<viewport[3]) // clipped by front and back, not sides or top or bottom.
      return;
    glRasterPos3d(rx,ry,rz);
    glBitmap(0, 0, 0, 0, winx-(viewport[2]+viewport[0])/2, winy-(viewport[3]+viewport[1])/2, NULL);
  }
  GLdouble rpos[4];
  glGetDoublev(GL_CURRENT_RASTER_POSITION,rpos);
}

void Text::draw(const double *override_colour, int selective_override){
  double x = 0;//vertices.front().get_x();
  double y = 0;//vertices.front().get_y();
  double z = 0;//vertices.front().get_z();
  glTexCoord2f(0,0);
  glVertex3f(x,y,z);
  glTexCoord2f(1,0);
  glVertex3f(x+2,y,z);
  glTexCoord2f(1,1);
  glVertex3f(x+2,y+2,z);
  glTexCoord2f(0,1);
  glVertex3f(x,y+2,z);
  return;
  SetRasterPosition(x,y,z);
  draw_main();
}

int hex2dec(std::string strval){
  int DEC_VALUE;
  sscanf(strval.c_str(), /*original string*/
	 "%x", /*format string; %x indicates hexadecimal*/
	 &DEC_VALUE); /*address to which the result is written*/
  return DEC_VALUE;
}

void SimpleText::draw_main(const double *override_colour, int selective_override){

	/*
  baseline = 0.0;
  double baselineskip = 0.0;

  int start_of_current_line = -1;

  //const char *s = text.c_str();
  std::string &s = text;

  int i;

  int font_id = LoadFont(); 
  MGFontInfo finfo = FontCache::GetFont(font_id);
  
  baselineskip = finfo.BaseLineSkip();

  std::string current_word = "";
  std::string previous_word = "";
  std::string tag = "";
  std::string tag_value = "";
  std::vector<std::string> tags;


  std::string slant_orig = slant;
  std::string weight_orig = weight;
  std::string fn_size_orig = fn_size;
  std::string family_orig = family;
  double yskip_orig = yskip;

  double raster_colour[4];
  double raster_colour_orig[4];
  glGetDoublev(GL_CURRENT_RASTER_COLOR,raster_colour);
  glGetDoublev(GL_CURRENT_RASTER_COLOR,raster_colour_orig);

  std::vector <std::string> slants;
  std::vector <std::string> weights;
  std::vector <std::string> fn_sizes;
  std::vector <double> yskips;
  std::vector <std::string> families;
  std::vector <std::vector<double> > raster_colours;

  std::vector <int> underlines;
  std::vector <int> strikeouts;
  underlines.push_back(0);
  strikeouts.push_back(0);

  slants.push_back(slant_orig);
  weights.push_back(weight_orig);
  fn_sizes.push_back(fn_size_orig);
  yskips.push_back(yskip_orig);
  families.push_back(family_orig);
  std::vector<double> raster_colour_orig_v;
  raster_colour_orig_v.push_back(raster_colour_orig[0]);
  raster_colour_orig_v.push_back(raster_colour_orig[1]);
  raster_colour_orig_v.push_back(raster_colour_orig[2]);
  raster_colour_orig_v.push_back(raster_colour_orig[3]);
  raster_colours.push_back(raster_colour_orig_v);

  unsigned int state = 0;
  unsigned int IN_TAG = 1;
  unsigned int block_depth = 0;
  int closing = 0;

  int just_done_subsup = 0;
  int doing_subsup = 0;
  int current_position=0;
  int blanking=0;

  
  void *glut_font;
  bool fast_font = false;
  if((glut_font=FontCache::isGlutFont(family)))
    fast_font = true;

  for (i = 0; i < int(text.size()); i++){
    if(i==current_position){
       blanking = 0;
    }
    if(text[i]=='<'){
      if(state&IN_TAG)
         printf("already in tag, syntax error\n");
      state |= IN_TAG;
      if(i==int(text.size()-1))
         printf("unclosed tag at end of string.\n");
      if(text[i+1]=='/'&&block_depth<1) 
         printf("closing block, with out being in block\n");
      if(text[i+1]!='/'){
         block_depth++;
	 closing = 0;
      }else{
	 closing = 1;
      }
      if(text[i+1]=='/'&&block_depth>0) 
         block_depth--;

      slant = slants.back();
      weight = weights.back();
      fn_size = fn_sizes.back();
      yskip = yskips.back();
      family = families.back();

      finfo = FontCache::GetFont(LoadFont());
      if(!getenv("DISABLE_BITMAP_FONTS")&&!fast_font){
        if(current_word.size() > 0 && !doing_subsup){
           just_done_subsup = 0;
        }
        if(doing_subsup&&just_done_subsup){
          for(unsigned int j = 0; j < previous_word.size(); j++){
	    int prev = (j>0)                      ? previous_word[j-1]:-1;
            BackspaceBitMapFont(previous_word[j],prev,finfo);
	  }
        }
      }

      for(unsigned int j = 0; j < current_word.size(); j++){
         if(!getenv("DISABLE_BITMAP_FONTS")&&!fast_font){
	   int prev = (j>0)                     ? current_word[j-1]:-1;
	   if(blanking)
	     BlankBitMapFont(current_word[j],prev,finfo);
	   else
	     BitMapFont(current_word[j],prev,finfo,underlines.back(),strikeouts.back());
	 }else{
	   if(!blanking){
	     //glutBitmapCharacter (glut_font, current_word[j]);
	   }else{
             //glBitmap (0,0,0,0,glutBitmapWidth (glut_font, current_word[j]),0,0);
	   }
         }
      }
 
      if(!getenv("DISABLE_BITMAP_FONTS")&&!fast_font){
        if(doing_subsup&&just_done_subsup){
	  if(previous_word.size()>current_word.size()){
            for(unsigned int j = 0; j < current_word.size(); j++){
	      int prev = (j>0)                      ? current_word[j-1]:-1;
              BackspaceBitMapFont(current_word[j],prev,finfo);
	    }
            for(unsigned int j = 0; j < previous_word.size(); j++){
	      int prev = (j>0)                      ? previous_word[j-1]:-1;
              BlankBitMapFont(previous_word[j],prev,finfo);
	    }
	  }
        }
      }

      if(current_word.size()>0) previous_word = current_word;
      current_word.erase(current_word.begin(),current_word.end());
    }

    if(text[i]=='>'){
      if(!(state&IN_TAG))
         printf("not in tag, syntax error\n");
      state ^= IN_TAG;
      if (tag.size()>0) tags.push_back(std::string(tag));
      if(closing){ 
	 if(0)//strcmp(tags.back().c_str()+1,(tags.end()-2)->c_str()))
           printf("%s does not close %s\n",tags.back().c_str(),(tags.end()-2)->c_str());
	 else{
	   if(tag=="/greek"||tag=="/font_family")
	     if(families.size()>0) families.pop_back();
	   if(tag=="/i"||tag=="/r"||tag=="/o")
	     if(slants.size()>0) slants.pop_back();
	   if(tag=="/b")
	     if(weights.size()>0) weights.pop_back();
	   if(tag=="/font_size"){
	     if(fn_sizes.size()>0) fn_sizes.pop_back();
	   }
	   if(tag=="/u")
	     underlines.pop_back();
	   if(tag=="/strike")
	     strikeouts.pop_back();
	   if(tag=="/sub"||tag=="/sup"){
	     if(fn_sizes.size()>0) fn_sizes.pop_back();
	     if(yskips.size()>0) yskips.pop_back();
	     just_done_subsup = 1;
	     doing_subsup = 0;
	   }
	   if(tag=="/colour"&&i>current_position){
	     current_position = i;
	     blanking = 1;
             i = start_of_current_line;
             double x = vertices.front().get_x();
             double y = vertices.front().get_y();
             double z = vertices.front().get_z();
	     if(raster_colours.size()>0) raster_colours.pop_back();
             std::vector<double> raster_colour_v = raster_colours.back();
	     SetColour(raster_colour_v[0],raster_colour_v[1],raster_colour_v[2],raster_colour_v[3]);
             SetRasterPosition(x,y,z);
	     if(tags.size()>0) tags.pop_back();
	   }
	   if(tags.size()>0) tags.pop_back();
	   if(tags.size()>0) tags.pop_back();
	 }
      }else{
	 int loc =  tag.find_first_of('=');
	 if(loc!=-1){
	   tag_value=tag.substr(loc+2,tag.size()-3-loc);
	   tag = tag.substr(0,loc);
	   tags.pop_back();
	   tags.push_back(tag);
	 }
	 if(tag=="br"){
           double x = vertices.front().get_x();
           double y = vertices.front().get_y();
           double z = vertices.front().get_z();
	   SetRasterPosition(x,y,z);
	   tags.pop_back();
           baseline += baselineskip; //yskip = yskips.back() = yskips.back() + 1;
	   start_of_current_line = i;
	 }
	 if(tag=="colour"&&i>current_position){
	   current_position = i;
	   blanking = 1;
	   if(tags.size()>0) tags.pop_back();
	   while(fn_sizes.size()>1)
	     fn_sizes.pop_back();
	   while(families.size()>1)
	     families.pop_back();
	   while(yskips.size()>1)
	     yskips.pop_back();
	   while(weights.size()>1)
	     weights.pop_back();
	   while(slants.size()>1)
	     slants.pop_back();
	   tag.erase(tag.begin(),tag.end());
           i = start_of_current_line;

           double x = vertices.front().get_x();
           double y = vertices.front().get_y();
           double z = vertices.front().get_z();

	   GLfloat colour[4] = {0.0, 0.0, 0.0, 1.0};
	   if(tag_value[0]=='#'){
	     colour[0] = (float)hex2dec(tag_value.substr(1,2))/255.0;
	     colour[1] = (float)hex2dec(tag_value.substr(3,2))/255.0;
	     colour[2] = (float)hex2dec(tag_value.substr(5,2))/255.0;
	   }else{
	     colour[0] = RGBReps::GetColour(tag_value)[0];
	     colour[1] = RGBReps::GetColour(tag_value)[1];
	     colour[2] = RGBReps::GetColour(tag_value)[2];
	   }
	   SetColour(colour[0],colour[1],colour[2],colour[3]);
           SetRasterPosition(x,y,z);
           glGetDoublev(GL_CURRENT_RASTER_COLOR,raster_colour);
           std::vector<double> raster_colour_v;
           raster_colour_v.push_back(raster_colour_orig[0]);
           raster_colour_v.push_back(raster_colour_orig[1]);
           raster_colour_v.push_back(raster_colour_orig[2]);
           raster_colour_v.push_back(raster_colour_orig[3]);
           raster_colours.push_back(raster_colour_v);
	 }
	 if(tag=="font_size")
	   fn_sizes.push_back(IntToString(FontCache::FindNextFontSize(GetFontFamily(),GetFontWeight(),GetFontSlant(), atoi(GetFontSize().c_str()),atoi(tag_value.c_str()))));
	 if(tag=="font_family")
	   families.push_back(tag_value);
	 if(tag=="greek")
	   families.push_back(std::string("symbol"));
	 if(tag=="sub")
	   yskips.push_back(yskips.back()+0.4);
	 if(tag=="sup")
	   yskips.push_back(yskips.back()-0.4);
	 if(tag=="sub"||tag=="sup"){
	   doing_subsup = 1;

         if(!getenv("DISABLE_BITMAP_FONTS")&&!fast_font)
	   fn_sizes.push_back(IntToString(FontCache::FindNextFontSize(GetFontFamily(),GetFontWeight(),GetFontSlant(), atoi(GetFontSize().c_str()),-2))); // Need something different from -2, a bit more intelligence.
	 }
	 if(tag=="b")
	   weights.push_back("bold");
	 if(tag=="r")
	   slants.push_back("r");
	 if(tag=="i")
	   slants.push_back("i");
	 if(tag=="o")
	   slants.push_back("o");
	 if(tag=="u")
	   underlines.push_back(1);
	 if(tag=="strike")
	   strikeouts.push_back(1);
      }

      tag.erase(tag.begin(),tag.end());
    }

    if(state&IN_TAG&&text[i]!='<'){
      if(1){//!blanking){
        tag.push_back(tolower(s[i]));
      }
    }else if(text[i]!='<'&&text[i]!='>'&&i>start_of_current_line){
      if(text[i]=='&'){
	i++;
	std::string special;
	while(text[i]!=';'&&i<int(text.size())&&text[i]!=' '&&text[i]!='<'){
	  special.push_back(text[i]);
	  i++;
	}
	if(text[i]!=';') i--;
	if(special=="")
	  current_word.push_back('&');
	if(special=="lt")
	  current_word.push_back('<');
	if(special=="gt")
	  current_word.push_back('>');
	if(special=="amp")
	  current_word.push_back('&');
      }else
	current_word.push_back(s[i]);
    }

  }

  slant = slants.back();
  weight = weights.back();
  fn_size = fn_sizes.back();
  yskip = yskips.back();
  family = family_orig;

  finfo = FontCache::GetFont(LoadFont());
  for(unsigned int j = 0; j < current_word.size(); j++){
    int prev = (j>0)                     ? current_word[j-1]:-1;

    if(!getenv("DISABLE_BITMAP_FONTS")&&!fast_font){
      BitMapFont(current_word[j],prev,finfo,underlines.back(),strikeouts.back());
    }
    //else
      //glutBitmapCharacter (glut_font, current_word[j]);
  }

  current_word.erase(current_word.begin(),current_word.end());

  slant = slant_orig;
  weight = weight_orig;
  fn_size = fn_size_orig;
  yskip = yskip_orig;
  family = family_orig;
  */

}  

std::string X112PSFontName(const std::string &name_in, const std::string &weight, const std::string &slant);

void PSBlank(char c){
}

void PSBackspace(char c){
}

int CalculateStringSize(const MGFontInfo &finfo, const std::string &str){
  if(!finfo.isValid()) return 0;
  int have_kerning = 0;
  if(finfo.HaveKerning()){
    have_kerning = 1;
  }
  int kern = 0;

#ifdef USE_FREETYPE
  FT_Vector kerning = {0,0};
#endif

  int current_advance = -int(ceil(finfo.LBearing(str[0])));;
  for(unsigned i=0;i<str.size();i++){
    int bearing = 0;
    bearing = int(ceil(finfo.LBearing(str[i])));
#ifdef USE_FREETYPE
    if(i>0&&have_kerning){
      kerning = finfo.Kerning(str[i-1],str[i]);
      //std::cout << "kerning of " << str[i-1] << ", " << str[i] << ": " << kerning.x/64.0 << "\n";
      kern = int(ceil(kerning.x/64.0));
    }
#endif
    current_advance += int(ceil(finfo.AdvanceX(str[i])
#ifdef USE_FREETYPE
                    +kerning.x/64.0
#endif
                    ));
  }
  return current_advance;
}

void SimpleText::DrawPSmain(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v, bool is_a_bill_board){

  /* Need Some Way of making sure we are using PS fonts. Bit baffled on that one ...
  // Start Fancy stuff
  int font_id = LoadFont(); 
  MGFontInfo finfo = FontCache::GetFont(font_id);

  std::vector<StringPos> str_pos;
  StringPos this_str_pos;

  std::string stripped;
  std::string psfname = X112PSFontName(GetFontFamily(),GetFontWeight(),GetFontSlant());
 
  baseline = 0.0;
  int start_of_current_line = -1;

  //const char *s = text.c_str();
  std::string &s = text;

  int i;

  std::string current_word = "";
  std::string previous_word = "";
  std::string tag = "";
  std::string tag_value = "";
  std::vector<std::string> tags;


  std::string slant_orig = slant;
  std::string weight_orig = weight;
  std::string fn_size_orig = fn_size;
  std::string family_orig = family;
  double yskip_orig = 0.0; // Could be yskip.
  double *colour_orig = colour;
  double *ccolour = colour;

  std::vector <std::string> slants;
  std::vector <std::string> weights;
  std::vector <std::string> fn_sizes;
  std::vector <double> yskips;
  std::vector <std::string> families;
  std::vector <double*> colours;

  std::vector <int> underlines;
  std::vector <int> strikeouts;
  underlines.push_back(0);
  strikeouts.push_back(0);

  slants.push_back(slant_orig);
  weights.push_back(weight_orig);
  fn_sizes.push_back(fn_size_orig);
  yskips.push_back(yskip_orig);
  families.push_back(family_orig);
  colours.push_back(colour_orig);

  unsigned int state = 0;
  unsigned int IN_TAG = 1;
  unsigned int block_depth = 0;
  int closing = 0;

  int just_done_subsup = 0;
  int doing_subsup = 0;
  int current_position=0;
  int blanking=0;
  double current_pos = 0.0;

  for (i = 0; i < int(text.size()); i++){
    if(i==current_position){
       blanking = 0;
    }
    if(text[i]=='<'){
      if(state&IN_TAG)
         printf("already in tag, syntax error\n");
      state |= IN_TAG;
      if(i==int(text.size()-1))
         printf("unclosed tag at end of string.\n");
      if(text[i+1]=='/'&&block_depth<1) 
         printf("closing block, with out being in block\n");
      if(text[i+1]!='/'){
         block_depth++;
	 closing = 0;
      }else{
	 closing = 1;
      }
      if(text[i+1]=='/'&&block_depth>0) 
         block_depth--;

      slant = slants.back();
      weight = weights.back();
      fn_size = fn_sizes.back();
      yskip = yskips.back();
      family = families.back();
      ccolour = colours.back();

      finfo = FontCache::GetFont(LoadFont());

      if(current_word.size() > 0 && !doing_subsup){
         just_done_subsup = 0;
      }
      if(doing_subsup&&just_done_subsup){
	//std::cout << "Need to go back to ...\n";
	//std::cout << "... beginning of " << previous_word << "\n"; 
	current_pos -= CalculateStringSize(finfo,previous_word);
	//std::cout << "current_pos backspaced to " << current_pos << "\n";
      }

      if(!blanking)
        stripped += current_word;

      if(current_word.size()>0){
         //std::cout << "Advance of " << current_word << " is " << CalculateStringSize(finfo,current_word) << "\n";
         if(blanking){
	   //std::cout << "Blanking ---" << current_word << "---\n";
	   current_pos += CalculateStringSize(finfo,current_word);
	   //std::cout << "current_pos blanked to " << current_pos << "\n";
	 }else{
	   this_str_pos.str = current_word;
           this_str_pos.col.clear();
           this_str_pos.col.push_back(ccolour[0]);
           this_str_pos.col.push_back(ccolour[1]);
           this_str_pos.col.push_back(ccolour[2]);
           if(str_pos.size()==0){
             this_str_pos.x = 0.0;
             this_str_pos.y = -yskips.back()*finfo.Size();
             //std::cout << "This is first word\n";
	   }else{
             this_str_pos.x = current_pos;
             this_str_pos.y = -yskips.back()*finfo.Size();
             //std::cout << "This is not first word\n";
	   }
	   current_pos += CalculateStringSize(finfo,current_word);
	   str_pos.push_back(this_str_pos);
	   //std::cout << "current_pos: " << current_pos << "\n";
	   //std::cout << "Writing ---" << current_word << "---\n";
	   //std::cout << "At " << str_pos.back().x << " " << str_pos.back().y << "\n";
	   //std::cout << "ccolour: " << ccolour[0] << " " << ccolour[1] << " " << ccolour[2] << "\n";
         } 
      }
 
      if(doing_subsup&&just_done_subsup){
	if(previous_word.size()>current_word.size()){
	  //std::cout << "Backspacing ---" << current_word << "---\n";
	  current_pos -= CalculateStringSize(finfo,current_word);
	  //std::cout << "current_pos backspaced to " << current_pos << "\n";
	  //std::cout << "Blanking ---" << previous_word << "---\n";
	  current_pos += CalculateStringSize(finfo,previous_word);
	  //std::cout << "current_pos blanked to " << current_pos << "\n";
	}
      }
      //std::cout << "current_word: " << current_word << "\n";
      //std::cout << "previous_word: " << previous_word << "\n";
      if(current_word.size()>0) previous_word = current_word;
      current_word.erase(current_word.begin(),current_word.end());
    }

    if(text[i]=='>'){
      if(!(state&IN_TAG))
         printf("not in tag, syntax error\n");
      state ^= IN_TAG;
      if (tag.size()>0) tags.push_back(std::string(tag));
      if(closing){ 
	 if(0)//strcmp(tags.back().c_str()+1,(tags.end()-2)->c_str()))
           printf("%s does not close %s\n",tags.back().c_str(),(tags.end()-2)->c_str());
	 else{
	   if(tag=="/greek"||tag=="/font_family")
	     if(families.size()>0) families.pop_back();
	   if(tag=="/i"||tag=="/r"||tag=="/o")
	     if(slants.size()>0) slants.pop_back();
	   if(tag=="/b")
	     if(weights.size()>0) weights.pop_back();
	   if(tag=="/font_size"){
	     if(fn_sizes.size()>0) fn_sizes.pop_back();
	   }
	   if(tag=="/u")
	     underlines.pop_back();
	   if(tag=="/strike")
	     strikeouts.pop_back();
	   if(tag=="/sub"||tag=="/sup"){
	     if(fn_sizes.size()>0) fn_sizes.pop_back();
	     if(yskips.size()>0) yskips.pop_back();
	     just_done_subsup = 1;
	     doing_subsup = 0;
	     // Make font size += 2
	   }
	   if(tag=="/colour"){
	     if(colours.size()>0) colours.pop_back();
	   }
	   if(tags.size()>0) tags.pop_back();
	   if(tags.size()>0) tags.pop_back();
	 }
      }else{
	 int loc =  tag.find_first_of('=');
	 if(loc!=-1){
	   tag_value=tag.substr(loc+2,tag.size()-3-loc);
	   tag = tag.substr(0,loc);
	   tags.pop_back();
	   tags.push_back(tag);
	 }
	 if(tag=="br"){
	   tags.pop_back();
           baseline += atof(fn_size_orig.c_str()); //yskip = yskips.back() = yskips.back() + 1;
	   start_of_current_line = i;
	 }
	 if(tag=="colour"){
	   if(tags.size()>0) tags.pop_back();
	   double *tcolour = new double[4]; // MEMORY LEAK
	   if(tag_value[0]=='#'){
	     //std::cout << "RGB: " << (float)hex2dec(tag_value.substr(1,2)) << " ";
	     //std::cout << (float)hex2dec(tag_value.substr(3,2)) << " ";
	     //std::cout << (float)hex2dec(tag_value.substr(5,2)) << "\n";
	     tcolour[0] = (float)hex2dec(tag_value.substr(1,2))/255.0;
	     tcolour[1] = (float)hex2dec(tag_value.substr(3,2))/255.0;
	     tcolour[2] = (float)hex2dec(tag_value.substr(5,2))/255.0;
	   }else{
	     tcolour[0] = RGBReps::GetColour(tag_value)[0];
	     tcolour[1] = RGBReps::GetColour(tag_value)[1];
	     tcolour[2] = RGBReps::GetColour(tag_value)[2];
	   }
           //std::cout << "tcolour: " << tcolour[0] << " " << tcolour[1] << " " << tcolour[2] << "\n";
           colours.push_back(tcolour);
	 }
	 if(tag=="font_size")
	   fn_sizes.push_back(IntToString(atoi(fn_sizes.back().c_str())+(atoi(tag_value.c_str()))));
	 if(tag=="font_family")
	   families.push_back(tag_value);
	 if(tag=="greek")
	   families.push_back(std::string("symbol"));
	 if(tag=="sub")
	   yskips.push_back(yskips.back()+0.4);
	 if(tag=="sup")
	   yskips.push_back(yskips.back()-0.4);
	 if(tag=="sub"||tag=="sup"){
	   doing_subsup = 1;
	   fn_sizes.push_back(IntToString(atoi(fn_sizes.back().c_str())-2));
	 }
	 if(tag=="b")
	   weights.push_back("bold");
	 if(tag=="r")
	   slants.push_back("r");
	 if(tag=="i")
	   slants.push_back("i");
	 if(tag=="o")
	   slants.push_back("o");
	 if(tag=="u")
	   underlines.push_back(1);
	 if(tag=="strike")
	   strikeouts.push_back(1);
      }

      tag.erase(tag.begin(),tag.end());
    }

    if(state&IN_TAG&&text[i]!='<'){
      if(1){//!blanking){
        tag.push_back(tolower(s[i]));
      }
    }else if(text[i]!='<'&&text[i]!='>'&&i>start_of_current_line){
      current_word.push_back(s[i]);
    }

  }

  slant = slants.back();
  weight = weights.back();
  fn_size = fn_sizes.back();
  yskip = yskips.back();
  family = family_orig;
  ccolour = colours.back();

  if(current_word.size()>0){
    StringPos this_str_pos;
    this_str_pos.str = current_word;
    if(str_pos.size()==0){
      this_str_pos.x = 0.0;
      this_str_pos.y = -yskips.back()*finfo.Size();
      //std::cout << "This is first word\n";
     }else{
      this_str_pos.x = current_pos;
      this_str_pos.y = -yskips.back()*finfo.Size();
      //std::cout << "This is not first word\n";
    }
    this_str_pos.col.push_back(ccolour[0]);
    this_str_pos.col.push_back(ccolour[1]);
    this_str_pos.col.push_back(ccolour[2]);
    current_pos += CalculateStringSize(finfo,current_word);
    str_pos.push_back(this_str_pos);
    //std::cout << "Writing ---" << current_word << "---\n";
    //std::cout << "At " << str_pos.back().x << " " << str_pos.back().y << "\n";
    //std::cout << "ccolour: " << ccolour[0] << " " << ccolour[1] << " " << ccolour[2] << "\n";
  }
  stripped += current_word;

  current_word.erase(current_word.begin(),current_word.end());

  slant = slant_orig;
  weight = weight_orig;
  fn_size = fn_size_orig;
  yskip = yskip_orig;
  family = family_orig;
  
  //for(unsigned itpos=0;itpos<str_pos.size();itpos++){
    //std::cout << "Writing ---" << str_pos[itpos].str << "---\n";
    //std::cout << "At " << str_pos[itpos].x << " " << str_pos[itpos].y << "\n";
    //std::cout << "With colour: " << str_pos[itpos].col[0] << " " << str_pos[itpos].col[1] << " " << str_pos[itpos].col[2] << "\n";
  //}
  // End Fancy stuff
  
  std::string psname = X112PSFontName(GetFontFamily(),GetFontWeight(),GetFontSlant());
  fp << "/" << psname << " ff " << GetFontSize() << " scf sf ";
  if(is_a_bill_board){
    p = vertices[0];
    for(unsigned itpos=0;itpos<str_pos.size();itpos++){
      fp << p.get_x()*xoff*2+str_pos[itpos].x << " " << p.get_y()*yoff*2+str_pos[itpos].y << " m ";
      fp << "gs 1 1 sc ("<< str_pos[itpos].str << ") " << str_pos[itpos].col[0] << " " << str_pos[itpos].col[1] << " " << str_pos[itpos].col[2] << " srgb sh gr %TEXT \n";
    }
  }
  if(!is_a_bill_board){
    p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
    Cartesian pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    Cartesian pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " m ";
    fp << "gs 1 1 sc ("<< mytext << ") " << 0 << " " << 0 << " " << 0 << " srgb sh gr %TEXT \n";
  }
  */


  Cartesian p;
  std::string mytext = StripTags();
  std::string psname = X112PSFontName(GetFontFamily(),GetFontWeight(),GetFontSlant());
  fp << "/" << psname << " ff " << GetFontSize() << " scf sf ";
  if(!is_a_bill_board){
    p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
    Cartesian pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    Cartesian pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " m ";
  }else{
    p = vertices[0];
    fp << p.get_x()*xoff*2 << " " << p.get_y()*yoff*2 << " m ";
  }
  fp << "gs 1 1 sc ("<< mytext << ") " << 0 << " " << 0 << " " << 0 << " srgb sh gr %TEXT \n";

}

void Text::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v)
{
  DrawPSmain(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,xoff,yoff,xscale,yscale,xscaleps,v,IsBillBoard());
}

std::string SimpleText::StripTags()
{

  return text;
	/*
  std::string stripped;
  std::string psfname = X112PSFontName(GetFontFamily(),GetFontWeight(),GetFontSlant());
 
  baseline = 0.0;
  int start_of_current_line = -1;

  //const char *s = text.c_str();
  std::string &s = text;

  int i;

  std::string current_word = "";
  std::string previous_word = "";
  std::string tag = "";
  std::string tag_value = "";
  std::vector<std::string> tags;


  std::string slant_orig = slant;
  std::string weight_orig = weight;
  std::string fn_size_orig = fn_size;
  std::string family_orig = family;
  double yskip_orig = yskip;
  double *colour_orig = colour;
  double *ccolour = colour;

  std::vector <std::string> slants;
  std::vector <std::string> weights;
  std::vector <std::string> fn_sizes;
  std::vector <double> yskips;
  std::vector <std::string> families;
  std::vector <double*> colours;

  std::vector <int> underlines;
  std::vector <int> strikeouts;
  underlines.push_back(0);
  strikeouts.push_back(0);

  slants.push_back(slant_orig);
  weights.push_back(weight_orig);
  fn_sizes.push_back(fn_size_orig);
  yskips.push_back(yskip_orig);
  families.push_back(family_orig);
  colours.push_back(colour_orig);

  unsigned int state = 0;
  unsigned int IN_TAG = 1;
  unsigned int block_depth = 0;
  int closing = 0;

  int just_done_subsup = 0;
  int doing_subsup = 0;
  int current_position=0;
  int blanking=0;

  for (i = 0; i < int(text.size()); i++){
    if(i==current_position){
       blanking = 0;
    }
    if(text[i]=='<'){
      if(state&IN_TAG)
         printf("already in tag, syntax error\n");
      state |= IN_TAG;
      if(i==int(text.size()-1))
         printf("unclosed tag at end of string.\n");
      if(text[i+1]=='/'&&block_depth<1) 
         printf("closing block, with out being in block\n");
      if(text[i+1]!='/'){
         block_depth++;
	 closing = 0;
      }else{
	 closing = 1;
      }
      if(text[i+1]=='/'&&block_depth>0) 
         block_depth--;

      slant = slants.back();
      weight = weights.back();
      fn_size = fn_sizes.back();
      yskip = yskips.back();
      family = families.back();
      ccolour = colours.back();

      if(current_word.size() > 0 && !doing_subsup){
         just_done_subsup = 0;
      }
      if(doing_subsup&&just_done_subsup){
        //for(unsigned int j = 0; j < previous_word.size(); j++)
          //PSBackspace(previous_word[j]);
      }

      if(!blanking)
        stripped += current_word;

      for(unsigned int j = 0; j < current_word.size(); j++){
         //if(blanking)
           //PSBlank(current_word[j]);
          //else
           //PSWrite(current_word[j],ccolour);
      }
 
      if(doing_subsup&&just_done_subsup){
	//if(previous_word.size()>current_word.size()){
          //for(unsigned int j = 0; j < current_word.size(); j++)
            //PSBackspace(current_word[j]);
          //for(unsigned int j = 0; j < previous_word.size(); j++)
            //PSBlank(previous_word[j]);
	//}
      }
      if(current_word.size()>0) previous_word = current_word;
      current_word.erase(current_word.begin(),current_word.end());
    }

    if(text[i]=='>'){
      if(!(state&IN_TAG))
         printf("not in tag, syntax error\n");
      state ^= IN_TAG;
      if (tag.size()>0) tags.push_back(std::string(tag));
      if(closing){ 
	 if(0)//strcmp(tags.back().c_str()+1,(tags.end()-2)->c_str()))
           printf("%s does not close %s\n",tags.back().c_str(),(tags.end()-2)->c_str());
	 else{
	   if(tag=="/greek"||tag=="/font_family")
	     if(families.size()>0) families.pop_back();
	   if(tag=="/i"||tag=="/r"||tag=="/o")
	     if(slants.size()>0) slants.pop_back();
	   if(tag=="/b")
	     if(weights.size()>0) weights.pop_back();
	   if(tag=="/font_size"){
	     if(fn_sizes.size()>0) fn_sizes.pop_back();
	   }
	   if(tag=="/u")
	     underlines.pop_back();
	   if(tag=="/strike")
	     strikeouts.pop_back();
	   if(tag=="/sub"||tag=="/sup"){
	     if(fn_sizes.size()>0) fn_sizes.pop_back();
	     if(yskips.size()>0) yskips.pop_back();
	     just_done_subsup = 1;
	     doing_subsup = 0;
	     // Make font size += 2
	   }
	   if(tag=="/colour"&&i>current_position){
	     current_position = i;
	     blanking = 1;
             i = start_of_current_line;
	     if(colours.size()>0) colours.pop_back();
	     if(tags.size()>0) tags.pop_back();
	   }
	   if(tags.size()>0) tags.pop_back();
	   if(tags.size()>0) tags.pop_back();
	 }
      }else{
	 int loc =  tag.find_first_of('=');
	 if(loc!=-1){
	   tag_value=tag.substr(loc+2,tag.size()-3-loc);
	   tag = tag.substr(0,loc);
	   tags.pop_back();
	   tags.push_back(tag);
	 }
	 if(tag=="br"){
	   tags.pop_back();
           baseline += atof(fn_size_orig.c_str()); //yskip = yskips.back() = yskips.back() + 1;
	   start_of_current_line = i;
	 }
	 if(tag=="colour"&&i>current_position){
	   current_position = i;
	   blanking = 1;
	   if(tags.size()>0) tags.pop_back();
	   while(fn_sizes.size()>1)
	     fn_sizes.pop_back();
	   while(families.size()>1)
	     families.pop_back();
	   while(yskips.size()>1)
	     yskips.pop_back();
	   while(weights.size()>1)
	     weights.pop_back();
	   while(slants.size()>1)
	     slants.pop_back();
	   tag.erase(tag.begin(),tag.end());
           i = start_of_current_line;

	   double tcolour[4] = {0.0, 0.0, 0.0, 1.0};
	   if(tag_value[0]=='#'){
	     tcolour[0] = (float)hex2dec(tag_value.substr(1,2))/255.0;
	     tcolour[1] = (float)hex2dec(tag_value.substr(3,2))/255.0;
	     tcolour[2] = (float)hex2dec(tag_value.substr(5,2))/255.0;
	   }else{
	     tcolour[0] = RGBReps::GetColour(tag_value)[0];
	     tcolour[1] = RGBReps::GetColour(tag_value)[1];
	     tcolour[2] = RGBReps::GetColour(tag_value)[2];
	   }
           colours.push_back(tcolour);
	 }
	 if(tag=="font_size")
	   fn_sizes.push_back(IntToString(atoi(fn_sizes.back().c_str())+(atoi(tag_value.c_str()))));
	 if(tag=="font_family")
	   families.push_back(tag_value);
	 if(tag=="greek")
	   families.push_back(std::string("symbol"));
	 if(tag=="sub")
	   yskips.push_back(yskips.back()+0.4);
	 if(tag=="sup")
	   yskips.push_back(yskips.back()-0.4);
	 if(tag=="sub"||tag=="sup"){
	   doing_subsup = 1;
	   fn_sizes.push_back(IntToString(atoi(fn_sizes.back().c_str())-2));
	 }
	 if(tag=="b")
	   weights.push_back("bold");
	 if(tag=="r")
	   slants.push_back("r");
	 if(tag=="i")
	   slants.push_back("i");
	 if(tag=="o")
	   slants.push_back("o");
	 if(tag=="u")
	   underlines.push_back(1);
	 if(tag=="strike")
	   strikeouts.push_back(1);
      }

      tag.erase(tag.begin(),tag.end());
    }

    if(state&IN_TAG&&text[i]!='<'){
      if(1){//!blanking){
        tag.push_back(tolower(s[i]));
      }
    }else if(text[i]!='<'&&text[i]!='>'&&i>start_of_current_line){
      current_word.push_back(s[i]);
    }

  }

  slant = slants.back();
  weight = weights.back();
  fn_size = fn_sizes.back();
  yskip = yskips.back();
  family = family_orig;

  //PSWrite(current_word);
  stripped += current_word;

  for(unsigned int j = 0; j < current_word.size(); j++){
      //PSWrite(current_word[j],ccolour);
  }

  current_word.erase(current_word.begin(),current_word.end());

  slant = slant_orig;
  weight = weight_orig;
  fn_size = fn_size_orig;
  yskip = yskip_orig;
  family = family_orig;


  return stripped;
  */
  
}

int SimpleText::GetFontDescent() const {
	return fn_descent;
}

void SimpleText::SetFontDescent(unsigned int fn_isize){
  fn_descent = fn_isize;
}

void SimpleText::SetFontSize(unsigned int fn_isize){
  fn_size = fn_isize;
  //SetFontDescent(fn_size); // Until we know
}

void SimpleText::SetFontFamily(const std::string &Family){
  family = Family;
}

void SimpleText::SetFontWeight(const std::string &Weight){
  weight = Weight;
}

void SimpleText::SetFontSlant(const std::string &Slant){
  slant = Slant;
}

BillBoardText::~BillBoardText(){
}

void BillBoardText::SetRasterPosition(double x_in, double y_in, double z_in) const {

  glLoadIdentity();
  double *projMatrix = GetProjectionMatrix();

  double l = (projMatrix[12] - 1.0)/projMatrix[0]; 
  double r = (projMatrix[12] + 1.0)/projMatrix[0];
  double t = (projMatrix[13] - 1.0)/projMatrix[5];
  double b = (projMatrix[13] + 1.0)/projMatrix[5];
  double n = 1;
  //double f = 1000;
  if(fabs(projMatrix[10])>1e-7){
    n = (projMatrix[14] + 1)/projMatrix[10];
    //f = (projMatrix[14] - 1)/projMatrix[10];
  }

  double x = (x_in-.5)*(r-l);
  double y = (y_in-.5)*(b-t);
  double z =  -n-3;
  
  glRasterPos3d(x,y,z);

  GLboolean valid=GL_TRUE;
  glGetBooleanv(GL_CURRENT_RASTER_POSITION_VALID,&valid);
  if(!valid){
    GLdouble rx,ry,rz;
    double curr_modelMatrix[16];
    double curr_projMatrix[16];
    GLint curr_viewport[4];
    glGetDoublev(GL_PROJECTION_MATRIX,curr_projMatrix);
    glGetDoublev(GL_MODELVIEW_MATRIX,curr_modelMatrix);
    glGetIntegerv(GL_VIEWPORT,curr_viewport);
    gluUnProject((curr_viewport[2]+curr_viewport[0])/2,(curr_viewport[3]+curr_viewport[1])/2,0.5,curr_modelMatrix,curr_projMatrix,curr_viewport,&rx,&ry,&rz);
    GLdouble winx,winy,winz;
    gluProject(x,y,z,curr_modelMatrix,curr_projMatrix,curr_viewport,&winx,&winy,&winz);
    glRasterPos3d(rx,ry,rz);
    glBitmap(0, 0, 0, 0, winx-(curr_viewport[2]+curr_viewport[0])/2, winy-(curr_viewport[3]+curr_viewport[1])/2, NULL);
  }

}

void BillBoardText::draw(const double *override_colour, int selective_override){

  glPushMatrix();
  double x = vertices.front().get_x();
  double y = vertices.front().get_y();
  double z = vertices.front().get_z();
  SetRasterPosition(x,y,z);
  draw_main();
  glPopMatrix();
}
