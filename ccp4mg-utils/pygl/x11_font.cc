/*
     pygl/x11_font.cc: CCP4MG Molecular Graphics Program
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


#ifdef USE_GLX
#include "font_info.h"
#include "x11_font.h"
#include <mgutil.h>
#include <iostream>
#include <vector>
#include <string.h>

void FreeX11FontPath(char **paths){
  XFreeFontPath(paths);
}

char **GetX11FontPath(int *npaths){
  Display *dpy = XOpenDisplay(0);
  if(!dpy) return 0;
  char **fp = XGetFontPath(dpy,npaths);
  XCloseDisplay(dpy);
  return fp;
}

int LoadAllX11Fonts(){

  int nfonts;

  Display *dpy = XOpenDisplay(0);
  if(!dpy) return -1;

  char **cfonts = XListFonts(dpy, "*--0-0-*iso8859-*", -1, &nfonts);
  if(nfonts==0||!cfonts){
    cfonts = XListFonts(dpy, "*--10-*-*iso8859-*", -1, &nfonts);
  }
  if(nfonts==0||!cfonts){
    cfonts = XListFonts(dpy, "*--11-*-*iso8859-*", -1, &nfonts);
  }
  if(nfonts==0||!cfonts){
    cfonts = XListFonts(dpy, "*--12-*-*iso8859-*", -1, &nfonts);
  }
  if(nfonts==0||!cfonts){
    cfonts = XListFonts(dpy, "*iso8859-*", -1, &nfonts);
  }
 
  MGFontInfo finfo;

  for(int i=0;i<nfonts;i++){
    std::vector<std::string> split_name = ParseX11Name(cfonts[i]);
    if(std::string(cfonts[i]).find(std::string("omegaserif"))==std::string::npos){
      if(FontCache::isFamilyCached(split_name[1])==-1){
        //std::cout << "Loading " << cfonts[i] << "\n"; std::cout.flush();
        finfo = LoadXFont(dpy,cfonts[i]);
        if(finfo.isValid()){
          //std::cout << "Adding " << cfonts[i] << "\n"; std::cout.flush();
          FontCache::AddFont(finfo);
        }
      }
    } 
  }	

  if(cfonts) XFreeFontNames(cfonts);
  XCloseDisplay(dpy);
  return 0;

}

void FillBitmap (Display *dpy, Window win, GC gc,
	       unsigned int width, unsigned int height,
	       int x0, int y0, char c, unsigned char *bitmap)
{
  
  XImage *image;
  unsigned int x, y;
  Pixmap pixmap;

  pixmap = XCreatePixmap (dpy, win, 8*width, height, 1);
  XSetForeground(dpy, gc, 0);
  XFillRectangle (dpy, pixmap, gc, 0, 0, 8*width, height);
  XSetForeground(dpy, gc, 1);
  XDrawString (dpy, pixmap, gc, x0, y0, &c, 1); // No unicode ??

  image = XGetImage (dpy, pixmap, 0, 0, 8*width, height, 1, XYPixmap);

  /* Fill the bitmap (X11 and OpenGL are upside down wrt each other).  */
  for (y = 0; y < height; y++)
    for (x = 0; x < 8*width; x++)
      if (XGetPixel (image, x, y))
	bitmap[width*(height - y - 1) + x/8] |= (1 << (7 - (x % 8)));
  
  XFreePixmap (dpy, pixmap);
  XDestroyImage (image);

}

MGFontInfo LoadXFont(Display *dpy, const char *name_in){


  //std::cout << "LoadXFont: " << name_in << "\n";
  if(!dpy) return MGFontInfo();
  MGFontInfo finfo;

  std::string name = std::string(name_in);

  GC gc;
  XGCValues values;
  unsigned long valuemask;

  XFontStruct *fs;

  fs = XLoadQueryFont(dpy, name.c_str());
  if(!fs) {
    //std::cout << "XLoadQueryFont failed\n";
    std::vector<std::string> split_name = ParseX11Name(name);
    if(split_name[1]==""||split_name[2]==""||split_name[3]==""||split_name[4]=="")
      return finfo;
    if(split_name[0]=="")
      split_name[0]="*";
    size_t t = split_name[4].find('*');
    if(t!=std::string::npos)
      split_name[4] = split_name[4].substr(0,t);
    std::string newname = std::string("-") + split_name[0] + std::string("-") + split_name[1] + std::string("-") + split_name[2] + std::string("-") + split_name[3] + std::string("-*--") + split_name[4] + std::string("-*-*-*-*-*-iso8859-*");
    //std::cout << "newname: " << newname << "\n";
    fs = XLoadQueryFont(dpy, newname.c_str());
    name = newname;
    if(!fs) 
      return finfo;
  }
  unsigned int first = fs->min_char_or_byte2;
  unsigned int last = fs->max_char_or_byte2;

  std::vector <unsigned char *> bm;
  std::vector <size_t> widths;
  std::vector <size_t> heights;
  std::vector <float> lbearings;
  std::vector <float> descents;
  std::vector <float> dxs;
  std::vector <float> dys;

  std::string size_str = ParseX11Name(name)[4];

  if(size_str=="0"){
    for(unsigned int c=first; c<last+1; c++){
      unsigned char *bmtmp = new unsigned char[0];
      bm.push_back(bmtmp);
      widths.push_back(0);
      heights.push_back(0);
      lbearings.push_back(0.0);
      descents.push_back(0.0);
      dxs.push_back(0.0);
      dys.push_back(0.0); 
    }
    //std::cout << "Returning empty stuff\n"; std::cout.flush();
    finfo = MGFontInfo(ParseX11Name(name),bm,widths,heights,lbearings,descents,dxs,dys,first,last);
    XFreeFontInfo( 0, fs, 0 );
    return finfo;
  }

  unsigned int max_width, max_height, max_bm_width, max_bm_height;

  /* Allocate a bitmap that can fit all characters.  */
  max_width = fs->max_bounds.rbearing - fs->min_bounds.lbearing;
  max_height = fs->max_bounds.ascent + fs->max_bounds.descent;
  max_bm_width = (max_width + 7) / 8;
  max_bm_height = max_height;


  values.foreground = BlackPixel (dpy, DefaultScreen (dpy));
  values.background = WhitePixel (dpy, DefaultScreen (dpy));
  values.font = fs->fid; // THIS IS IMPORTANT BIT, it decides what XDrawString does !!!
  valuemask = GCForeground | GCBackground | GCFont;
  gc = XCreateGC (dpy, XCreatePixmap (dpy, RootWindow(dpy,DefaultScreen(dpy)), 10, 10, 1), valuemask, &values);

  unsigned int bm_width, bm_height;
  size_t width, height;

  float dx, dy;
  XCharStruct *ch;
  int x, y;
  
  for(unsigned int c=first; c<last+1; c++){
    
    if (fs->per_char && (c >= fs->min_char_or_byte2) && (c <= fs->max_char_or_byte2))
      ch = &fs->per_char[c-fs->min_char_or_byte2];
    else
      ch = &fs->max_bounds;
    
    width = ch->rbearing - ch->lbearing;
    height = ch->ascent + ch->descent;
    dx = ch->width;
    dy = 0;
    
    /* X11's starting point.  */
    x = - ch->lbearing;
    y = ch->ascent;
    
    /* Round the width to a multiple of eight.  We will use this also
       for the pixmap for capturing the X11 font.  This is slightly
       inefficient, but it makes the OpenGL part real easy.  */
    bm_width = (width + 7) / 8;
    bm_height = height;
    
    bm.push_back(new  unsigned char[max_bm_width * max_bm_height]);
    widths.push_back(width);
    heights.push_back(height);
    lbearings.push_back(ch->lbearing);
    descents.push_back(ch->descent);
    dxs.push_back(dx);
    dys.push_back(dy);

    if ((c >= fs->min_char_or_byte2) && (c <= fs->max_char_or_byte2)
	&& (bm_width > 0) && (bm_height > 0))
      {
	memset(bm.back(), '\0', bm_width * bm_height);
	FillBitmap (dpy, RootWindow(dpy,DefaultScreen(dpy)),
		    gc, bm_width, bm_height, x, y, c, bm.back());
      }
  }

  finfo = MGFontInfo(ParseX11Name(name),bm,widths,heights,lbearings,descents,dxs,dys,first,last);

  XFreeFontInfo( 0, fs, 0 );
  XFreeGC (dpy, gc);

  //std::cout << "XLoadFont gives: " << finfo.fn_size << " " << finfo.family << " " << finfo.weight << " " << finfo.slant << " " <<  finfo.bm.size() << " " <<  finfo.widths.size() << " " <<  finfo.heights.size() << " " <<  finfo.lbearings.size() << " " <<  finfo.descents.size() << " " <<  finfo.dxs.size() << " " <<  finfo.dys.size() << "\n";

  return finfo;
  
}

MGFontInfo LoadXFont(Display *dpy, const std::string &foundry, const std::string &family,  const std::string &weight,  const std::string &slant,  const std::string &encoding, int size){
  std::string name = std::string("-") + foundry + std::string("-") + family + std::string("-") + weight + std::string("-") + slant + std::string("-*--") + IntToString(size) + std::string("*-") + encoding;
  return LoadXFont(dpy,name.c_str());
}

MGFontInfo LoadXFont(const std::string &foundry, const std::string &family,  const std::string &weight,  const std::string &slant,  const std::string &encoding, int size){
  std::string name = std::string("-") + foundry + std::string("-") + family + std::string("-") + weight + std::string("-") + slant + std::string("-*--") + IntToString(size) + std::string("*-") + encoding;
  return LoadXFont(name.c_str());
}

MGFontInfo LoadXFont(const char *name){
  Display *dpy = XOpenDisplay(0);
  if(!dpy) return MGFontInfo();
  MGFontInfo finfo = LoadXFont(dpy,name);
  XCloseDisplay(dpy);
  return finfo;
}

int FindNextX11FontSize(const MGFontInfo &finfo, int step){
  if(!finfo.isValid()) return 0;

  Display *dpy = XOpenDisplay(0);
  if(!dpy) return 0;
  int nfonts;
  std::string name;
  int this_size;
  int isize;
  int i;
  int nstep = 0;

  if((finfo.fn_size+step)>0)
    return finfo.Size();

  this_size = finfo.Size();
  isize = this_size;

  std::string newname = std::string("-") + finfo.foundry + std::string("-") + finfo.family + std::string("-") + finfo.weight + std::string("-") + finfo.slant + std::string("-normal--") + std::string("0-*-*-*-*-*-iso8859-*");
  XFontStruct *fs = XLoadQueryFont(dpy, newname.c_str());
  if(fs){
    XFreeFontInfo( 0, fs, 0 );
    return isize+step;
  } 

  if(step<0){
    for(i=isize-1;i>0;i--){
      this_size = i;
      name = std::string("*") + finfo.family + std::string("*") + finfo.weight + std::string("*-") + finfo.slant + std::string("-*--") + IntToString(this_size) + std::string("-*");
      XListFonts(dpy, name.c_str(), 100, &nfonts);
      if(nfonts>0)
        nstep--;
      if(nstep==step){
        XCloseDisplay(dpy);
        return this_size;
      }
    }
  }
  if(step>0){
    for(i=isize+1;i<100;i++){
      this_size = i;
      name = std::string("*") + finfo.family + std::string("*") + finfo.weight + std::string("*-") + finfo.slant + std::string("-*--") + IntToString(this_size) + std::
string("-*");

      XListFonts(dpy, name.c_str(), 100, &nfonts);
      if(nfonts>0)
        nstep++;
      if(nstep==step){
        XCloseDisplay(dpy);
        return this_size;
      }
    }
  }

  XCloseDisplay(dpy);
  return finfo.Size();

}
#endif
