/*
     pygl/font_info.h: CCP4MG Molecular Graphics Program
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


#ifndef _CCP4MG_FONT_CACHE_
#define _CCP4MG_FONT_CACHE_
#ifdef USE_FREETYPE
#include <ft2build.h>
#include FT_FREETYPE_H
#endif
#include <vector>
#include <string>
#include <iostream>

/* Probably linux specific */
#include <wchar.h>
typedef  wchar_t WIDE_CHAR_T;

enum { FONT_BITMAP_MONO, FONT_BITMAP_GRAYSCALE_FT, FONT_BITMAP_GRAYSCALE };

class MGFontInfo {
 public:
  std::vector<unsigned char *> bm;
  std::string foundry;
  std::string family;
  std::string weight;
  std::string slant;
  std::string encoding;
  std::vector <size_t> widths;
  std::vector <size_t> heights;
  std::vector <float> lbearings;
  std::vector <float> descents;
  std::vector <float> dxs;
  std::vector <float> dys;
  int fn_size;
#ifdef USE_FREETYPE
  FT_Face face;
  FT_Vector Kerning(int left,int right) const;
#endif
  int HaveKerning() const ;
  double BaseLineSkip() const ;
  double UnderLinePosition()  const ;
  double UnderLineWidth()  const ;
  double StrikeThroughPosition()  const ;
  WIDE_CHAR_T base;
  WIDE_CHAR_T last;
  int isValid() const;
  std::string Foundry() const { return foundry; };
  std::string Family() const { return family; };
  std::string Weight() const { return weight; };
  std::string Slant() const { return slant; };
  std::string Encoding() const { return encoding; };
  int Size() const { return fn_size; };
  size_t Width(WIDE_CHAR_T i) const { return (i<=last) ? widths[i-base]:Width('?'); };
  size_t Height(WIDE_CHAR_T i) const { return (i<=last) ? heights[i-base]:Height('?'); };
  float LBearing(WIDE_CHAR_T i) const { return (i<=last) ? lbearings[i-base]:LBearing('?'); };
  float Descent(WIDE_CHAR_T i) const { return (i<=last) ? descents[i-base]:Descent('?'); };
  float AdvanceX(WIDE_CHAR_T i) const { return (i<=last) ? dxs[i-base]:AdvanceX('?'); };
  float AdvanceY(WIDE_CHAR_T i) const { return (i<=last) ? dys[i-base]:AdvanceY('?'); };
  size_t Width(const std::string &s) const { return Width(s[0]); };
  size_t Height(const std::string &s) const { return Height(s[0]); };
  float LBearing(const std::string &s) const { return LBearing(s[0]); };
  float Descent(const std::string &s) const { return Descent(s[0]); };
  float AdvanceX(const std::string &s) const { return AdvanceX(s[0]); };
  float AdvanceY(const std::string &s) const { return AdvanceY(s[0]); };
   
  unsigned char *BitMap(WIDE_CHAR_T i) const;
  unsigned char *FullBitMap(WIDE_CHAR_T i) const;
  int FullBitMapWidth(WIDE_CHAR_T i) const;

  int format;

  /* 
   * The following 3 lines specify an OpenGL texture ID associated with each
   * character. Note that although this implies a usage in an OpenGL context, it
   * does not in any way depend on OpenGL for the general usage of the Font Cache  
   * or for the building of this library.
   */
  std::vector <unsigned int> texture_id;
  void SetTextureID(const std::vector <unsigned int> &texture_id_in) { texture_id = texture_id_in; };
  std::vector <unsigned int> GetTextureID(void) const { return texture_id; };

  MGFontInfo();
  MGFontInfo(const std::vector<std::string> &name, const std::vector<unsigned char *> &bm_in, const std::vector <size_t> &widths_in, const std::vector <size_t> &heights_in, const std::vector <float> &lbearings_in, const std::vector <float> &descents_in, const std::vector <float> &dxs_in, const std::vector <float> &dys_in, WIDE_CHAR_T base_in, WIDE_CHAR_T last_in, int format_in=FONT_BITMAP_MONO
#ifdef USE_FREETYPE
               ,FT_Face face=0
#endif
               );
   ~MGFontInfo();
};

class FontCache {
  static std::vector<MGFontInfo> finfos;
  static std::vector<std::string> glut_fonts;
  static std::vector<void*> glut_font_ids;
public:
  static void* isGlutFont(const std::string &family_in);
  static int NumberOfFonts(); 
  static int isFamilyCached(const std::string &family);
  static int isFaceCached(const std::string &foundry, const std::string &family,  const std::string &weight,  const std::string &slant,  const std::string &encoding);
  static int isCached(const std::string &foundry, const std::string &family,  const std::string &weight,  const std::string &slant,  const std::string &encoding, int size);
  static int LoadFont(const std::string &foundry, const std::string &family,  const std::string &weight,  const std::string &slant,  const std::string &encoding, int size, int truetype=1, int use_uncached_families=1);
  static MGFontInfo GetFont(int i);
  static void SetFont(const MGFontInfo &finfo_in, int i);
  static MGFontInfo GetFont(const std::string &foundry, const std::string &family,  const std::string &weight,  const std::string &slant,  const std::string &encoding, int size, int truetype=1);
  static MGFontInfo GetFont(const std::string &name);
  static int LoadFont(const std::string &name, int truetype=1);
  static int AddFont(const MGFontInfo &finfo);
  static std::vector<std::string> GetFamilyList();
  static std::vector<std::string> GetSortedFamilyList();
  static int LoadAllFonts();
  static int LoadAllGlutFonts();
  static int FindNextFontSize(const std::string &family, const std::string &weight, const std::string &slant, int size, int step);
};

std::vector<std::string> ParseX11Name(const std::string &name);
std::vector<std::string> ParseX11Name(const char *name);

#endif /* _CCP4MG_FONT_CACHE_ */
