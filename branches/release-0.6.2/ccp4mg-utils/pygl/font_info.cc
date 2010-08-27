/*
     pygl/font_info.cc: CCP4MG Molecular Graphics Program
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


#include "font_info.h"
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <map>
#include <strings.h>
#ifdef USE_GLX
#include "x11_font.h"
#endif
#if defined(_WIN32) || defined(__WIN32__)
#include "win_font.h"
#endif
#include "mgutil.h"
#ifdef __APPLE_CC__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#ifdef _USE_GLUT_
#include <GL/glut.h>
#endif


#ifdef USE_FREETYPE
#include "freetype_font.h"
  extern FT_Error Freetype_Set_Pixel_Sizes( FT_Face face, int fn_size1, int fn_size2 );
  extern FT_Error Freetype_Get_Kerning( FT_Face face, int left, int right, FT_Vector *kerning);
#endif

std::vector<MGFontInfo> FontCache::finfos;
std::vector<std::string> FontCache::glut_fonts;
std::vector<void*> FontCache::glut_font_ids;

int FontCache::FindNextFontSize(const std::string &family, const std::string &weight, const std::string &slant, int size, int step){
  MGFontInfo finfo = GetFont("",family,weight,slant,"",size);
  if(!finfo.isValid()) return 0;
#ifdef USE_FREETYPE
  return FindNextFreeTypeFontSize(finfo, step);
#endif
#ifdef USE_GLX
  return FindNextX11FontSize(finfo, step);
#endif
#if defined(_WIN32) || defined(__WIN32__)
  return FindNextWinFontSize(finfo, step);
#endif
  return size;
}

int MGFontInfo::isValid() const {
  if(FontCache::isGlutFont(family))
    return 1;
  if(family==""||weight==""||slant==""||bm.size()==0||widths.size()==0||heights.size()==0||lbearings.size()==0||descents.size()==0||dxs.size()==0||dys.size()==0)
    return 0;
  else
    return 1;
}

double MGFontInfo::BaseLineSkip() const {
  if(FontCache::isGlutFont(family))
    return 15;
#ifdef USE_FREETYPE
  if(face) {
    Freetype_Set_Pixel_Sizes( face, fn_size, fn_size );
    return face->size->metrics.height/64.0;
  }
#endif
  return fn_size;
}

double MGFontInfo::UnderLinePosition() const {
  if(FontCache::isGlutFont(family))
    return 12;
#ifdef USE_FREETYPE
  if(face) {
    Freetype_Set_Pixel_Sizes( face, fn_size, fn_size );
    if(face->underline_position!=0)
      return -face->underline_position/64.0;
  }
#endif
  if(fn_size<10) return 1;
  return (fn_size+10)/10 - 1;
}

double MGFontInfo::UnderLineWidth() const {
  if(FontCache::isGlutFont(family))
    return 2;
#ifdef USE_FREETYPE
  if(face) {
    Freetype_Set_Pixel_Sizes( face, fn_size, fn_size );
    if(face->underline_thickness!=0)
      return face->underline_thickness/64.0;
  }
#endif
  if(fn_size<10) return 1;
  return (fn_size+10)/10 - 1;
}

double MGFontInfo::StrikeThroughPosition() const {
  if(FontCache::isGlutFont(family))
    return 7;
#ifdef USE_FREETYPE
  if(face) {
    Freetype_Set_Pixel_Sizes( face, fn_size, fn_size );
    return face->size->metrics.height/128.0;
  }
#endif
  return fn_size/2;
}

int MGFontInfo::HaveKerning() const {
  if(FontCache::isGlutFont(family))
    return 0;
#ifdef USE_FREETYPE
  if(!face)
    return 0;
  Freetype_Set_Pixel_Sizes( face, fn_size, fn_size );
  return FT_HAS_KERNING(face);
#else
  return 0;
#endif  
}

#ifdef USE_FREETYPE
FT_Vector MGFontInfo::Kerning(int left,int right) const {
  FT_Vector kerning;
  Freetype_Set_Pixel_Sizes( face, fn_size, fn_size );
  Freetype_Get_Kerning(face,left,right,&kerning);
  return kerning;
}
#endif

int MGFontInfo::FullBitMapWidth(WIDE_CHAR_T i) const { 
  if(FontCache::isGlutFont(family))
    return 9;
  if(format==FONT_BITMAP_MONO)
    return (i<=last) ? ((widths[i-base] + 7) / 8)*8:FullBitMapWidth(bm.size()-1);
  else 
    return (i<=last) ? widths[i-base]:FullBitMapWidth(bm.size()-1);
  return 0;
}
unsigned char *MGFontInfo::FullBitMap(WIDE_CHAR_T i) const { 
  if(FontCache::isGlutFont(family))
    return 0;

  if(i-base>=int(bm.size()))
    return FullBitMap(bm.size()-1);

  if(format==FONT_BITMAP_MONO){
    unsigned val,p;
    int bm_width = (widths[i-base] + 7) / 8;
    unsigned char *bitmap = new unsigned char[8*bm_width*heights[i-base]];
    for(unsigned ii=0;ii<heights[i-base];ii++){
      for(int jj=0;jj<bm_width;jj++){
	int index = (ii*bm_width+jj)*8;
	val = bm[i-base][ii*bm_width+jj];
	p = ((val & 128)/128)*255;
	bitmap[index]   = p;
	p = ((val & 64)/64)*255;
	bitmap[index+1] = p;
	p = ((val & 32)/32)*255;
	bitmap[index+2] = p;
	p = ((val & 16)/16)*255;
	bitmap[index+3]  = p;
	p = ((val & 8)/8)*255;
	bitmap[index+4] = p;
	p = ((val & 4)/4)*255;
	bitmap[index+5] = p;
	p = ((val & 2)/2)*255;
	bitmap[index+6] = p;
	p = ((val & 1))*255;
	bitmap[index+7] = p;
      }
    }
    return bitmap;
  }else if(format==FONT_BITMAP_GRAYSCALE_FT){
    return bm[i-base];
  }else{
    return bm[i-base];
  }
  return 0;
}

void* FontCache::isGlutFont(const std::string &family_in){
  for(unsigned i=0;i<glut_fonts.size();i++){
    if(family_in==glut_fonts[i])
       return glut_font_ids[i];
  }
  return 0;
}

int FontCache::LoadAllGlutFonts(){
/*
#ifdef USE_FREETYPE
  FT_Face face_in;
#endif
  WIDE_CHAR_T base_in = 32;
  WIDE_CHAR_T last_in = 255;
  std::vector<unsigned char *> bm_in;
  std::vector <size_t> widths_in;
  std::vector <size_t> heights_in;
  std::vector <float> lbearings_in;
  std::vector <float> descents_in;
  std::vector <float> dxs_in;
  std::vector <float> dys_in;
  std::vector<std::string> name;

  glut_fonts.clear();
  glut_fonts.push_back(std::string("Fast 8x13"));
  glut_font_ids.push_back(GLUT_BITMAP_8_BY_13);
  glut_fonts.push_back(std::string("Fast 9x15"));
  glut_font_ids.push_back(GLUT_BITMAP_9_BY_15);
  glut_fonts.push_back(std::string("Fast Times Roman 10"));
  glut_font_ids.push_back(GLUT_BITMAP_TIMES_ROMAN_10);
  glut_fonts.push_back(std::string("Fast Times Roman 24"));
  glut_font_ids.push_back(GLUT_BITMAP_TIMES_ROMAN_24);
  glut_fonts.push_back(std::string("Fast Helvetica 10"));
  glut_font_ids.push_back(GLUT_BITMAP_HELVETICA_10);
  glut_fonts.push_back(std::string("Fast Helvetica 12"));
  glut_font_ids.push_back(GLUT_BITMAP_HELVETICA_12);
  glut_fonts.push_back(std::string("Fast Helvetica 18"));
  glut_font_ids.push_back(GLUT_BITMAP_HELVETICA_18);

  name.push_back(std::string("default-foundry"));
  name.push_back(glut_fonts[0]);
  name.push_back(std::string("medium"));
  name.push_back(std::string("r"));
  name.push_back(std::string("15"));
  name.push_back(std::string("iso8859"));
  MGFontInfo f(name,bm_in,widths_in,heights_in,lbearings_in,descents_in,dxs_in,dys_in,base_in,last_in,FONT_BITMAP_MONO
#ifdef USE_FREETYPE
  ,face_in
#endif
  );
  AddFont(f);
  for(unsigned i=1;i<glut_fonts.size();i++){
    name[1] = glut_fonts[i];
    MGFontInfo f1(name,bm_in,widths_in,heights_in,lbearings_in,descents_in,dxs_in,dys_in,base_in,last_in,FONT_BITMAP_MONO
#ifdef USE_FREETYPE
    ,face_in
#endif
    );
    AddFont(f1);
  }
*/
  return 0;
}

int FontCache::LoadAllFonts(){
  int error=0;
#ifdef USE_FREETYPE
  error = LoadAllFreeTypeFonts();
#endif
#if defined(USE_GLX)
  error = LoadAllX11Fonts();
#endif
#if defined(_WIN32) || defined(__WIN32__)
  error = LoadAllWinFonts();
#endif
  return error;
}

int FontCache::NumberOfFonts(){ 
  return finfos.size();
};

int FontCache::LoadFont(const std::string &foundry_in, const std::string &family,  const std::string &weight_in,  const std::string &slant_in,  const std::string &encoding_in, int size, int truetype, int use_uncached_families){
  int error=0;


  //std::cout << "foundry_in: --" << foundry_in << "--\n";

  std::string foundry = foundry_in;
  std::string encoding = encoding_in;
  std::string weight = weight_in;
  std::string slant = slant_in;

  //printf("Request: %s %s %s %s %s %d\n",foundry.c_str(),family.c_str(),weight.c_str(),slant.c_str(),encoding.c_str(),size);

  MGFontInfo finfo;

  if(FontCache::isGlutFont(family)){
    if((error=FontCache::isFaceCached(foundry,family,"medium","r",encoding))>-1){
      return error;
    }
  }
  if((error=FontCache::isCached(foundry,family,weight,slant,encoding,size))>-1){
    return error;
  }

  if(slant_in=="i"){
    //std::cout << "See if \"i\" cached\n";
    if((error=FontCache::isCached(foundry,family,weight,"i",encoding,size))>-1)
      return error;
    if((error=FontCache::isFaceCached(foundry,family,weight,"i",encoding))>-1)
      slant = "i";
    //else
      //std::cout << "\"i\" is not cached\n"; std::cout.flush();
  }

  if(slant_in=="o"){
    //std::cout << "See if \"o\" cached\n";
    if((error=FontCache::isCached(foundry,family,weight,"o",encoding,size))>-1)
      return error;
    if((error=FontCache::isFaceCached(foundry,family,weight,"o",encoding))>-1)
      slant = "o";
    //else
      //std::cout << "\"o\" is not cached\n"; std::cout.flush();
  }

  //std::cout << "slant is " << slant << "\n";

  if(use_uncached_families==0){
    if((error=FontCache::isFamilyCached(family))==-1){
       if(NumberOfFonts()>0)
         return 0;
       else
         return -1;
    }else{
      //std::cout << "Family is cached\n"; std::cout.flush();
      //if(FontCache::isFaceCached(foundry,family,weight,slant,encoding)==-1)
        //return error;
    }
  }

  std::vector<std::string> weights;
  std::vector<std::string> slants;
  std::vector<std::string>::iterator weight_iter;
  std::vector<std::string>::iterator slant_iter;
  if(weight.length()==0){
    weights.push_back("medium");
    weights.push_back("bold");
    weights.push_back("regular");
    weight="normal";
  }else{
    weights.push_back(weight);
    weights.push_back("medium");
    weights.push_back("bold");
    weights.push_back("regular");
  }
  if(slant.length()==0){
    slants.push_back("r");
    slants.push_back("i");
    slants.push_back("o");
    slant="r";
  }else{
    slants.push_back(slant);
    if(slant=="o"){
      slants.push_back("i");
      slants.push_back("r");
    }
    if(slant=="i"){
      slants.push_back("o");
      slants.push_back("r");
    }
    if(slant=="r"){
      slants.push_back("r");
      slants.push_back("i");
      slants.push_back("o");
    }
  }

#ifdef USE_FREETYPE
  //std::cout << "Attempt freetype load\n"; std::cout.flush();
  if((error=FontCache::isFamilyCached(family))>-1){
    //std::cout << "Family is cached\n"; std::cout.flush();
    finfo = GetFont(error);
    if(finfo.format==FONT_BITMAP_GRAYSCALE_FT)
      truetype = 1;
    else
      truetype = 0;
    //if(truetype)
      //std::cout << "and is a freetype\n";
  }
  //MGFontInfo ftmp = GetFont(error);
  //std::cout << "Cached face family: " << ftmp.family << "\n";
  //std::cout << "Cached face weight: " << ftmp.weight << "\n";
  //std::cout << "Cached face slant: " << ftmp.slant << "\n";
  //std::cout << "Cached face size: " << ftmp.Size() << "\n";
  if(truetype){
  //std::cout << "looking for freetype font\n";
  if(foundry_in=="")
    foundry="default-foundry";
  if(encoding_in=="")
    encoding="iso8859";
  int loaded = -1;
  //std::cout << weights.size() << " " << slants.size() << "\n";
  int done = 0;
  weight_iter=weights.begin();
  while(weight_iter!=weights.end()&&!done){
    slant_iter=slants.begin();
    while(slant_iter!=slants.end()&&!done){
      //std::cout << "Trying: " << family << " " << *weight_iter << " " << *slant_iter << "\n"; std::cout.flush();
      int dummy;
      if((dummy=FontCache::isCached(foundry,family,*weight_iter,*slant_iter,encoding,size))==-1){
        //std::cout << "calling LoadFreeTypefont\n"; std::cout.flush();
	finfo = LoadFreeTypeFont(foundry,family,*weight_iter,*slant_iter,encoding,size);std::cout.flush();
        //std::cout << "called LoadFreeTypefont\n"; std::cout.flush();
	if(finfo.foundry.size()>0||finfo.family.size()>0||finfo.widths.size()>0||finfo.heights.size()>0||finfo.lbearings.size()>0||finfo.descents.size()>0||finfo.dxs.size()>0||finfo.dys.size()>0){
	  loaded = AddFont(finfo);
          //printf("%s %s %s %s %s %d being added\n",finfo.foundry.c_str(),finfo.family.c_str(),finfo.weight.c_str(),finfo.slant.c_str(),finfo.encoding.c_str(),finfo.fn_size);
	  //std::cout << "Loaded " << loaded << "\n";
          if((*weight_iter==weight_in||(*weight_iter=="medium"&&weight_in==""))
           &&(*slant_iter==slant_in)||(*slant_iter=="r"&&slant_in=="")||!use_uncached_families)
            done = 1;
	}else{
	  //std::cout << "Didn't Load " << "\n";
        }
      }else{
        //std::cout << "Hmm, thinks it's cached now: " << dummy << "\n";
	return dummy;
      }
      slant_iter++;
    }
    weight_iter++;
  }
  //std::cout << "Finished loop\n"; std::cout.flush();
  if(loaded>-1){
    //std::cout << "returning loaded " << loaded << "\n"; std::cout.flush();
    MGFontInfo tmpfnt = GetFont(loaded);
    std::string slant_req = slant_in;
    std::string weight_req = weight_in;
    if(slant_req=="") slant_req = "r";
    if(weight_req=="") weight_req = "medium";
    //std::cout << "Requested " << slant_req << ", got " << tmpfnt.Slant() << "\n";
    //std::cout << "Requested " << weight_req << ", got " << tmpfnt.Weight() << "\n";
    return loaded;
  }
  }
#endif

#ifdef USE_GLX
  //std::cout << "USE_GLX foundry_in: --" << foundry_in << "--\n";
  if(foundry_in==""||foundry_in=="default-foundry")
    foundry="*";
  if(encoding_in=="")
    encoding="*";
  //std::cout << "USE_GLX foundry: --" << foundry << "--\n";
  weight_iter=weights.begin();
  while(weight_iter!=weights.end()){
    slant_iter=slants.begin();
    while(slant_iter!=slants.end()){
      //std::cout << "GLX Trying: " << family << " " << *weight_iter << " " << *slant_iter << "\n"; std::cout.flush();
      int dummy;
      if((dummy=FontCache::isCached(foundry,family,*weight_iter,*slant_iter,encoding,size))==-1){
	finfo = LoadXFont(foundry,family,*weight_iter,*slant_iter,encoding,size);
	if(finfo.foundry.size()>0||finfo.family.size()>0||finfo.widths.size()>0||finfo.heights.size()>0||finfo.lbearings.size()>0||finfo.descents.size()>0||finfo.dxs.size()>0||finfo.dys.size()>0)
	  return AddFont(finfo);
      }else{
	return dummy;
      }
      slant_iter++;
    }
    weight_iter++;
  }
#endif

#if defined(_WIN32) || defined(__WIN32__)
  if(foundry_in=="")
    foundry="*";
  if(encoding_in=="")
    encoding="iso8859";
  weight_iter=weights.begin();
  while(weight_iter!=weights.end()){
    slant_iter=slants.begin();
    while(slant_iter!=slants.end()){
      finfo = LoadWinFont(foundry,family,*weight_iter,*slant_iter,encoding,size);
      if(finfo.foundry.size()>0||finfo.family.size()>0||finfo.widths.size()>0||finfo.heights.size()>0||finfo.lbearings.size()>0||finfo.descents.size()>0||finfo.dxs.size()>0||finfo.dys.size()>0)
	return AddFont(finfo);
      slant_iter++;
    }
    weight_iter++;
  }
#endif
  return error;
}

int FontCache::LoadFont(const std::string &name, int truetype){
  std::vector<std::string> name_split = ParseX11Name(name);
  if(name_split.size()<6)
    return 0;
  return LoadFont(name_split[0],name_split[1],name_split[2],name_split[3],name_split[4],atoi(name_split[5].c_str()),truetype);
}

std::vector<std::string> FontCache::GetSortedFamilyList(){
  std::vector<std::string> families = GetFamilyList();
  std::vector<std::string> non_glut_families;
  std::vector<std::string> glut_families;
  std::vector<std::string>::iterator fam_iter=families.begin();
  while(fam_iter!=families.end()){
    if(FontCache::isGlutFont(*fam_iter))
       glut_families.push_back(*fam_iter);
    else
       non_glut_families.push_back(*fam_iter);
    fam_iter++;
  }
  std::sort(glut_families.begin(),glut_families.end());
  std::sort(non_glut_families.begin(),non_glut_families.end());
  families.clear();
  families.insert(families.end(),glut_families.begin(),glut_families.end());
  families.insert(families.end(),non_glut_families.begin(),non_glut_families.end());
  return families;
}

std::vector<std::string> FontCache::GetFamilyList(){

  std::vector<std::string> families;
  std::map<std::string,int> families_map;

  std::vector<MGFontInfo>::const_iterator finfo_iter=finfos.begin();
  while(finfo_iter!=finfos.end()){
    families_map[finfo_iter->family] = 1;
    finfo_iter++;
  }
  std::map<std::string,int>::iterator map_iter=families_map.begin();
  while(map_iter!=families_map.end()){
    families.push_back(map_iter->first);
    map_iter++;
  }

  return families;

}

int FontCache::AddFont(const MGFontInfo &finfo){
  finfos.push_back(finfo);
  return finfos.size()-1;
}

int FontCache::isFamilyCached(const std::string &family){
  std::vector<std::string> fams = GetFamilyList();
  std::vector<std::string>::iterator fiter=fams.begin();
  int i=0;
  while(fiter!=fams.end()){
    if(strcasecmp(fiter->c_str(),family.c_str())==0)
      return i;
    fiter++; i++;
  }
  return -1;
}

void FontCache::SetFont(const MGFontInfo &finfo_in, int i){
  if(i>int(finfos.size()-1))
    return;
  finfos[i] = finfo_in;
}

MGFontInfo FontCache::GetFont(int i){
  MGFontInfo finfo;
  if(finfos.size()<1)
    return finfo;
  if(i>int(finfos.size()-1))
    return finfos[0];
  return finfos[i];
}

MGFontInfo FontCache::GetFont(const std::string &name){
  MGFontInfo finfo;

  std::vector<std::string> name_split = ParseX11Name(name);
  if(name_split.size()<6)
    return finfo;
  
  return GetFont(name_split[0],name_split[1],name_split[2],name_split[3],name_split[4],atoi(name_split[5].c_str()));
}

MGFontInfo FontCache::GetFont(const std::string &foundry_in, const std::string &family,  const std::string &weight,  const std::string &slant,  const std::string &encoding_in, int size, int truetype){

  MGFontInfo finfo;
  std::string foundry=foundry_in;
  std::string encoding=encoding_in;
  if(foundry=="")
    foundry="default-foundry";
  if(encoding=="")
    encoding="iso8859";

  int i;
  if(FontCache::isGlutFont(family)){
    if((i=FontCache::isFaceCached(foundry,family,"medium","r",encoding))>-1){
      return finfos[i];
    }
  }
  //std::cout << "Cached?: " << foundry << "-" << family << "-" << weight << "-" << "-" << slant << "-" << encoding << "-" << size << "\n"; std::cout.flush();
  if((i=FontCache::isCached(foundry,family,weight,slant,encoding,size))>-1){
    //std::cout << "cached "  << i << "\n"; std::cout.flush();
    return finfos[i];
  }else{
    if(slant=="o")
      if((i=FontCache::isCached(foundry,family,weight,"i",encoding,size))>-1)
        return finfos[i];
    if(slant=="i")
      if((i=FontCache::isCached(foundry,family,weight,"o",encoding,size))>-1)
        return finfos[i];
    //std::cout<< "Attempting LoadFont\n"; std::cout.flush();
    if((i=LoadFont(foundry,family,weight,slant,encoding,size,truetype))>-1){
      //std::cout << "Loaded Font " << i << "\n"; std::cout.flush();
      return finfos[i];
      /*
	}else{
	std::cout << "hmm, i comes back as " << i << "\n";
      */
    }
  }
  return finfo;
  
}

int FontCache::isCached(const std::string &foundry, const std::string &family,  const std::string &weight,  const std::string &slant,  const std::string &encoding, int size){

  std::vector<MGFontInfo>::const_iterator finfo=finfos.begin();
  int i=0;
  while(finfo!=finfos.end()){
    if((strcasecmp(finfo->family.c_str(),family.c_str())==0)&&finfo->weight==weight&&finfo->slant==slant&&finfo->fn_size==size){
      if(i>10){
        std::iter_swap(finfos.begin(), finfos.begin()+i);
        return 0;
      }else return i;
    }
    finfo++; i++;
  }
  return -1;
}

int FontCache::isFaceCached(const std::string &foundry, const std::string &family,  const std::string &weight,  const std::string &slant,  const std::string &encoding){

  //std::cout << "Looking for: " << foundry << " " << family << " " << weight << " " << slant << " " << encoding << "\n";
  std::vector<MGFontInfo>::const_iterator finfo=finfos.begin();
  int i=0;
  while(finfo!=finfos.end()){
    //std::cout << "Checking: " << i << " " << finfo->foundry << " " << finfo->family << " " << finfo->weight << " " << finfo->slant << " " << finfo->encoding << " " << finfo->fn_size << " \n";
    if((strcasecmp(finfo->family.c_str(),family.c_str())==0)&&finfo->weight==weight&&finfo->slant==slant){
      //std::cout << "Font Face: " << i << " " << finfo->family << " " << finfo->weight << " " << finfo->slant << " " << finfo->encoding << " " << finfo->fn_size << " matches\n";
      return i;
    }
    finfo++; i++;
  }
  //std::cout << "No match\n";
  return -1;
}


MGFontInfo::MGFontInfo(){
  base = 32;
  last = 256;
  format = -1;
#ifdef USE_FREETYPE
  face = 0;
#endif
}

MGFontInfo::~MGFontInfo(){
}

unsigned char *MGFontInfo::BitMap(WIDE_CHAR_T i) const { 
  if(i-base>=int(bm.size()))
    return BitMap(bm.size()-1);
  return bm[i-base];
}

MGFontInfo::MGFontInfo(const std::vector<std::string> &name, const std::vector<unsigned char *> &bm_in, const std::vector <size_t> &widths_in, const std::vector <size_t> &heights_in, const std::vector <float> &lbearings_in, const std::vector <float> &descents_in, const std::vector <float> &dxs_in, const std::vector <float> &dys_in, WIDE_CHAR_T base_in, WIDE_CHAR_T last_in, int format_in
#ifdef USE_FREETYPE
               ,FT_Face face_in
#endif
		   ){
  if(!(name.size()>1&&FontCache::isGlutFont(name[1]))&&(name.size()==0||name[1].size()==0||widths_in.size()==0||heights_in.size()==0||lbearings_in.size()==0||descents_in.size()==0||dxs_in.size()==0||dys_in.size()==0)){
    return;
  }
  foundry = name[0];
  family = name[1];
  weight = name[2];
  slant = name[3];
  encoding = name[5];
  bm = bm_in;
  fn_size = atoi(name[4].c_str());
  widths = widths_in;
  heights = heights_in;
  lbearings = lbearings_in;
  descents = descents_in;
  dxs = dxs_in;
  dys = dys_in;
  base = base_in;
  last = last_in;
  format = format_in;
#ifdef USE_FREETYPE
  face = face_in;
#endif
}

std::vector<std::string> ParseX11Name(const std::string &name){
  std::vector<std::string> names;

  size_t t = name.find('-');
  if(t==std::string::npos){
     return names;
  }
  std::string foundry = name.substr(1,name.substr(1).find('-'));
  names.push_back(foundry);

  std::string rest = name.substr(name.substr(1).find('-')+1);
  std::string family = rest.substr(1,rest.substr(1).find('-'));
  names.push_back(family);

  rest = rest.substr(rest.substr(1).find('-')+1);
  std::string weight = rest.substr(1,rest.substr(1).find('-'));
  names.push_back(weight);
  
  rest = rest.substr(rest.substr(1).find('-')+1);
  std::string slant = rest.substr(1,rest.substr(1).find('-'));
  names.push_back(slant);

  rest = rest.substr(rest.substr(1).find('-')+1);
  rest = rest.substr(rest.substr(1).find('-')+1);

  rest = rest.substr(rest.substr(1).find('-')+1);
  std::string size = rest.substr(1,rest.substr(1).find('-'));
  names.push_back(size);

  for(int i=0;i<6;i++)
    rest = rest.substr(rest.substr(1).find('-')+1);
  rest = rest.substr(1);
  names.push_back(rest);

  /*
  for(unsigned i=0;i<names.size();i++)
    std::cout << names[i] << "\n";
  */
  return names;
}

std::vector<std::string> ParseX11Name(const char *name){
  return ParseX11Name(std::string(name));
}
