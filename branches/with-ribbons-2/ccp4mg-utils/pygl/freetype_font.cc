/*
     pygl/freetype_font.cc: CCP4MG Molecular Graphics Program
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


#ifdef USE_FREETYPE
#include "freetype_font.h"
#include "font_info.h"

#include <mgutil.h>

#include <ctype.h>
#include <sys/types.h>
#include <dirent.h>
#include <regex.h>

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <stdio.h>
#include <string.h>

#ifdef USE_GLX
extern char **GetX11FontPath(int *npaths);
extern void FreeX11FontPath(char **paths);
#endif

#include <ft2build.h>
#include FT_FREETYPE_H
#include FT_ERRORS_H
#include FT_MULTIPLE_MASTERS_H                                   
#include FT_GLYPH_H                                              
#include FT_SFNT_NAMES_H
#include FT_TRUETYPE_IDS_H
#include FT_TRUETYPE_TABLES_H

#ifdef _USE_DL_FREETYPE_
#include "freetype_dl.h"
#if defined (linux)
#define LIBFT_SOLIBRARY "libfreetype.so.6"
#elif defined (__APPLE_CC__)
#define LIBFT_SOLIBRARY "/usr/X11R6/lib/libfreetype.dylib"
#elif  defined(_WIN32) || defined(__WIN32__)
#define LIBFT_SOLIBRARY "libfreetype.dll"
#else
#define LIBFT_SOLIBRARY "libfreetype.so"
#endif
static int have_ft=0;
#endif

FT_Library library;
static int library_initialized=0;

FT_Error Freetype_Set_Pixel_Sizes( FT_Face face, int fn_size1, int fn_size2 ){
  return FT_Set_Pixel_Sizes( face, fn_size1, fn_size2 );
}
FT_Error Freetype_Get_Kerning(FT_Face face, int left, int right, FT_Vector *kerning){
  return FT_Get_Kerning( face,  FT_Get_Char_Index( face, left), FT_Get_Char_Index( face, right),  ft_kerning_default, kerning );
}

int InitFreetypeLibrary(void){
#ifdef _USE_DL_FREETYPE_
  if(!have_ft){
    if(init_ft(LIBFT_SOLIBRARY)){
      have_ft = 1;
    }
  }
  if(!have_ft){
    printf("Freetype library not found\n");
    return 0;
  }
#endif
  int error = FT_Init_FreeType( &library ); 
  if ( error ){
  }
  if(!error) library_initialized=1;
  return error;
}

#if defined(_WIN32) || defined(__WIN32__)
#define ENVVARSEPARATOR ";"
#else
#define ENVVARSEPARATOR ":"
#endif

#if defined (__APPLE_CC__)
#include <string>
#include <Carbon/Carbon.h>
#include FT_MAC_H

OSStatus GetFontFamilyResource(FMFont font, Handle* oHandle, SInt16 *orsrcFRefNum, Str255 fontFamilyName) {

    InitFreetypeLibrary();

    SInt16 rsrcFRefNum;
    Handle rsrcHandle;
    FSSpec rsrcFSSpec;
    FSRef rsrcFSRef;
    HFSUniStr255 forkName;
    OSStatus status;

    rsrcFRefNum = -1;
    rsrcHandle = 0;
    status = noErr;


    /* Get the font family name to use with the Resource
        Manager when grabbing the 'FOND' resource. */

    status = FMGetFontContainer(font, &rsrcFSSpec);
    require(status == noErr, FMGetFontContainer_Failed);

        /* Open the resource fork of the file. */
    rsrcFRefNum = FSpOpenResFile(&rsrcFSSpec, fsRdPerm);

        /* If the font is based on the ".dfont" file format,
        we need to open the data fork of the file. */
    if ( rsrcFRefNum == -1 ) {

            /* The standard fork name is required to open
            the data fork of the file. */
        status = FSGetDataForkName(&forkName);
        require(status == noErr, FSGetDataForkName_Failed);

            /* The file specification (FSSpec) must be converted
            to a file reference (FSRef) to open the data fork
            of the file. */
        status = FSpMakeFSRef(&rsrcFSSpec, &rsrcFSRef);
        require(status == noErr, FSpMakeFSRef_Failed);

        status = FSOpenResourceFile(&rsrcFSRef,
            forkName.length, forkName.unicode,
            fsRdPerm, &rsrcFRefNum);
        require(status == noErr, FSOpenResourceFile_Failed);
    }

    UseResFile(rsrcFRefNum);

        /* On Mac OS X, the font family identifier may not
        match the resource identifier after resolution of
        conflicting and duplicate fonts. */

    rsrcHandle = Get1NamedResource(FOUR_CHAR_CODE('FOND'), fontFamilyName);

FSOpenResourceFile_Failed:
FSpMakeFSRef_Failed:
FSGetDataForkName_Failed:
FMGetFontContainer_Failed:

    if ( oHandle != 0 ){
        *oHandle = rsrcHandle;
        *orsrcFRefNum = rsrcFRefNum;
    }

    return status;
}

#endif

#define GETLINESIZE 1024

std::string EatWhiteSpace(const std::string &in){
  std::string ret;
  std::string::const_iterator s = in.begin();
  while(isspace(*s)&&s!=in.end()){
    s++;
  }
  while(s!=in.end()){
    ret += *s;
    s++;
  }
  return ret;
}

std::vector<std::string> Getttfontdirs(){

  std::vector<std::string> ttfontdirs;

  char *ttfontdirenvcp = getenv("TTFONTDIR");
  if(ttfontdirenvcp){

    std::string ttfontdirenv =  std::string(ttfontdirenvcp);
    int i;
    if((i=ttfontdirenv.find(ENVVARSEPARATOR))>-1){
      ttfontdirs.push_back(ttfontdirenv.substr(0,i));
      ttfontdirenv = ttfontdirenv.substr(i+1,ttfontdirenv.size());
      while((i=ttfontdirenv.find(ENVVARSEPARATOR))>-1){
        ttfontdirs.push_back(ttfontdirenv.substr(0,i));
        ttfontdirenv = ttfontdirenv.substr(i+1,ttfontdirenv.size());
      }
      ttfontdirs.push_back(ttfontdirenv);
    }
    else
      ttfontdirs.push_back(ttfontdirenv); 
  }
  
#ifdef USE_GLX
  int npaths;
  char **fp = GetX11FontPath(&npaths);
  for(int i=0;i<npaths;i++){
    if(fp[i]&&strlen(fp[i])>0){
      if(fp[i][0]=='/'){
        ttfontdirs.push_back(std::string(fp[i])); 
      }else{
        /* Look for font server config file */
        std::ifstream fsfile("/usr/X11R6/lib/X11/fs/config");

        char buf[GETLINESIZE];

        bool finished = false;

        std::string outstr;

        while(fsfile.good()&&!finished){
          fsfile.getline(buf,GETLINESIZE);
          if(buf){
            std::string bufstr(buf);
            if(bufstr.size()){
              if(bufstr.substr(0,9) == "catalogue"&&bufstr.find("=")!=std::string::npos){
                int eq = bufstr.find("=");
                bufstr = bufstr.substr(eq+1);
                bufstr = EatWhiteSpace(bufstr);
                if(bufstr[0]=='/'&&bufstr.find(":unscaled")==std::string::npos){
                  if(bufstr[bufstr.size()-1]==',')
                    outstr = bufstr.substr(0,bufstr.size()-1);
                  else
                    outstr = bufstr;
                  ttfontdirs.push_back(outstr); 
                }
                while(fsfile.good()&&bufstr.find(",")!=std::string::npos){
                  fsfile.getline(buf,GETLINESIZE);
                  bufstr = std::string(buf);
                  bufstr = EatWhiteSpace(bufstr);
                  if(bufstr[0]=='/'&&bufstr.find(":unscaled")==std::string::npos){
                    if(bufstr[bufstr.size()-1]==',')
                      outstr= bufstr.substr(0,bufstr.size()-1);
                    else
                      outstr = bufstr;
                    ttfontdirs.push_back(EatWhiteSpace(outstr)); 
                  }
                }
                finished = true;
              }
            }
          }
        }
        fsfile.close();

      }
    }
  }
  FreeX11FontPath(fp);
#endif

  std::map<std::string,int> fontmap;
  int cur_idx = 1;
  for(unsigned i=0;i<ttfontdirs.size();i++)
    if(fontmap[ttfontdirs[i]]==0)
      fontmap[ttfontdirs[i]] = cur_idx++;
  std::vector<std::string> ttfontdirs_uniq(fontmap.size());
  std::map<std::string,int>::iterator map_iter = fontmap.begin();
  while(map_iter!=fontmap.end()){
    ttfontdirs_uniq[map_iter->second-1] = map_iter->first;
    map_iter++;
  }
  ttfontdirs = ttfontdirs_uniq;

  return ttfontdirs;
}

std::string GetT1PostScriptFamilyName(const char* fontName){

  const char * const kFontName = "/FamilyName";
  unsigned int headerSize = 5000;
  char buffer[5000];
  FILE * fontFILE = fopen(fontName, "rb");
  if(!fontFILE){
    std::cout << "error opening file " << fontName << " in GetT1PostScriptFamilyName\n"; std::cout.flush();
    return std::string("");
  }

  fread(buffer, 1, headerSize, fontFILE);
  char * tok = NULL;

  for(unsigned int i = 0; i < headerSize; ++i){
    if(strncmp(buffer+i, kFontName, strlen(kFontName)) == 0){
      tok = strtok(buffer + i + strlen(kFontName), ")\r\n\t");
      break;
    }
  }
  fclose(fontFILE);
  if(tok){
    std::string tokstr = std::string(tok);
    unsigned int pos = tokstr.find('(');
    if(pos != std::string::npos && tokstr.length()>1){
      return tokstr.substr(pos+1);
    }
  }
  return std::string("");
}

const char *GetFreeTypeWeight (TT_OS2 *os2){
    switch (os2->usWeightClass) {
    case 100: return "thin";
    case 200: return "extra light";
    case 300: return "light";
    case 400: return "medium";
    case 500: return "medium";
    case 600: return "semi bold";
    case 700: return "bold";
    case 800: return "extra bold";
    case 900:  return "black";
    }

    return 0;
}

const char *ExtractFreeTypeName (FT_UInt name_len, FT_Byte* name, bool unicode)
{
    int          name_pos = 0;
    FT_UInt      i;
    int          increment;
    char *name_buffer = new char[1024];

    if (unicode) {
        name_len =  name_len < 512 ? name_len : 512;
        increment = 2;
        i = 1;
    } else {
        name_len =  name_len < 256 ? name_len : 256;
        increment = 1;
        i = 0;
    }

    for (; i < name_len; i += increment) {
        if (!isascii (name[i])) {
            name[i] = '$';
        }
        name_buffer[name_pos++] = (name[i] == '-') ? '_' : name[i];
    }

    name_buffer[name_pos] = 0;

    return name_buffer;
}

void Tokenize(const std::string &str,std::vector<std::string> &tokens,
                      const std::string &delimiters = " "){

    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

int cmp_nocase(const std::string &s1, const std::string &s2){
//======================================================================
// compare two strings without caring about case
// taken from Stroustrup Special Ed, 20.3.8
//======================================================================
  std::string::const_iterator p1=s1.begin();
  std::string::const_iterator p2=s2.begin();
  while (p1!=s1.end() && p2!=s2.end()) {
    if (toupper(*p1) != toupper(*p2))
      return (toupper(*p1)<toupper(*p2)) ? -1 : 1;
    ++p1; ++p2;
  }
  return (s2.size()==s1.size()) ? 0 : (s1.size()<s2.size()) ? -1:1;
}

int match_family(const char *test_family, const char *req_family){

  /* The current family we are checking */
  std::vector<std::string> test_tokens;
  Tokenize(std::string(test_family), test_tokens);

  /* The family we have asked for */
  std::vector<std::string> req_tokens;
  Tokenize(std::string(req_family), req_tokens);

  int match_ok = 0;
  for(unsigned i=0;i<req_tokens.size();i++){
    match_ok = 0;
    for(unsigned j=0;j<test_tokens.size();j++){
      if(req_tokens[i]==test_tokens[j]){
	//if(cmp_nocase(req_tokens[i],test_tokens[j])==0){
	match_ok = 1;
      }
    }
    if(!match_ok) return 0;
  }
  return match_ok;

}

int match_regexp(const char *string, const char *pattern){
    int    status;
    regex_t    re;

    //printf("%s\n",string);
    //printf("%s\n",pattern);

    if (regcomp(&re, pattern, REG_EXTENDED|REG_NOSUB|REG_ICASE) != 0) {
        return(0);      /* Report error. */
    }
    status = regexec(&re, string, (size_t) 0, NULL, 0);
    regfree(&re);
    if (status != 0) {
        return 0;      /* Report error. */
    }
    return 1;
}

const char *FindFreeTypeFontPath(const char* fn, const char *weight, const char *slant){

  std::vector<std::string> ttfontdirs = Getttfontdirs();
  if(ttfontdirs.size()==0)
    return 0;

  DIR *dirp;
  std::vector<std::string>::iterator dir=ttfontdirs.begin();
  while(dir!=ttfontdirs.end()){
    if(*dir=="")
      dirp = opendir(".");
    else
      dirp = opendir(dir->c_str());

    if(dirp){
      struct dirent *dirent;
      while( (dirent = readdir(dirp)) != NULL){
	std::string fstr = std::string(dirent->d_name);
	if(fstr.size()>3&&fstr.substr(fstr.size()-4,fstr.size()-1)==std::string(".ttf")
	  ||fstr.size()>3&&fstr.substr(fstr.size()-4,fstr.size()-1)==std::string(".otf")
	  ||fstr.size()>3&&fstr.substr(fstr.size()-4,fstr.size()-1)==std::string(".pfa")
	  ||fstr.size()>3&&fstr.substr(fstr.size()-4,fstr.size()-1)==std::string(".pfb")){
#if defined(_WIN32) || defined(__WIN32__)
	  std::string fullfstr = *dir+std::string("\\")+fstr;
#else
	  std::string fullfstr = *dir+std::string("/")+fstr;
#endif
	  FT_Library library;
	  FT_Face face;
	  int error = FT_Init_FreeType( &library ); 
	  if ( error ){
	    printf("Error initializing freetype library\n");
	  }
	  error = FT_New_Face( library, fullfstr.c_str(), 0, &face );
          /* Should probably use face->family_name here as well. */
	    
	  int names = FT_Get_Sfnt_Name_Count (face);
	  
	  for(int i=0;i<names;i++){
	    FT_SfntName aname;
	    error = FT_Get_Sfnt_Name( face, i, &aname );
	    if ((aname.name_id == TT_NAME_ID_FONT_FAMILY) || (aname.name_id == TT_NAME_ID_PS_NAME)){
	      if ( error ){
		printf("Error getting name %d\n",i);
	      }
	      TT_OS2 *os2 = (TT_OS2 *) FT_Get_Sfnt_Table(face, ft_sfnt_os2);
	      if(match_regexp((const char*)aname.string,fn)){
		const char *name = ExtractFreeTypeName(aname.string_len,aname.string,false);
                if(name){
		if(match_family(name,fn)&&match_regexp(weight,GetFreeTypeWeight(os2))){
		  if(os2->fsSelection&(1 << 0)&&(!strncmp(slant,"i",1)||!strncmp(slant,"o",1))
		     ||!(os2->fsSelection&(1 << 0))&&(strncmp(slant,"i",1)&&strncmp(slant,"o",1))){
		    //printf("index: %d of %d\n",i,names);
		    //printf("%s -- %s -- %s\n",name,GetFreeTypeWeight(os2),slant);
		    //printf("%s\n",fullfstr.c_str());
		    delete name;
                    closedir(dirp);
                    if(face) FT_Done_Face(face);
		    return fullfstr.c_str();
		  }
		}
                }
		delete name;
	      }
	    }
	  }
	  if(fstr.size()>3&&fstr.substr(fstr.size()-4,fstr.size()-1)==std::string(".pfa")||
	     fstr.size()>3&&fstr.substr(fstr.size()-4,fstr.size()-1)==std::string(".pfb")){
	    const char *psname = FT_Get_Postscript_Name(face);
	    if(psname&&strlen(psname)){
	      std::string psnamestr = std::string(psname);
	      std::string psfamily;
	      std::string psweight;
	      std::string psslant;
	      unsigned int pos = psnamestr.find('-');
	      //psfamily = psnamestr;
	      psfamily = GetT1PostScriptFamilyName(fullfstr.c_str());
	      psweight = "medium";
	      psslant = "r";
	      if(pos != std::string::npos){
		psfamily = psnamestr.substr(0,psnamestr.find('-'));
		std::string psweightslant = psnamestr.substr(psnamestr.find('-')+1);
		unsigned int posbold = psweightslant.find("Bold");
		if(posbold != std::string::npos)
		  psweight = "bold";
		unsigned int posital = psweightslant.find("Ital");
		unsigned int posobli = psweightslant.find("Obli");
		if(posital != std::string::npos || posobli != std::string::npos)
		  psslant = "i";
	      }
	      //if(match_regexp(psfamily.c_str(),fn)&&match_regexp(psweight.c_str(),weight)&&match_regexp(psslant.c_str(),slant)){
	      if(match_family(GetT1PostScriptFamilyName(fullfstr.c_str()).c_str(),fn)&&match_regexp(psweight.c_str(),weight)&&match_regexp(psslant.c_str(),slant)){
		std::string afm = fullfstr.substr(0,fullfstr.rfind('.'))+std::string(".afm");
		std::string pfm = fullfstr.substr(0,fullfstr.rfind('.'))+std::string(".pfm");
		FT_Attach_File( face, afm.c_str() );
		FT_Attach_File( face, pfm.c_str() );
		/*
		if(FT_HAS_KERNING(face))
		  printf("Have kerning for %s\n",fstr.c_str());
		else
		  printf("Do not have kerning for %s\n",fstr.c_str());
		*/
                closedir(dirp);
                if(face) FT_Done_Face(face);
		return fullfstr.c_str();
	      }
	    }
	  }
          if(face) FT_Done_Face(face);
	}
      }
    }
    
    closedir(dirp);
    dir++;
  }
  return 0;
}

static inline FT_BitmapGlyph LoadFreetypeFontChar(const FT_Face face, const wchar_t &c){

  int error;
  int glyph_index;

  FT_Glyph        glyph;

  glyph_index = FT_Get_Char_Index(face,  (face->charmap->encoding == ft_encoding_symbol) ? c | 0xf000: c);

  error = FT_Load_Glyph( face, glyph_index, 0 );
  if ( error ){
    printf("Error loading glyph\n");
    return 0;
  }
  
  error = FT_Get_Glyph( face->glyph, &glyph );
  if ( error ){
    printf("Error getting glyph\n");
    return 0;
  }

  FT_Vector origin;

  error = FT_Glyph_To_Bitmap(&glyph,ft_render_mode_normal,&origin,1);
  if ( error ){
    //printf("Warning:cannot get glyph for \"%c\" in FT_Glyph_To_Bitmap\n",c);
    return 0;
  }

  return (FT_BitmapGlyph)glyph;

}

FT_BitmapGlyph LoadFreetypeFontChar(const FT_Face face, const int &width, const int &height, const wchar_t &c){
  int error;
  error = FT_Set_Pixel_Sizes( face, width, height );
  if ( error ){
    printf("Error setting size\n");
    return 0;
  }
  return LoadFreetypeFontChar(face,c);
}

MGFontInfo LoadFreeTypeFont(const std::string &foundry, const std::string &family,  const std::string &weight,  const std::string &slant,  const std::string &encoding, int size){

  int font_index;
  MGFontInfo finfo;

  if((font_index=FontCache::isCached(foundry,family,weight,slant,encoding,size))!=-1)
    return FontCache::GetFont(font_index);

  /* else do something ...*/
  int error=0;
  if(!library_initialized){
    error = InitFreetypeLibrary();
    if(error) return finfo;
  }  

  if((font_index=FontCache::isFaceCached(foundry,family,weight,slant,encoding))!=-1){
    MGFontInfo ftmp = FontCache::GetFont(font_index);
    //std::cout << "got index: " << font_index << "\n";
    if(ftmp.face){
      //std::cout << "Returning load of cached face\n"; std::cout.flush();
      return LoadFreeTypeFont(ftmp.face,ftmp.foundry,ftmp.family,ftmp.weight,ftmp.slant,ftmp.encoding,size);
    }//else
      //std::cout << "Not returning load of cached face\n"; std::cout.flush();
  }
  /*
    Need to convert family, etc to font path.
  */
  const char *fontpath = FindFreeTypeFontPath(family.c_str(),weight.c_str(),slant.c_str());
  return LoadFreeTypeFont(fontpath,size);

}

MGFontInfo LoadFreeTypeFont(const char *fontpath, int size){

  //std::cout << "MGFontInfo LoadFreeTypeFont(const char *fontpath, int size) 1\n"; std::cout.flush();
  MGFontInfo finfo;
  if(!fontpath) {
    return finfo;
  }

  int error=0;
  if(!library_initialized){
    error = InitFreetypeLibrary();
    if(error) return finfo;
  }  

  FT_Face face;
  error = FT_New_Face( library, fontpath, 0, &face );
  if(error) return finfo;
  if(!face) return finfo;
  //std::cout << "face->family_name " << face->family_name << "\n";
  std::string family_name = face->family_name;

  std::string foundry="default-foundry";
  std::string family;
  std::string weight;
  std::string slant;
  std::string encoding="iso8859";
  //std::string fstr = std::string(fontpath);
  
  if(strlen(fontpath)>3&&!strcmp(fontpath+(strlen(fontpath)-4),".pfa")||
     !strcmp(fontpath+(strlen(fontpath)-4),".pfb")){
    //if(fstr.size()>3&&fstr.substr(fstr.size()-4,fstr.size()-1)==std::string(".pfa")||
    // fstr.size()>3&&fstr.substr(fstr.size()-4,fstr.size()-1)==std::string(".pfb")){
    const char *psname = FT_Get_Postscript_Name(face);
    if(psname&&strlen(psname)){
      std::string psnamestr = std::string(psname);
      unsigned int pos = psnamestr.find('-');
      family = GetT1PostScriptFamilyName(fontpath);
      if(family=="") family=family_name;
      //std::cout << "Postscript Name " << family << "\n";
      if(family_name!=family)
        std::cout << "ps and face name not equal!!\n";
      weight = "medium";
      slant = "r";
      if(pos != std::string::npos){
	std::string psweightslant = psnamestr.substr(psnamestr.find('-')+1);
	unsigned int posbold = psweightslant.find("Bold");
	if(posbold != std::string::npos)
	  weight = "bold";
	unsigned int posital = psweightslant.find("Ital");
	unsigned int posobli = psweightslant.find("Obli");
	if(posital != std::string::npos || posobli != std::string::npos)
	  slant = "i";
      }
      std::string fullfstr(fontpath);
      std::string afm = fullfstr.substr(0,fullfstr.rfind('.'))+std::string(".afm");
      std::string pfm = fullfstr.substr(0,fullfstr.rfind('.'))+std::string(".pfm");
      FT_Attach_File( face, afm.c_str() );
      FT_Attach_File( face, pfm.c_str() );
      /*
      char afm[1024];
      char pfm[1024];
      afm[0] = '\0'; pfm[0] = '\0';
      strcat(afm,fontpath);
      strcat(pfm,fontpath);
      afm[strlen(fontpath)-3] = '\0'; pfm[strlen(fontpath)-3] = '\0';
      strcat(afm,"afm");
      strcat(pfm,"pfm");
      FT_Attach_File( face, afm );
      FT_Attach_File( face, pfm );
      */
    }
  }else{
    /* So this loop can go, to be replaced by face->family_name.
     * The os2 stuff for bold/italic will remain and is independent
     * of the aname.
     */
    int names = FT_Get_Sfnt_Name_Count (face);
    for(int i=0;i<names;i++){
      FT_SfntName aname;
      error = FT_Get_Sfnt_Name( face, i, &aname );
      if ((aname.name_id == TT_NAME_ID_FONT_FAMILY) || (aname.name_id == TT_NAME_ID_PS_NAME)){
	if ( error ){
	  fprintf(stderr,"Error getting name %d\n",i);
	}
	TT_OS2 *os2 = (TT_OS2 *) FT_Get_Sfnt_Table(face, ft_sfnt_os2);
	family = std::string(ExtractFreeTypeName(aname.string_len,aname.string,false));
        //std::cout << "TTF Name " << family << "\n";
        //if(family_name!=family)
          //std::cout << "ttf and face name not equal!!\n";
	if(os2) {
	  const char *wp = GetFreeTypeWeight(os2);
	  if(wp){
	    weight = std::string(wp);
	    if(os2->fsSelection&(1 << 0))
	      slant = "i";
	    else
	      slant = "r";
	  } else return MGFontInfo();
	} else return MGFontInfo();
        if(family!=std::string("")&&weight!=std::string(""))
	  break;
      }
    }
  }

  return LoadFreeTypeFont(face,foundry,family_name,weight,slant,encoding,size);

}

MGFontInfo LoadFreeTypeFont(FT_Face face, const std::string &foundry, const std::string &family,
		          const std::string &weight, const std::string &slant, const std::string &encoding,
			  int size){

  MGFontInfo finfo;
  if(!face) return finfo;

  //std::cout << "LoadFreeTypeFont(FT_Face face, ...): " << family << " " << weight << " " << slant << " " << size << "\n"; std::cout.flush();
  int width, height;
  width = height = size;

  int error=0;
  if(!library_initialized){
    error = InitFreetypeLibrary();
    if(error) return finfo;
  }  

  if ( error ){
    printf("Error loading font: %x\n",error);
    return finfo;
  }

  std::vector <unsigned char *> bm;
  std::vector <size_t> widths;
  std::vector <size_t> heights;
  std::vector <float> lbearings;
  std::vector <float> descents;
  std::vector <float> dxs;
  std::vector <float> dys;

  WIDE_CHAR_T base = 32;
  WIDE_CHAR_T last = 256;

  error = FT_Set_Pixel_Sizes( face, width, height );
  if ( error ){
    printf("Error setting size\n");
  }

  for(int e=0; e < face->num_charmaps; e++) {
    //fprintf(stderr, "found encoding pid=%d eid=%d\n", face->charmaps[e]->platform_id, face->charmaps[e]->encoding_id);
    if( FT_Set_Charmap(face, face->charmaps[e]) ) {
      fprintf(stderr, "**** Cannot set charmap in FreeType ****\n");
    }
  }

  int nerrors = 0;
  for(WIDE_CHAR_T c=base;c<=last;c++){
    if(size==0){
      unsigned char *bmtmp = new unsigned char[0];
      bm.push_back(bmtmp);
      widths.push_back(0);
      heights.push_back(0);
      lbearings.push_back(0.0);
      descents.push_back(0.0);
      dxs.push_back(0.0);
      dys.push_back(0.0); 
      break;
    }
    FT_BitmapGlyph bitmap = LoadFreetypeFontChar(face,c);
    if(bitmap){
      unsigned char *bmtmp = new unsigned char[bitmap->bitmap.width*bitmap->bitmap.rows];
      for(int ii=0;ii<bitmap->bitmap.rows;ii++)
        for(int jj=0;jj<bitmap->bitmap.width;jj++)
	  bmtmp[(bitmap->bitmap.rows-1-ii)*bitmap->bitmap.width + jj] = bitmap->bitmap.buffer[ii*bitmap->bitmap.width + jj];
      bm.push_back(bmtmp);
      widths.push_back(bitmap->bitmap.width);
      heights.push_back(bitmap->bitmap.rows);
      lbearings.push_back((float)(face->glyph->metrics.horiBearingX/64.0));
      descents.push_back((float)(face->glyph->metrics.height/64.0 - face->glyph->metrics.horiBearingY/64.0));
      dxs.push_back((float)(face->glyph->advance.x/64.0));
      dys.push_back(0.0);
    }else{
      nerrors++;
      unsigned char *bmtmp = new unsigned char[0];
      bm.push_back(bmtmp);
      widths.push_back(0);
      heights.push_back(0);
      lbearings.push_back(0.0);
      descents.push_back(0.0);
      dxs.push_back(0.0);
      dys.push_back(0.0);
    }
  }
  if(nerrors>5)
    return MGFontInfo();

  std::vector<std::string> name;
  name.push_back(foundry);
  name.push_back(family);
  name.push_back(weight);
  name.push_back(slant);
  name.push_back(IntToString(size));
  name.push_back(encoding);

  return MGFontInfo(name,bm,widths,heights,lbearings,descents,dxs,dys,base,last,FONT_BITMAP_GRAYSCALE_FT,face);

}

MGFontInfo LoadFreeTypeFont(const char *name){
  MGFontInfo finfo;
  std::vector<std::string> name_split = ParseX11Name(std::string(name));
  if(name_split.size()<5)
    return finfo;
  return LoadFreeTypeFont(name_split[0],name_split[1],name_split[2],name_split[3],name_split[5],atoi(name_split[4].c_str()));
}

int LoadAllFreeTypeFonts(){

#if defined (__APPLE_CC__)
  FT_Face face;
  OSStatus status;
  OptionBits options = kFMUseGlobalScopeOption;

  FMFontFamilyIterator  famIter;
  status = FMCreateFontFamilyIterator( NULL, NULL, options, &famIter );
  FMFontFamily          family   = NULL;
  FT_Long face_index;

  face_index = 0;
  while( status == 0  ) {
    status = FMGetNextFontFamily( &famIter, &family );
    if( status == 0 ) {
      int                           stat2;
      FMFontFamilyInstanceIterator  instIter;
      Str255                        famNameStr;
      char                          famName[256];


      /* get the family name */
      FMGetFontFamilyName( family, famNameStr );
      CopyPascalStringToC( famNameStr, famName );

      /* iterate through the styles */
      FMCreateFontFamilyInstanceIterator( family, &instIter );

      face_index = 0;
      stat2 = 0;
      std::map<std::string,int> mymap;
      mymap.clear();
      while( stat2 == 0  ) {
        FMFontStyle  style;
        FMFontSize   size;
        FMFont       font;


        stat2 = FMGetNextFontFamilyInstance( &instIter, &font,
                                              &style, &size );
        if( stat2 == 0 && size == 0 ) {
          char fullName[256];

          /* build up a complete face name */
          strcpy( fullName, famName );
          int isa = 1;
          for(unsigned i=0;i<strlen(famName);i++){
            if(!isascii(famName[i]))
              isa = 0;
          }
          if(isa){
          std::string family = std::string(fullName);
          std::string weight = "medium";
          std::string slant = "r";
          if ( style & bold ){
            strcat( fullName, " Bold" );
            weight = "bold";
          }
          if ( style & italic ){
            strcat( fullName, " Italic" );
            slant = "i";
          }
          /*
          FSSpec   pathSpec;
          FMGetFontContainer( font, &pathSpec );
          int error = FT_New_Face_From_FSSpec(library,&pathSpec,face_index,&face);
          if(!error){
            const char *psname = FT_Get_Postscript_Name(face);
            printf("psname: %s\n",psname);
            mymap[std::string(fullName)] = face_index+1;
          }else printf("error: --%d--\n",error);
          */
          Handle fond=0;
          SInt16 resfile;
          OSStatus stat3 = GetFontFamilyResource(font,&fond,&resfile,famNameStr);
          //printf("Trying fullName: %s %d %d\n",fullName,face_index,mymap[std::string(fullName)]); fflush(stdout);
	  bool stuffed = (style==0&&face_index>0);
          if(!stuffed&&stat3==0&&fond&&mymap[std::string(fullName)]==0){
            int error = FT_New_Face_From_FOND(library,fond,face_index,&face);
            CloseResFile(resfile);
            ReleaseResource(fond);
            if(!error){
              MGFontInfo finfo; 
              std::string foundry = "default-foundry";
              std::string encoding = "iso8859-1";
              finfo = LoadFreeTypeFont(face,foundry,family,weight,slant,encoding,0);
	      if(finfo.foundry.size()>0&&finfo.family.size()>0&&finfo.widths.size()>0&&finfo.heights.size()>0&&finfo.lbearings.size()>0&&finfo.descents.size()>0&&finfo.dxs.size()>0&&finfo.dys.size()>0){
	        if(FontCache::isFaceCached(finfo.foundry,finfo.family,finfo.weight,finfo.slant,finfo.encoding)==-1){
                  //std::cout << "Mac Font " << family << "\n";
                  FontCache::AddFont(finfo);
                }
              }
              const char *psname = FT_Get_Postscript_Name(face);
              //printf("psname: %s\n",psname);
              mymap[std::string(fullName)] = face_index+1;
            }
          } //else {std::cout << "Skipping...\n"; std::cout.flush(); }
        }
      }
      face_index++;
    }

    FMDisposeFontFamilyInstanceIterator( &instIter );
    }
  }

  FMDisposeFontFamilyIterator( &famIter );

#endif
  
  std::vector<std::string> ttfontdirs = Getttfontdirs();
  if(ttfontdirs.size()==0)
    return 0;

  MGFontInfo finfo;
  
  DIR *dirp;
  std::vector<std::string>::iterator dir=ttfontdirs.begin();
  while(dir!=ttfontdirs.end()){
    if(*dir=="")
      dirp = opendir(".");
    else
      dirp = opendir(dir->c_str());

    if(dirp){
      struct dirent *dirent;
      while( (dirent = readdir(dirp)) != NULL){
	std::string fstr = std::string(dirent->d_name);
	if(fstr.size()>3&&fstr.substr(fstr.size()-4,fstr.size()-1)==std::string(".ttf")
	  ||fstr.size()>3&&fstr.substr(fstr.size()-4,fstr.size()-1)==std::string(".otf")
	  ||fstr.size()>3&&fstr.substr(fstr.size()-4,fstr.size()-1)==std::string(".pfa")
	  ||fstr.size()>3&&fstr.substr(fstr.size()-4,fstr.size()-1)==std::string(".pfb")){
#if defined(_WIN32) || defined(__WIN32__)
	  std::string fullfstr = *dir+std::string("\\")+fstr;
#else
	  std::string fullfstr = *dir+std::string("/")+fstr;
#endif
	  finfo = LoadFreeTypeFont(fullfstr.c_str(),0);
	  /*
	  std::cout << "Consider:\n";
          std::cout << finfo.foundry << "\n";
          std::cout << finfo.family << "\n";
          std::cout << finfo.weight << "\n";
          std::cout << finfo.slant << "\n";
          std::cout << finfo.encoding << "\n";
	  */
	  if(finfo.foundry.size()>0&&finfo.family.size()>0&&finfo.widths.size()>0&&finfo.heights.size()>0&&finfo.lbearings.size()>0&&finfo.descents.size()>0&&finfo.dxs.size()>0&&finfo.dys.size()>0){
	    if(FontCache::isFaceCached(finfo.foundry,finfo.family,finfo.weight,finfo.slant,finfo.encoding)==-1){
	      /*
	      std::cout << "Adding:\n";
	      std::cout << finfo.foundry << "\n";
	      std::cout << finfo.family << "\n";
	      std::cout << finfo.weight << "\n";
	      std::cout << finfo.slant << "\n";
	      std::cout << finfo.encoding << "\n";
	      */
	      FontCache::AddFont(finfo);
	    }
	    /*
	  }else{
	  std::cout << "Not added\n";
	    */
	  }
	}
      }
    }
    closedir(dirp);
    dir++;
  }
  return 0;
}

int FindNextFreeTypeFontSize(const MGFontInfo &finfo, int step){
  if(!finfo.isValid()) return 0;
  return finfo.Size()+step;
}


#endif
