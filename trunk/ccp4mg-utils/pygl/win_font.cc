/*
     pygl/win_font.cc: CCP4MG Molecular Graphics Program
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


#if defined(_WIN32) || defined(__WIN32__)

#include <windows.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>

#include "win_font.h"
#include "font_info.h"

#include <mgutil.h>
#include <ppmutil.h>

std::vector<std::string> font_families;

BOOL CALLBACK EnumFamCallBack(LPLOGFONT lplf, LPNEWTEXTMETRIC lpntm, DWORD FontType, LPVOID aFontCount) 
{ 

    int far * aiFontCount = (int far *) aFontCount; 
 
    // Record the number  TrueType
    // fonts in the font_families array. 
 
    if (FontType & RASTER_FONTTYPE) {
        aiFontCount[0]++; 
    }
    else if (FontType & TRUETYPE_FONTTYPE) {
        aiFontCount[2]++; 
	font_families.push_back(std::string(lplf->lfFaceName));
    }
    else {
        aiFontCount[1]++; 
    }

    if (aiFontCount[0] || aiFontCount[1] || aiFontCount[2]) 
        return TRUE; 
    else 
        return FALSE; 

    UNREFERENCED_PARAMETER( lpntm ); 
} 
unsigned char* GetPixelData(HBITMAP hBMP, HDC hDC, int height){ 
/*
 16bpp (555)   The blue mask is 0x001F, the green mask is 0x03E0, and the red mask is 0x7C00.  
 16bpp (565)   The blue mask is 0x001F, the green mask is 0x07E0, and the red mask is 0xF800.  
 32bpp (8888*) The blue mask is 0x000000FF, the green mask is 0x0000FF00, and the red mask is 0x00FF0000 
 */
    LPBITMAPINFO pbi=NULL;

    if ((pbi = (LPBITMAPINFO)(new char[sizeof(BITMAPINFOHEADER) +
        256 * sizeof(RGBQUAD)])) == NULL) { 
          std::cout << "pbi alloc\n"; std::cout.flush(); return NULL; }

    unsigned char *lpvBits=0;

    ZeroMemory(&pbi->bmiHeader, sizeof(BITMAPINFOHEADER));
    pbi->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);

    pbi->bmiHeader.biCompression = BI_RGB;
    pbi->bmiHeader.biSizeImage = 0;

    // Get info but first de-select hBMP because GetDIBits requires it:
    HBITMAP hbmp_tmp=CreateCompatibleBitmap(hDC,(height+7)/8,height);
    HGDIOBJ obj = SelectObject(hDC, hbmp_tmp);

    if (!GetDIBits(hDC, hBMP, 0, height, NULL, pbi, DIB_RGB_COLORS))
         { std::cout << "getdibits 1\n"; std::cout.flush();  
          delete pbi; return NULL; }


    // Reserve memory for bitmap bits:
    if ((lpvBits = new unsigned char[pbi->bmiHeader.biSizeImage]) == NULL)
         { std::cout << "lpvBits alloc\n"; std::cout.flush(); 
            delete pbi; return NULL; }

    // Have GetDIBits convert hBMP to a DIB (device-independent bitmap):
    if (!GetDIBits(hDC, hBMP, 0, height, lpvBits, pbi,
        DIB_RGB_COLORS))  { std::cout << "getdibits 2\n"; std::cout.flush(); 
        delete pbi; delete lpvBits;  return NULL; }

    DeleteObject(obj);

    unsigned char *pixels = new unsigned char[pbi->bmiHeader.biSizeImage];
    if (pbi->bmiHeader.biBitCount == 32) {
      int ii=0;
      for(int i=0;i<pbi->bmiHeader.biHeight; i++){
        for(int j=0;j<pbi->bmiHeader.biWidth; j++){
           pixels[ii++] = 255-lpvBits[(i*pbi->bmiHeader.biWidth+j)*4];
        }
      }
    }else if(pbi->bmiHeader.biBitCount == 24) {
      int ii=0;
      for(int i=0;i<pbi->bmiHeader.biHeight; i++){
        for(int j=0;j<pbi->bmiHeader.biWidth; j++){
           pixels[ii++] = 255-lpvBits[(i*pbi->bmiHeader.biWidth+j)*3];
        }
      }
    }else if(pbi->bmiHeader.biBitCount == 16) {
      int ii=0;
      for(int i=0;i<pbi->bmiHeader.biHeight; i++){
        for(int j=0;j<pbi->bmiHeader.biWidth; j++){
           pixels[ii++] = 255-lpvBits[(i*pbi->bmiHeader.biWidth+j)*2];
        }
      }
    }
    delete pbi;
    delete lpvBits;
    return pixels;
}


HFONT GetHFONT(const std::string &family, const std::string &weight, const std::string &slant, const int size){

  int fw; 
  BOOLEAN italic;

  if(weight=="bold")
    fw = FW_BOLD;
  else
    fw = FW_NORMAL;

  if(slant=="i"||slant=="o")
    italic = TRUE;
  else
    italic = FALSE;

  int antialias_flag = ANTIALIASED_QUALITY;

  OSVERSIONINFO osinfo;
  osinfo.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
  GetVersionEx(&osinfo);

  if(osinfo.dwMajorVersion==5 && osinfo.dwMinorVersion > 0) antialias_flag = 5;

  HFONT font=0;
  font = CreateFont(-size,			//Height
                      0,			// Width Of Font
                      0,			// Angle Of Escapement
                      0,			// Orientation Angle
                      fw,			// Font Weight
                      italic,			// Italic
                      FALSE,			// Underline
                      FALSE,			// Strikeout
                      ANSI_CHARSET,		// Character Set Identifier
                      OUT_TT_PRECIS,		// Output Precision
                      CLIP_DEFAULT_PRECIS,	// Clipping Precision
                      antialias_flag,	// Output Quality		
                      FF_DONTCARE|DEFAULT_PITCH,// Family And Pitch
                      family.c_str());		// Font Name
  return font;

}

MGFontInfo LoadWinFont(const std::string &foundry, const std::string &family,  const std::string &weight,  const std::string &slant,  const std::string &encoding, int isize, HDC hdc, HDC hMemDC){

  MGFontInfo finfo;

  TEXTMETRIC tm;
  ABC abc[256];

  std::vector <unsigned char *> bm;
  std::vector <size_t> widths;
  std::vector <size_t> heights;
  std::vector <float> lbearings;
  std::vector <float> descents;
  std::vector <float> dxs;
  std::vector <float> dys;
  HFONT font = GetHFONT(family,weight,slant,isize);
  if(font){
    HGDIOBJ obj = SelectObject(hMemDC, font);
    HBRUSH brush = CreateSolidBrush(RGB(0xff,0xff,0xff));
    GetCharABCWidths(hMemDC, 0, 255, abc); 
    GetTextMetrics(hMemDC, &tm);

    //WIDE_CHAR_T base = 32;
    //WIDE_CHAR_T last = 255;
    int base = 32;
    int last = 254;
    SetTextColor(hMemDC, RGB(0, 0, 0)) ;
    SetBkColor(hMemDC, RGB(0xff, 0xff, 0xff)) ;

    for(int i=base;i<=last;i++){
      int bm_width = (abc[i].abcA+abc[i].abcB+abc[i].abcC + 7) / 8;
      int bm_height = tm.tmHeight;
      HBITMAP hbmp=CreateCompatibleBitmap(hdc,8*bm_width,bm_height);
      HGDIOBJ obj2 = SelectObject(hMemDC, hbmp);
      RECT r = {0,0,8*bm_width,bm_height};
      FillRect(hMemDC, &r, brush) ;
      char text_out[2] = {i,'\0'};
      TextOut(hMemDC,0,0,text_out,strlen(text_out));
      unsigned char *pixel_data = GetPixelData(hbmp, hMemDC, bm_height);
      DeleteObject(hbmp);
      DeleteObject(obj2);
/*
 * abcA
 * Specifies the "A" spacing of the character. The "A" spacing is the 
 * distance to add to the current position before drawing the character glyph.
 * This is lbearing in XFont speak?
 *
 * abcB
 * Specifies the "B" spacing of the character. The "B" spacing is the 
 * width of the drawn portion of the character glyph.
 * This is width in XFont speak?
 *
 * abcC
 * Specifies the "C" spacing of the character. The "C" spacing is the 
 * distance to add to the current position to provide white space to 
 * the right of the character glyph.
 * This is what? in XFont speak?
 *
 */ 
      if (pixel_data) {
        bm.push_back(pixel_data);
        widths.push_back(bm_width*8);
        heights.push_back(tm.tmHeight);
        lbearings.push_back(-abc[i].abcA);
        descents.push_back(0);
        dxs.push_back(abc[i].abcB+abc[i].abcC);
        dys.push_back(0);
        //std::cout << "Loaded " <<  foundry << " " << family << std::endl;
      } else {
        std::cout << "Error reading font " << foundry << " " << family << std::endl;
      }
    }
    DeleteObject(brush);
    DeleteObject(obj);

    std::vector<std::string> name;
    name.push_back("default-foundry");
    name.push_back(family);
    name.push_back(weight);
    name.push_back(slant);
    name.push_back(IntToString(isize));
    name.push_back(encoding);
    finfo = MGFontInfo(name,bm,widths,heights,lbearings,descents,dxs,dys,base,last,FONT_BITMAP_GRAYSCALE);
  }
  return finfo;
}

MGFontInfo LoadWinFont(const std::string &foundry, const std::string &family,  const std::string &weight,  const std::string &slant,  const std::string &encoding, int isize){
  HDC hdc = GetDC(NULL) ;
  HDC hMemDC = CreateCompatibleDC(hdc);
  MGFontInfo finfo = LoadWinFont(foundry,family,weight,slant,encoding,isize,hdc,hMemDC);
  DeleteDC(hMemDC);
  DeleteDC(hdc);
  return finfo;
}

MGFontInfo LoadWinFont(const char *name){
  MGFontInfo finfo;
  std::vector<std::string> name_split = ParseX11Name(std::string(name));
  if(name_split.size()<6)
    return finfo;
  return LoadWinFont(name_split[0],name_split[1],name_split[2],name_split[3],name_split[4],atoi(name_split[5].c_str()));
}

int LoadAllWinFonts(){
  HDC hdc = GetDC(NULL) ;
  HDC hMemDC = CreateCompatibleDC(hdc);

  int aFontCount[] = { 0, 0, 0 }; 
 
  EnumFontFamilies(hdc, (LPCTSTR) NULL, 
        (FONTENUMPROC) EnumFamCallBack, (LPARAM) aFontCount);

  for(unsigned j=0;j<font_families.size();j++){
     //std::cout << "Loading " << font_families[j] << " medium r " << 12 << "\n";
     MGFontInfo finfo = LoadWinFont("", font_families[j],  "medium",  "r",  "", 22,hdc,hMemDC);
     if(finfo.isValid())
        FontCache::AddFont(finfo);
     /*
     std::cout << "Loading " << font_families[j] << " bold r " << 12 << "\n";
     finfo = LoadWinFont("", font_families[j],  "bold",  "r",  "", 12,hdc,hMemDC);
     if(finfo.isValid())
        FontCache::AddFont(finfo);
     std::cout << "Loading " << font_families[j] << " medium i " << 12 << "\n";
     finfo = LoadWinFont("", font_families[j],  "medium",  "i",  "", 12,hdc,hMemDC);
     if(finfo.isValid())
        FontCache::AddFont(finfo);
     std::cout << "Loading " << font_families[j] << " bold i " << 12 << "\n";
     finfo = LoadWinFont("", font_families[j],  "bold",  "i",  "", 12,hdc,hMemDC);
     if(finfo.isValid())
        FontCache::AddFont(finfo);
     */
  }
  DeleteDC(hMemDC);
  DeleteDC(hdc);

  return 0;
}

int FindNextWinFontSize(const MGFontInfo &finfo, int step){
  if(!finfo.isValid()) return 0;
  if((finfo.Size()+step)>0)
    return finfo.Size();
  return finfo.Size()+step;
}
#endif
