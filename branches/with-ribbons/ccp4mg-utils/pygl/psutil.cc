/*
     pygl/psutil.cc: CCP4MG Molecular Graphics Program
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


#include <psutil.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#if defined (_WIN32) && !defined (__CYGWIN__) 
#include <windows.h>
#else
#include <sys/time.h>
#include <unistd.h>
#include <pwd.h>
#include <sys/types.h>
#endif

void PSTail(std::ofstream &fp){
  fp << "$F2psEnd\nrs\nshowpage\n";
}

void PSHeader(std::ofstream &fp, const std::string &filename){
  PSHeader(fp,filename,"");
}

void PSHeader(std::ofstream &fp, const std::string &filename, const std::string &papersize){

  std::string mypapersize = "A4";

  mypapersize = papersize;

  if(mypapersize=="")
    mypapersize = "A4";

  std::string bounding_box = "";
  std::string page_size_feature = "";

  if(mypapersize.substr(0,12)==std::string("BoundingBox:")){
    bounding_box = std::string("%%") + mypapersize + std::string("\n");
  } else {
    page_size_feature = std::string("%%IncludeFeature: *PageSize ") + mypapersize + std::string("\n");
  }

  std::string paper_scaling_string = "";

  if(mypapersize=="a0"||mypapersize == "A0"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n1.122749 0.000000 translate\n4.002843 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 595.560976 0 rlineto 0 842.000000 rlineto -595.560976 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="a1"||mypapersize == "A1"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n1.001552 0.000000 translate\n2.831279 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 595.707491 0 rlineto 0 842.000000 rlineto -595.707491 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="a2"||mypapersize == "A2"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n0.769260 -0.000000 translate\n2.003105 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 595.768067 0 rlineto 0 842.000000 rlineto -595.768067 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="a3"||mypapersize == "A3"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n0.189430 0.000000 translate\n1.414489 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 595.267842 0 rlineto 0 842.000000 rlineto -595.267842 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="a5"||mypapersize == "A5"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n0.271378 0.000000 translate\n0.706651 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n595.768067 0 rlineto 0 842.000000 rlineto -595.768067 0 rlineto\nclosepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="a6"||mypapersize == "A6"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n419.527559 0.589148 translate\n90 rotate\n0.498251 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 597.364865 0 rlineto 0 842.000000 rlineto -597.364865 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="a7"||mypapersize == "A7"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n297.239992 0.000000 translate\n90 rotate\n0.352544 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 595.000000 0 rlineto 0 844.256757 rlineto -595.000000 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="a8"||mypapersize == "A8"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n209.177794 0.000000 translate\n90 rotate\n0.247734 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 595.000000 0 rlineto 0 846.730769 rlineto -595.000000 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="a9"||mypapersize == "A9"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n147.401575 0.360222 translate\n90 rotate\n0.175061 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 599.115385 0 rlineto 0 842.000000 rlineto -599.115385 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="a10"||mypapersize == "A10"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n104.588897 0.000000 translate\n90 rotate\n0.123867 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 595.000000 0 rlineto 0 846.730769 rlineto -595.000000 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="2a0"||mypapersize=="2a"||mypapersize=="2A"||mypapersize == "2A0"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n0.585782 0.000000 translate\n5.662558 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 595.206897 0 rlineto 0 842.000000 rlineto -595.206897 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="4a0"||mypapersize=="4a"||mypapersize=="4A"|| mypapersize == "4A0"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n2.245497 0.000000 translate\n8.005686 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto 595.560976 0 rlineto 0 842.000000 rlineto -595.560976 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="b0"||mypapersize == "B0"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n1.127798 0.000000 translate\n4.760319 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 595.473833 0 rlineto 0 842.000000 rlineto -595.473833 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="b1"||mypapersize == "B1"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n0.494885 0.000000 translate\n3.366563 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 595.294000 0 rlineto 0 842.000000 rlineto -595.294000 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="b2"||mypapersize == "B2"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n0.563899 0.000000 translate\n2.380160 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 595.473833 0 rlineto 0 842.000000 rlineto -595.473833 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="b3"||mypapersize == "B3"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n0.000000 0.652683 translate\n1.681731 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 595.000000 0 rlineto 0 842.776204 rlineto -595.000000 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="b4"||mypapersize == "B4"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n0.782726 0.000000 translate\n1.188397 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 596.317280 0 rlineto 0 842.000000 rlineto -596.317280 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="b5"||mypapersize == "B5"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n0.000000 1.329187 translate\n0.838483 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 595.000000 0 rlineto 0 845.170455 rlineto -595.000000 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="b6"||mypapersize == "B6"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n0.892139 -0.000000 translate\n0.592515 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 598.011364 0 rlineto 0 842.000000 rlineto -598.011364 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="b7"||mypapersize == "B7"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n0.000000 0.664593 translate\n0.419242 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 595.000000 0 rlineto 0 845.170455 rlineto -595.000000 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="b8"||mypapersize == "B8"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n0.000000 0.371601 translate\n0.295375 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 595.000000 0 rlineto 0 844.516129 rlineto -595.000000 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="b9"||mypapersize == "B9"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n0.265958 0.000000 translate\n0.208727 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 597.548387 0 rlineto 0 842.000000 rlineto -597.548387 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="b10"||mypapersize == "B10"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n0.000000 0.185800 translate\n0.147687 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 595.000000 0 rlineto 0 844.516129 rlineto -595.000000 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="2b0"||mypapersize == "2B0"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n0.989769 0.000000 translate\n6.733125 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 595.294000 0 rlineto 0 842.000000 rlineto -595.294000 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="4b0"||mypapersize == "4B0"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n2.255597 0.000000 translate\n9.520639 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n 595.473833 0 rlineto 0 842.000000 rlineto -595.473833 0 rlineto\n closepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="letter"||mypapersize == "Letter"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n26.166271 0.000000 translate\n0.940618 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n650.636364 0 rlineto 0 842.000000 rlineto -650.636364 0 rlineto\nclosepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="legal"||mypapersize == "Legal"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n0.000000 70.971429 translate\n1.028571 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n595.000000 0 rlineto 0 980.000000 rlineto -595.000000 0 rlineto\nclosepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="tabloid"||mypapersize == "Tabloid"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n0.000000 51.610084 translate\n1.331092 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n595.000000 0 rlineto 0 919.545455 rlineto -595.000000 0 rlineto\nclosepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="statement"||mypapersize == "Statement"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n0.000000 25.805042 translate\n0.665546 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n595.000000 0 rlineto 0 919.545455 rlineto -595.000000 0 rlineto\nclosepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="executive"||mypapersize == "Executive"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n15.605701 0.000000 translate\n0.855107 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n631.500000 0 rlineto 0 842.000000 rlineto -631.500000 0 rlineto\nclosepath}put initclip\nPStoPSxform concat\n");
  }else if(mypapersize=="folio"||mypapersize == "Folio"){
    paper_scaling_string = std::string("userdict/PStoPSsaved save put\nPStoPSmatrix setmatrix\n0.000000 34.971429 translate\n1.028571 dup scale\nuserdict/PStoPSmatrix matrix currentmatrix put\nuserdict/PStoPSclip{0 0 moveto\n595.000000 0 rlineto 0 910.000000 rlineto -595.000000 0 rlineto\nclosepath}put initclip\nPStoPSxform concat\n");
  }


 std::string paper_scaling_defs  = 
  std::string("/PStoPSmatrix matrix currentmatrix def\n/PStoPSxform matrix def/PStoPSclip{clippath}def\n/defaultmatrix{PStoPSmatrix exch PStoPSxform exch concatmatrix}bind def\n/initmatrix{matrix defaultmatrix setmatrix}bind def\n/initclip[{matrix currentmatrix PStoPSmatrix setmatrix\n[{currentpoint}stopped{$error/newerror false put{newpath}}\n{/newpath cvx 3 1 roll/moveto cvx 4 array astore cvx}ifelse]\n{[/newpath cvx{/moveto cvx}{/lineto cvx}\n{/curveto cvx}{/closepath cvx}pathforall]cvx exch pop}\nstopped{$error/errorname get/invalidaccess eq{cleartomark\n$error/newerror false put cvx exec}{stop}ifelse}if}bind aload pop\n/initclip dup load dup type dup/operatortype eq{pop exch pop}\n{dup/arraytype eq exch/packedarraytype eq or\n{dup xcheck{exch pop aload pop}{pop cvx}ifelse}\n{pop cvx}ifelse}ifelse\n{newpath PStoPSclip clip newpath exec setmatrix} bind aload pop]cvx def\n/initgraphics{initmatrix newpath initclip 1 setlinewidth\n0 setlinecap 0 setlinejoin []0 setdash 0 setgray\n10 setmiterlimit}bind def\n");

  std::string psheader = std::string("%%Orientation: Portrait\n%%Pages: 1\n") + 
bounding_box +
std::string("%%BeginSetup\n") +
page_size_feature + std::string("%%EndSetup\n%%Magnification: 1.0000\n%%EndComments\n/$F2psDict 200 dict def\n$F2psDict begin\n$F2psDict /mtrx matrix put\n/col-1 {0 setgray} bind def\n/col0 {0.000 0.000 0.000 srgb} bind def\n/col1 {0.000 0.000 1.000 srgb} bind def\n/col2 {0.000 1.000 0.000 srgb} bind def\n/col3 {0.000 1.000 1.000 srgb} bind def\n/col4 {1.000 0.000 0.000 srgb} bind def\n/col5 {1.000 0.000 1.000 srgb} bind def\n/col6 {1.000 1.000 0.000 srgb} bind def\n/col7 {1.000 1.000 1.000 srgb} bind def\n/col8 {0.000 0.000 0.560 srgb} bind def\n/col9 {0.000 0.000 0.690 srgb} bind def\n/col10 {0.000 0.000 0.820 srgb} bind def\n/col11 {0.530 0.810 1.000 srgb} bind def\n/col12 {0.000 0.560 0.000 srgb} bind def\n/col13 {0.000 0.690 0.000 srgb} bind def\n/col14 {0.000 0.820 0.000 srgb} bind def\n/col15 {0.000 0.560 0.560 srgb} bind def\n/col16 {0.000 0.690 0.690 srgb} bind def\n/col17 {0.000 0.820 0.820 srgb} bind def\n/col18 {0.560 0.000 0.000 srgb} bind def\n/col19 {0.690 0.000 0.000 srgb} bind def\n/col20 {0.820 0.000 0.000 srgb} bind def\n/col21 {0.560 0.000 0.560 srgb} bind def\n/col22 {0.690 0.000 0.690 srgb} bind def\n/col23 {0.820 0.000 0.820 srgb} bind def\n/col24 {0.500 0.190 0.000 srgb} bind def\n/col25 {0.630 0.250 0.000 srgb} bind def\n/col26 {0.750 0.380 0.000 srgb} bind def\n/col27 {1.000 0.500 0.500 srgb} bind def\n/col28 {1.000 0.630 0.630 srgb} bind def\n/col29 {1.000 0.750 0.750 srgb} bind def\n/col30 {1.000 0.880 0.880 srgb} bind def\n/col31 {1.000 0.840 0.000 srgb} bind def\n\nend\n") + paper_scaling_defs + std::string("save\n\n/cp {closepath} bind def\n/ef {eofill} bind def\n/gr {grestore} bind def\n/gs {gsave} bind def\n/sa {save} bind def\n/rs {restore} bind def\n/l {lineto} bind def\n/m {moveto} bind def\n/rm {rmoveto} bind def\n/n {newpath} bind def\n/s {stroke} bind def\n/sh {show} bind def\n/slc {setlinecap} bind def\n/slj {setlinejoin} bind def\n/slw {setlinewidth} bind def\n/srgb {setrgbcolor} bind def\n/rot {rotate} bind def\n/sc {scale} bind def\n/sd {setdash} bind def\n/ff {findfont} bind def\n/sf {setfont} bind def\n/scf {scalefont} bind def\n/sw {stringwidth} bind def\n/tr {translate} bind def\n/tnt {dup dup currentrgbcolor\n  4 -2 roll dup 1 exch sub 3 -1 roll mul add\n  4 -2 roll dup 1 exch sub 3 -1 roll mul add\n  4 -2 roll dup 1 exch sub 3 -1 roll mul add srgb}\n  bind def\n/shd {dup dup currentrgbcolor 4 -2 roll mul 4 -2 roll mul\n  4 -2 roll mul srgb} bind def\n/$F2psBegin {$F2psDict begin /$F2psEnteredState save def} def\n/$F2psEnd {$F2psEnteredState restore end} def\n%%Page: 1 1\n") + paper_scaling_string + std::string("10 setmiterlimit\n 1.00000 1.00000 sc\n%\n% CCP4MG Lines Begin\n%\n0.100 slw\n$F2psBegin\n");
  fp << "%!PS-Adobe-2.0\n"; fp << "%%Title: " << filename << "\n"; fp << "%%Creator: CCP4MG Version 1.0\n"; 
#if defined (_WIN32) && !defined (__CYGWIN__) 
#define LENLEN 128
  LPSYSTEMTIME lpSystemTime = new SYSTEMTIME();
  GetLocalTime(lpSystemTime);
  
  fp << "%%CreationDate: " << lpSystemTime->wHour << ":" << lpSystemTime->wMinute << ":" << lpSystemTime->wSecond << " " << lpSystemTime->wDay << "/" << lpSystemTime->wMonth << "/" << lpSystemTime->wYear << "\n";
  DWORD len;
  char userName[LENLEN];
  GetUserName(userName,&len);

  fp << "%%For: " << userName << "\n";

#else
  struct timeval tp;
  struct timezone tzp;
  std::string date;
  time_t t1;
  
  gettimeofday(&tp,&tzp);
  t1 = tp.tv_sec;
  date = ctime(&t1);

  char name[100];
  gethostname(name,100);
  struct passwd *pwd = getpwuid(getuid());

  fp << "%%CreationDate: " << date;
  fp << "%%For:" << pwd->pw_name << "@" << name << "(" << pwd->pw_gecos << ")\n";
#endif
  fp << psheader;

}
 
