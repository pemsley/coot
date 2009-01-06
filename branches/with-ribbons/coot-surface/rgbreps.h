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
//   by writing to the Free Software Foundation, Inc., 59 Temple Place,
//   Suite 330, Boston, MA 02111-1307 USA


#ifndef _CCP4MG_RGBREPS_
#define _CCP4MG_RGBREPS_


#include <map>
#include <vector>
#include <string>
#include <fstream>

enum ColourRepresentation { NOCOLOURREP , RGBCOLOURREP , HSVCOLOURREP };

class RGBReps {
  static std::vector<std::string> colours_name_map;
  static std::map<std::string,std::vector<double> > string_colours;
  static std::vector<std::vector<double> > colours;
  static std::vector<std::vector<double> >hsv;

  static void BuildRGBReps(std::ifstream &infile);
 public:
  RGBReps(std::ifstream &infile);
  RGBReps(std::string &infile);
  RGBReps(const char *infile);
  RGBReps(){};
  ~RGBReps(){};

  static void SetFile(std::ifstream &infile);
  static void SetFile(std::string &infile);
  static void SetFile(const char *infile);

  static std::vector<std::vector<double> > GetColours();
  static std::vector<double> GetColour(unsigned int i);
  static double* GetColourP(unsigned int i);
  static std::vector<double> GetColourHSV(unsigned int i);
  static std::vector<std::vector<double> > GetColourDefinitions( std::vector<int> codes,int mode);
  static std::string GetNameOfColour(unsigned int i);
  static std::vector<double> GetColour(const std::string &name);
  static unsigned int GetNumberOfColours();
  static void AddColour(std::string colour, double r, double g, double b, double a=255.0);
  static void ChangeColour(std::string colour, double r, double g, double b, double a=255.0);
  static void DeleteColour(std::string colour);
  static int GetColourNumber(std::string name ); 
  static std::vector<int> GetColourNumbers(std::vector<std::string> names ); 
  static std::vector<double> hsvtorgb (std::vector<double> hsv);
  static std::vector<double> rgbtohsv (std::vector<double> rgb );
  static double* hsvtorgb (double *hsv);
  static double* rgbtohsv (double *rgb );

};

#endif
