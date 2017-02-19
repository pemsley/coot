/* coot-surface/coot-surface.hh
 * 
 * Copyright 2005 The University of Oxford
 * Author: Martin Noble, Jan Gruber
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#include "clipper/core/nxmap.h"
#include "surface/CXXSurface.h"

namespace coot {

class CColourScheme {
  friend class CMolColour;
 public:
  CColourScheme();
  ~CColourScheme();
  int SetSchemeInt (  const std::vector<int>& ityp , const std::vector<std::string>& cols);
  int SetSchemeInt ( int n, int typ[] , char *cols[] );
  int SetSchemeFloat (  const std::vector<float>& rtyp , const std::vector<std::string>& cols);
  int SetSchemeFloat ( int n, mmdb::realtype typ[] , char *cols[] );
  int SetSchemeString (  const std::vector<std::string>& chtyp , const std::vector<std::string>& cols);
  int SetSchemeString ( int n, char *typ[] , char *cols[] );
  void SetBlendMode (int modein ) { blend_mode = modein; }

  std::vector<int> GetSchemeInt ( );
  std::vector<float> GetSchemeFloat ( );
  std::vector<std::string>GetSchemeString ( );
  std::vector<std::string> GetSchemeColours ( );
  std::vector<int> GetSchemeCodes ( );
  void Print(void);
  double GetFColour(double value);
  std::vector<double> GetRGB(double value);
  std::string GetMode() { return mode; }

 protected:
  int defColour;
  int blend_mode;
  int nTypes;
  std::vector<int> itypes;
  std::vector<float> ranges;
  std::vector<std::string> strtypes;
  std::vector<std::string> colours;
  std::vector<int> codes;
  std::string mode;

};

   class surface {

      CXX_mot::CXXSurface *theSurface;
      int iEval;
      
   public:

      void fill_from(mmdb::Manager *mol, int selHnd, float col_scale, bool assign_charges);
      void fill_surface(mmdb::Manager *mol, int selHnd_selection, int SelHnd_all, float col_scale,
			bool assign_charges);
      void draw(double *override_colour, int selective_override);
      void transparent_draw(float opacity);
      void evaluateElectrostaticPotential(mmdb::Manager *theManager, int selHnd, float col_scale);
      int evaluatePhiAndColourWithDefaultScheme(mmdb::Manager *theManager, const int selHnd, float col_scale);
      int evaluatePhiAndColourWithScheme(mmdb::Manager *theManager,
					 const int selHnd,
					 CColourScheme &colourScheme);
      int interpolateIntoMap(const std::string &coordinateType,
			     const std::string &scalarType,
			     const clipper::NXmap<double> &aMap);
      int colourByScalarValue(const std::string &scalarType,
			      CColourScheme &colourScheme);
   };
}
