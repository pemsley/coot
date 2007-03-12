
#include "clipper/core/nxmap.h"
#include "CXXSurface.h"

namespace coot {

class CColourScheme {
  friend class CMolColour;
 public:
  CColourScheme();
  ~CColourScheme();
  int SetSchemeInt (  const std::vector<int>& ityp , const std::vector<std::string>& cols);
  int SetSchemeInt ( int n, int typ[] , char *cols[] );
  int SetSchemeFloat (  const std::vector<float>& rtyp , const std::vector<std::string>& cols);
  int SetSchemeFloat ( int n, realtype typ[] , char *cols[] );
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

      CXXSurface *theSurface;
      int iEval;
      
   public:

      void fill_from(CMMDBManager *mol, int selHnd);
      void draw(double *override_colour, int selective_override);
      void evaluateElectrostaticPotential(CMMDBManager *theManager, int selHnd);
      int evaluatePhiAndColourWithDefaultScheme(CMMDBManager *theManager, const int selHnd);
      int evaluatePhiAndColourWithScheme(CMMDBManager *theManager,
					 const int selHnd,
					 CColourScheme &colourScheme);
      int interpolateIntoMap(const std::string &coordinateType,
			     const std::string &scalarType,
			     const clipper::NXmap<double> &aMap);
      int colourByScalarValue(const std::string &scalarType,
			      CColourScheme &colourScheme);
   };
}
