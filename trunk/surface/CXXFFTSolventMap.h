#include "rfftw.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <stdio.h>

#ifndef  __CXXException__
#include "CXXException.h"
#endif

class SolventMap {

 protected:
  
  int dim[3];
  fftw_real *Grid;
  fftw_real *distanceGrid;
  fftw_real *SolidGrid;
  fftw_complex *FFTGrid;  
  fftw_real probeRadius, gridSpacing;
  fftw_real angstroemMinMax[6];
  fftw_real originOfGrid[3];
  int intOriginOfGrid[3];
  fftwnd_plan realToComplexPlan, complexToRealPlan;
  int countRapid;
  void positionGrid ();
  void findGridDim();
  void optimizeGridDim();

  void prepareForRapidFFT();
  void cleanUpForRapidFFT();  

 public:
	  SolventMap();
   void makeDualContactMap();
   SolventMap (fftw_real dGrid, fftw_real rProbe, fftw_real xmin, fftw_real xmax, 
                fftw_real ymin, fftw_real ymax, fftw_real zmin, fftw_real zmax);
  ~SolventMap ();
  int atomAt3d( fftw_real x, fftw_real y, fftw_real z, fftw_real AtomRadius);
  int atomAt3dv( fftw_real xyz[3], fftw_real AtomRadius);
  int atomAt3f( float x, float y, float z, float AtomRadius);
  int atomAt3fv( float xyz[3], float AtomRadius);
  int convoluteSolidProbe (fftw_real choosenProbeRadius, int rapidFlag, int smoothFlag, float smoothRadius );
  void makeDistMap (int sampleNr); 
  fftw_real getSpacing(); 
  fftw_real getProbeRadius(); 
  int getFloatOrigin(fftw_real *theOrigin);
  int getIntOrigin(int *theIntOrigin);
  int getExtent(int *theExtent);
  void setDistanceGrid3i (int i, int j, int k, float value);
  void setDistanceGrid3iv (int *ijk, float value);
  void setGrid3i (int i, int j, int k, float value);
  void setGrid3iv (int *ijk, float value);
  void setSolidGrid3i (int i, int j, int k, float value);
  void setSolidGrid3iv (int *ijk, float value);
  fftw_real getDistanceAt3i (int i, int j, int k);
  fftw_real getDistanceAt3iv (int *ijk);
  fftw_real getSolidAt3i (int i, int j, int k);
  fftw_real getSolidAt3iv (int *ijk);
  int dumpXSlice(int type, int SliceNr, int broad);
};



























