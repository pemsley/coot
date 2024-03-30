/*
 *  creator.h
 *  lpbSolver
 *
 *  Created by gruber on Fri Jul 02 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <iomanip>


#ifndef  __STRING_H
#include <string>
#endif

#ifndef  __MMDB_Manager__
#include "mmdb2/mmdb_manager.h"
#include "mmdb2/mmdb_tables.h"
#endif
#ifndef  __CXXException__
#include "CXXException.h"
#endif

#include "CXXCoord.h"
#ifndef  __CXXException__
#include <CXXException.h>
#endif
#include "CXXChargeTable.h"
class CXXSpace;

#include "clipper/clipper.h"

using namespace std;

class CXXCreator {

private:
	
	// parameters - values set in constructor
	
	double gridSpacing; // in A
	double probeRadius; // in A
	double zeroSpace;
	double temp; // in K
	double ionicStrength; // in mol / dm3
	
	int antiAliasPara;
	double scale ;
	
	
	double externalDielectric; 
	double internalDielectric;
	
	// Knowledge about charges of residues
	CXXChargeTable theChargeTable;

	// Matter
	mmdb::Manager*  theMMDBManager;
	int  nSelAtoms;
	mmdb::Atom** SelAtom;
	
	// private read accessors
	double lookUpCharge(int atomNr);
	CXXCoord<CXXCoord_ftype>getAtomCoord(int atomNr);
	double getAtomRadius(int atonNr);
	string getAtomElement(int atomNr);
	string getAtomName(int atomNr);
	string getAtomResidueName(int atomNr);
	
		// private write accessors
	int addAtomToSpace();
	
	// non public stuff thats needs to be done
	double fractionOfGridPointCoveredByAtom(CXXCoord<CXXCoord_ftype>thePoint, CXXCoord<CXXCoord_ftype>xGridVector,CXXCoord<CXXCoord_ftype>yGridVector, CXXCoord<CXXCoord_ftype>zGridVector, int atomNr);
	double getGridVolumeOfAtom(CXXCoord<CXXCoord_ftype>gridOrigin, CXXCoord<CXXCoord_ftype>xGridVector, CXXCoord<CXXCoord_ftype>yGridVector, CXXCoord<CXXCoord_ftype>zGridVector, int atomNr);
	int distributeAtomCharge(CXXCoord<CXXCoord_ftype>gridOrigin, CXXCoord<CXXCoord_ftype>xGridVector, CXXCoord<CXXCoord_ftype>yGridVector, CXXCoord<CXXCoord_ftype>zGridVector, int atomNr);
	void init(); //Initialise various parameters, the charge table etc
	int selectAllAtoms(); //Minimal routine to return selection handle 
						  //of a selection containing all atoms
public:
	
	CXXCreator (mmdb::pstr thePdb);   // reads a pdb from file to make creator 
	CXXCreator (mmdb::Manager* theManager); //can also be created from an MMDBManager ...
	CXXCreator (mmdb::Manager* theManager, int selHnd, int context_selHnd=-1); //can also be created from an MMDBManager + selHnd

	int setParameters( double IonicStrength, double Temperature, double gridSpacing);
	
	int createSpace();  // calculates box size, makes space (reserves memory), init for FFTW / solventMap stuff
	
	// adds all the atom shperes to the dielectric map and makes solvent envelope (solventMap stuff)
	// then generates dielectric grids (i, j and k directed) then assignes and smoothes the dielectric maps
	int introduceMatter(double dielectricOfSolvent, double dielectricOfProtein);
	
	// Finds selfconsistent solution for the electrostatic potential based on what we have done so far ..
	int evolve(int optionSOR, double convergenceCriterion);

	// space is discrete ...
	CXXSpace *space;
	
	//Return the grid spacing
	double getGridSpacing();
	
	// force a calculation from defined selection with default parameters
	int calculate();

	// after calculation, coerce thecreators data into a clipper map for contouring, interpolating etc
	clipper::NXmap<double> coerceToClipperMap(clipper::Cell &cell);

};
	
