/* 
 * 
 * Copyright 2004 by The University of Oxford
 * Author: Martin Noble, Jan Gruber
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */
/*
 *  creator.h
 *  lpbSolver
 *
 *  Created by gruber on Fri Jul 02 2004.
 *  
 *
 */

#include <iostream>
#include <iomanip>


#ifndef  __STRING_H
#include <string>
#endif

#ifndef  __MMDB_Manager__
#include "mmdb_manager.h"
#include "mmdb_tables.h"
#endif
#ifndef  __CXXException__
#include <CXXException.h>
#endif

#include <CXXCoord.h>
#ifndef  __CXXException__
#include <CXXException.h>
#endif
#include <CXXChargeTable.h>
#include <CXXSpace.h>

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/core/nxmap.h>

using namespace clipper;


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
	PCMMDBManager  theMMDBManager;
	int  nSelAtoms;
	PPCAtom SelAtom;
	
	// private read accessors
	double lookUpCharge(int atomNr);
	CXXCoord getAtomCoord(int atomNr);
	double getAtomRadius(int atonNr);
	string getAtomElement(int atomNr);
	string getAtomName(int atomNr);
	string getAtomResidueName(int atomNr);
	int addAtomSolvationParameters(); // WARNING - need this to check if the solvent map does yet excist...
	
		// private write accessors
	int addAtomToSpace();
	
	// non public stuff thats needs to be done
	double fractionOfGridPointCoveredByAtom(CXXCoord thePoint, CXXCoord xGridVector,CXXCoord yGridVector, CXXCoord zGridVector, int atomNr);
	double getGridVolumeOfAtom(CXXCoord gridOrigin, CXXCoord xGridVector, CXXCoord yGridVector, CXXCoord zGridVector, int atomNr);
	int distributeAtomCharge(CXXCoord gridOrigin, CXXCoord xGridVector, CXXCoord yGridVector, CXXCoord zGridVector, int atomNr);
	void init(); //Initialise various parameters, the charge table etc
	int selectAllAtoms(); //Minimal routine to return selection handle 
						  //of a selection containing all atoms
public:
	
	CXXCreator (pstr thePdb);   // reads a pdb from file to make creator 
	CXXCreator (PCMMDBManager theManager); //can also be created from an MMDBManager ...
	CXXCreator (PCMMDBManager theManager, int selHnd); //can also be created from an MMDBManager + selHnd

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
	
