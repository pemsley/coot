/*
 * MoleculesToTriangles/CXXSurface/CXXSpace.h
 *
 * Copyright 2009 by Martin Noble, University of Oxford
 * Author: Martin Noble
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
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

#include "CXXCoord.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include "CXXFFTSolventMap.h"

#ifndef  __MMDB_Manager__
#include "mmdb2/mmdb_manager.h"
#include "mmdb2/mmdb_tables.h"
#endif

#ifndef  __CXXException__
#include <CXXException.h>
#endif

class CXXSpace:public SolventMap {
	
private:
	
	
	double *chargeGrid;		// charge density grid
	double *solvationGrid;  // empirical atomic solvation parameters
	double *potentialGrid;  // electrostatic potential grid 
	CXXCoord<CXXCoord_ftype>*dielGrid;		// dielectric constant between i,j,k and i+1,j,k charge grid points in x()
							// dielectric constant between i,j,k and i,j+1,k charge grid points in y()
							// dielectric constant between i,j,k and i,j,k+1 charge grid points in z()
							// of coord ...
	double *epsilonKappaSq;		// grid containg both the dielectric around a point and the debey-hueckel parmeter
	
	double dielectricBoundary; // constat for dielectric outside box ...
	double temp;
	double saltConc;
	
	int setSolvationGrid(int i, int j, int k, double value);
	int setChargeGrid(int i, int j, int k, double value);
	int setDielGrid(int i, int j, int k, int direction,  double value);
	int setEpsilonKappaSq(int i, int j, int k, double value);
	
public:
	
	CXXSpace() ;
	CXXSpace(float probeRadius, float gridSpacing,
			 float xMin, float xMax, float yMin, float yMax, float zMin, float zMax);

	int setSolventParameters(double saltConc, double temp);
	int addDielectricSphere(double x, double y, double z, double r);
	int addGridCharge(int i, int j, int k, double charge);
	int setPotential(int i, int j, int k, double value);
	
	
	CXXCoord<CXXCoord_ftype>getOrigin();
	CXXCoord<CXXCoord_ftype>getSpaceSpanningVectorX();
	CXXCoord<CXXCoord_ftype>getSpaceSpanningVectorZ();
	CXXCoord<CXXCoord_ftype>getSpaceSpanningVectorY();
	
	int getDimI();
	int getDimJ();
	int getDimK();
	
	double getGridCharge(int i, int j, int k);
	double getGridSolvationParameter(int i, int j, int k);
	double getBoundaryMap(int i, int j, int k);
	double getDielGrid(int i, int j, int k, int directionn);	
	double getPotential(int i, int j, int k);	
	double getEpsilonKappaSq(int i, int j, int k);	
	
	int defineBoundaryConditions(double dielectricValue);
	int introduceMedium(double dielctricInMedium, double dielectricInProtein, double probeRadius);
	int dumpSpaceSlice(int dimension, int gridType, int sliceNumber);
	
	
	
};
