/*
 *  CXXSpace.cpp
 *  lpbSolver
 *
 *  Created by gruber on Thu Jul 15 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include "CXXSpace.h"

CXXSpace::CXXSpace() { // dummy constructor
}

CXXSpace::CXXSpace(float probeRadius, float gridSpacing, 
				   float xMin, float xMax, float yMin, float yMax, 
				   float zMin, float zMax):SolventMap(gridSpacing, probeRadius, 
													  xMin, xMax, yMin, yMax, zMin, zMax) {
					   
					   saltConc = 0.0;
					   temp = 300;
					   
					   // WARNING dummy value - if this is not set by defineBondaryConditions will raise and Exception...
					   dielectricBoundary = -1;
					   
					   // constructor calls base class constructor (solvent map constructor) generating fftw optimised conservative
					   // discrete grid - now need to make sure the chargeGrid uses the same dimensions and is consistent with the
					   // optimised solvnetMap grid.
					   try {				   
						   chargeGrid = new double[dim[0]*dim[1]*dim[2]];
						   solvationGrid = new double[dim[0]*dim[1]*dim[2]];
						   potentialGrid = new double[dim[0]*dim[1]*dim[2]];
						   dielGrid = new CXXCoord<CXXCoord_ftype>[dim[0]*dim[1]*dim[2]];
						   epsilonKappaSq = new double[dim[0]*dim[1]*dim[2]];
						   
						   if (chargeGrid == 0 | potentialGrid == 0 | dielGrid == 0) {
							   CXXException theException = CXXException(" ERROR: (CXXSpace::CXXSpace()) :Could not reserve suffiecent memory !\n");
							   throw theException;
						   }
						   
						   
						   // set all values in chargeGrid / potentialGrid and dielectric Grids to 0 		
						   for (int i = 0; i < dim[0]; i++) {
							   for (int j  = 0; j < dim[1]; j++) {
								   for (int k  = 0; k < dim[2]; k++) {
									   setChargeGrid(i, j, k, 0.);
									   setPotential(i,j,k,0.);
									   setDielGrid(i,j,k,0,0.);
									   setDielGrid(i,j,k,1,0.);
									   setDielGrid(i,j,k,2,0.);
									   setSolvationGrid(i,j,k,0);
								   }
							   }
						   }
						   
						   
						   
						   // If we got here no problems where identified so now give a simple status report about what is happening
						   std::cout << "Generated Space:	\nReal origin		x: " << originOfGrid[0] << " y: " << originOfGrid[1] << " z: " << originOfGrid[2]; 
						   std::cout << "					\nGrid dimensions   i: " << dim[0] << " j: " << dim[1] << " k: " << dim[2]; 
						   std::cout << "					\nGrid spacing: " << gridSpacing << " \nProbeRadius: " << probeRadius << "\n"; 
						   
					   }	
					   catch (CXXException theExeption) {
						   theExeption.Report();
						   throw theExeption;
					   }
				   }


int CXXSpace::setChargeGrid(int i, int j, int k, double value) {
	
	if((i + j*dim[0] + k*dim[0]*dim[1]) >= dim[0]*dim[1]*dim[2]){
		CXXException theException = CXXException("ERROR in: CXXSpace::setChargeGrid - index error");
		throw theException;
	}	
	chargeGrid[i+j*dim[0]+ k*dim[0]*dim[1]]= value;
	return 0;
}

int CXXSpace::setDielGrid(int i, int j, int k, int direction,  double value) {
	
	if((i + j*dim[0] + k*dim[0]*dim[1]) >= dim[0]*dim[1]*dim[2]){
		CXXException theException = CXXException("ERROR in: CXXSpace::setChargeGrid - index error");
		throw theException;
	}	
	switch (direction) {
		case 0: // looking into direction i
			dielGrid[i+j*dim[0]+ k*dim[0]*dim[1]].setX(value);
			break;
			
		case 1: // looking into direction j
			dielGrid[i+j*dim[0]+ k*dim[0]*dim[1]].setY(value);
			break;
			
		case 2: // looking into direction k
			dielGrid[i+j*dim[0]+ k*dim[0]*dim[1]].setZ(value);
			break;		
	}
	return 0;
}

int CXXSpace::setPotential(int i, int j, int k, double value) {
	
	if((i + j*dim[0] + k*dim[0]*dim[1]) >= dim[0]*dim[1]*dim[2]){
		
		CXXException theException = CXXException("ERROR in: CXXSpace::setPotential - index error");
		throw theException;
	}	
	potentialGrid[i+j*dim[0]+ k*dim[0]*dim[1]] = value;
	return 0;
}

int CXXSpace::setSolvationGrid(int i, int j, int k, double value) {
	
	if((i + j*dim[0] + k*dim[0]*dim[1]) >= dim[0]*dim[1]*dim[2]){
		
		CXXException theException = CXXException("ERROR in: CXXSpace::setPotential - index error");
		throw theException;
	}	
	solvationGrid[i+j*dim[0]+ k*dim[0]*dim[1]] = value;
	return 0;
}

int CXXSpace::setEpsilonKappaSq(int i, int j, int k, double value) {
	
	if((i + j*dim[0] + k*dim[0]*dim[1]) >= dim[0]*dim[1]*dim[2]){
		
		CXXException theException = CXXException("ERROR in: CXXSpace::setepsilonKappaSq - index error");
		throw theException;
	}	
	epsilonKappaSq[i+j*dim[0]+ k*dim[0]*dim[1]] = value;
	return 0;
}


int CXXSpace::addDielectricSphere(double x, double y, double z, double r) { //mmdbManager thePdb) {
	
	// adds an atom to the solvent map
	int returnInt = atomAt3d(  x,  y, z, r);
	return returnInt;
	
}


int CXXSpace::defineBoundaryConditions(double value) {
	
	dielectricBoundary = value;
	return 0;
}



int CXXSpace::addGridCharge(int i, int j, int k,  double gridCharge) { 
	
	if((i + j*dim[0] + k*dim[0]*dim[1]) >= dim[0]*dim[1]*dim[2]){
		CXXException theException = CXXException("ERROR in: CXXSpace::addGridCharge - index error");
		throw theException;
	}
	
	double charge = getGridCharge(i,j,k) + gridCharge;	
	setChargeGrid(i, j, k, charge);	
	return 0;
	
}

double CXXSpace::getGridCharge(int i, int j, int k) { 
	
	if (i < 0 | j < 0 | k < 0 | i == dim[0] | j == dim[1] | k == dim[2]) {
		return 0; // WARNING ugly  - need  proper boundary value ...
	}
	else {
	
	if((i + j*dim[0] + k*dim[0]*dim[1]) >= dim[0]*dim[1]*dim[2]){
		
		CXXException theException = CXXException("ERROR in: CXXSpace::getGridCharge - index error");
		throw theException;
		
	}	
	double charge = chargeGrid[i+j*dim[0]+ k*dim[0]*dim[1]];	
	return charge;
	}
}


double CXXSpace::getGridSolvationParameter(int i, int j, int k) { 
	
	if (i < 0 | j < 0 | k < 0 | i == dim[0] | j == dim[1] | k == dim[2]) {
		return 0; // WARNING ugly  - need  proper boundary value ...
	}
	else {
	
	if((i + j*dim[0] + k*dim[0]*dim[1]) >= dim[0]*dim[1]*dim[2]){
		
		CXXException theException = CXXException("ERROR in: CXXSpace::getGridCharge - index error");
		throw theException;
		
	}	
	double parameter = solvationGrid[i+j*dim[0]+ k*dim[0]*dim[1]];	
	return parameter;
	}
}


double CXXSpace::getEpsilonKappaSq(int i, int j, int k) { 
	
	if (i < 0 | j < 0 | k < 0 | i == dim[0] | j == dim[1] | k == dim[2]) {
		return epsilonKappaSq[0]; // WARNING ugly  - need  proper boundary value ...
	}
	else {
	if((i + j*dim[0] + k*dim[0]*dim[1]) >= dim[0]*dim[1]*dim[2]){
		
		CXXException theException = CXXException("ERROR in: CXXSpace::getEpsilonKappa - index error");
		throw theException;
		
	}	
		double charge = epsilonKappaSq[i+j*dim[0]+ k*dim[0]*dim[1]];	
	return charge;
	}
}


double CXXSpace::getPotential(int i, int j, int k) {	
	
	if (i < 0 | j < 0 | k < 0 | i == dim[0] | j == dim[1] | k == dim[2]) {
		return 0; // WARNING later have alternative boundary conditions ....
	}
	else {
		if((i + j*dim[0] + k*dim[0]*dim[1]) >= dim[0]*dim[1]*dim[2]){
			
			CXXException theException = CXXException("ERROR in: CXXSpace::getPotential - index error");
			throw theException;
			
		}	
		double value = potentialGrid[i+j*dim[0]+ k*dim[0]*dim[1]];	
		return value;
		
	}
}

double CXXSpace::getBoundaryMap(int i, int j, int k) { 
	
	// this special case only comes into effect when the dielectric grid is generated becuase this step 
	// involves averaging over box around i j k ± 1 - then should give boindary value ...
	
	if (i < 0 | j < 0 | k < 0 | i == dim[0] | j == dim[1] | k == dim[2]) {
		return dielectricBoundary;
	}
	else {
		
		// WARNIG two different mappings from ijk into array used in solventMap and CXX classes ...
		if((i*dim[1]*dim[2] + j*dim[2] + k) > dim[0]*dim[1]*dim[2]){
			CXXException theException = CXXException("ERROR in: CXXSpace::getBoundaryMap - index error");
			throw theException;
		}	
		if(dielectricBoundary == -1) {
			CXXException theException = CXXException("ERROR in: CXXSpace::getBoundaryMap - no boundary conditions set");
			throw theException;
		}	
		else {
			
			double inOut = getSolidAt3i(i,j,k);
			return inOut;
		}
	}
}

double CXXSpace::getDielGrid(int i, int j, int k, int direction) {
	
	if(i<0||j<0||k<0) return 78;
	if((i + j*dim[0] + k*dim[0]*dim[1]) >= dim[0]*dim[1]*dim[2]){
		
		CXXException theException = CXXException("ERROR in: CXXSpace::getGridCharge - index error");
		throw theException;
		return 78;
	}	
	else switch (direction) {
		case 0: // looking into direction i
			return dielGrid[i+j*dim[0]+ k*dim[0]*dim[1]].x();
			break;
			
		case 1: // looking into direction j
			return dielGrid[i+j*dim[0]+ k*dim[0]*dim[1]].y();
			break;
			
		case 2: // looking into direction k
			return dielGrid[i+j*dim[0]+ k*dim[0]*dim[1]].z();
			break;
		default:
			return 78;
	}
}

CXXCoord<CXXCoord_ftype>CXXSpace::getOrigin() {
	
	CXXCoord<CXXCoord_ftype>theOrigin;
	theOrigin.setX( originOfGrid[0]);
	theOrigin.setY( originOfGrid[1]);
	theOrigin.setZ( originOfGrid[2]);
	
	return theOrigin;
}

CXXCoord<CXXCoord_ftype>CXXSpace::getSpaceSpanningVectorX() {
	
	// for now space is always orthogonal and equidistant with grid spacing ...
	CXXCoord<CXXCoord_ftype>span(gridSpacing, 0., 0.);
	return span;
}

CXXCoord<CXXCoord_ftype>CXXSpace::getSpaceSpanningVectorY() {
	
	// for now space is always orthogonal and equidistant with grid spacing ...
	CXXCoord<CXXCoord_ftype>span(0.,gridSpacing, 0.);
	return span;
}

CXXCoord<CXXCoord_ftype>CXXSpace::getSpaceSpanningVectorZ() {
	
	// for now space is always orthogonal and equidistant with grid spacing ...
	CXXCoord<CXXCoord_ftype>span( 0., 0., gridSpacing);
	return span;
}

int CXXSpace::getDimI() { return dim[0];}
int CXXSpace::getDimJ() { return dim[1];}
int CXXSpace::getDimK() { return dim[2];}

int CXXSpace::setSolventParameters(double salt, double temperature) {
	
	saltConc = salt;
	temp = temperature;
	
	return 0;
}

int CXXSpace::introduceMedium(double dielectricInMedium, double dielectricInProtein, double probeRadius) { 
	
	try {
		
		// this calculates the envelope of the protein giving values in the interior a flag of zero in
		// solidMap while points with solvent contact are larger than zero...
		std::cout << "\nNow generating solvent envelope for protein\n";
		
		// Generating solvent map on identical grid as charge grid first 
		// this grid has values of zero inside protein and otherwise of 1
		
		// Smooth dielectric by method:
		// i)   assign all points currently 0 (inside) the dielectric value of the protein interiour
		// ii)	assign all points currently 1 (outside) the dielectric value of the solvent
		
		convoluteSolidProbe(probeRadius, 0,0,0);
		for (int i = 0; i < dim[0]*dim[1]*dim[2]; i++) {    
			if ( SolidGrid[i] > 0.01 )
				SolidGrid[i] = dielectricInMedium;
			else
				SolidGrid[i] = dielectricInProtein;
		} 
		
		// iii) loop over all grid points of the three dielectric grids and assign to each point in these
		//		three grids the the geometric average of the ten nearest solvent grid points
		std::cout << "Dielectric smoothing in progress\n";
		
		double dielLookI, dielLookJ, dielLookK;
		
//#pragma omp parallel for default(none) shared(dielLookI, dielLookJ, dielLookK)
		for (int i = 0; i < dim[0]; i++) {	
			for (int j  = 0; j < dim[1]; j++) {
				for (int k  = 0; k < dim[2]; k++) {
					
					// ten nearest points for i,j,k looking into i direction (that is looking i -> i+1):
					// i,j,k i,j+1,k i,j-1,k i,j,k-1 i,j,k+1 i+1,j,k i+1,j+1,k i+1,j-1,k i+1,j,k-1 i+1,j,k+1
					
					dielLookI = 10*1/( 1/getBoundaryMap(i,j,k)	+   1/getBoundaryMap(i+1,j,k) +
									   1/getBoundaryMap(i,j+1,k)  +   1/getBoundaryMap(i+1,j+1,k) +
									   1/getBoundaryMap(i,j,k+1)  +   1/getBoundaryMap(i+1,j,k+1) + 
									   1/getBoundaryMap(i,j-1,k)  +   1/getBoundaryMap(i+1,j-1,k) +
									   1/getBoundaryMap(i,j,k-1)  +   1/getBoundaryMap(i+1,j,k-1));
					
					setDielGrid(i,j,k,0,dielLookI);
					
					// now do the same for j  directions ...
					
					dielLookJ = 10*1/( 1/getBoundaryMap(i,j,k)	+   1/getBoundaryMap(i,j+1,k) +
									   1/getBoundaryMap(i+1,j,k)  +   1/getBoundaryMap(i+1,j+1,k) +
									   1/getBoundaryMap(i,j,k+1)  +   1/getBoundaryMap(i,j+1,k+1) + 
									   1/getBoundaryMap(i-1,j,k)  +   1/getBoundaryMap(i-1,j+1,k) +
									   1/getBoundaryMap(i,j,k-1)  +   1/getBoundaryMap(i,j+1,k-1));
					
					setDielGrid(i,j,k,1,dielLookJ);
					
					// finally do the same for k directions ...
					
					dielLookK = 10*1/( 1/getBoundaryMap(i,j,k)	+   1/getBoundaryMap(i,j,k+1) +
									   1/getBoundaryMap(i+1,j,k)  +   1/getBoundaryMap(i+1,j,k+1) +
									   1/getBoundaryMap(i,j+1,k)  +   1/getBoundaryMap(i,j+1,k+1) + 
									   1/getBoundaryMap(i-1,j,k)  +   1/getBoundaryMap(i-1,j,k+1) +
									   1/getBoundaryMap(i,j-1,k)  +   1/getBoundaryMap(i,j-1,k+1));
					
					
					setDielGrid(i,j,k,2,dielLookK);
				}
			}
		}
		
		// now use the dielecric grid to calculate the epsilonKappaSq
		// epsilonKappaSq stores for each point pint in the solvent the sum of all surrounding dielectric values 
		// multiplied by solventDielectric*kappaSpuared 
		// where:  kappa = Sqrt( 8*Pi*bulcSaltConc*e*e / (Avogadro*BoltzmannConst*Temperature))
		

		
		// kappa = sqrt(8*e*e*Na*I/(1000*k*T))
		// where I ionic strngth in M and T absolute temperature in K
		// constants = sqrt(8*e*e*Na/(1000*k))
		double constants = 5304.75324359;
		// => kappa = constants*sqrt(I/T)
		double kappa = constants*sqrt(saltConc/temp);
		// kappa factore in linear PBE k*k*h*h where h in m !
		double h = gridSpacing*(1e-10);
		double kappaFactor = kappa*kappa*h*h;
		
		double solventDielectric = 78;	// MAGIC NUBERS HERE

//#pragma omp parallel for default(none) shared(solventDielectric, kappaFactor)
		for (int i = 0; i < dim[0]; i++) {	
			for (int j  = 0; j < dim[1]; j++) {
				for (int k  = 0; k < dim[2]; k++) {
					double value =   getDielGrid(i-1,j,k,0) + getDielGrid(i,j,k,0)
							+ getDielGrid(i,j-1,k,1) + getDielGrid(i,j,k,1)
							+ getDielGrid(i,j,k-1,2) + getDielGrid(i,j,k,2);
					
					// only the points inside get values other then just the dielectrics ....
					if (getBoundaryMap(i,j,k) == solventDielectric) {
						value = value + kappaFactor;																
					}
					setEpsilonKappaSq(i,j,k,value);
				}
			}
		}
		
		
		
	}
	
	catch( CXXException theException) {
		theException.Report();
		throw theException;
		return 1;
	}
	// status lines
	std::cout << "Solvent envelope generatd.\nDielectric inside protein: " << dielectricInProtein << "\nDielectric in solvent:  " << dielectricInMedium << "\n" ;
	
	return 0;
}



int CXXSpace::dumpSpaceSlice(int dimension, int gridType, int sliceNumber){
	
	/*  dumps a slice through space perpendicular to <dimension> with property of <type>. Space is sliced 
	at position <sliceNumber> along <dimension> axis. The dump is a textfile called dump.dat:
	
dimension: 
	
	0	x
	1   y
	3   z
	
	looking at property:
	
type:
	
	1   charge
	2   dielectric value
	3   potential
	
	*/
	// first make sure that the slice exists in the dimension specified
	// make sure this direction is one of the three axis ..
	int i, j, k;
	//	fstream dump ("dumpfile",ios::out|ios::app );
	
	std::fstream dumper ("dumpfile.txt",ios::out|ios::app );
	
	dumper.width(5);
	dumper.precision(2);
	
	//dumpXSlice(1,sliceNumber, 3);
	
	
	if (dimension > 2 | dimension < 0) {
		
		CXXException theException("ERROR CXXSpace::dumpSpaceSlice() - dimension can only be 0,1 or 2!");
		throw theException;
	}
	
	
	// make sure have that many slices in this direction 
	switch (dimension) {
		case 0: // x-slice
			if (sliceNumber > dim[0]) {
				CXXException theException("ERROR CXXSpace::dumpSpaceSlice() 0 - outside of grid boundary ");
				throw theException;
			}
			// now dump x slice...
			dumper << "\nAtomGrid: Slicing through 3D grid at x= " << sliceNumber <<". Directions: (horizontal/vetical)  <-> (z/y): \n\n"; 
			for (k = 0; k < dim[2]; k++) {
				for (j = 0; j < dim[1]; j++) {					
					
					switch (gridType) {
						case 1: 		
							// dump charge grid							
							dumper << getGridCharge(sliceNumber, j, k) << " ";	
							
							if (j == (dim[1]-1))
								dumper << "\n\n";
								break;
						case 2:
							// dump boundaryMap grid
							
							dumper << setprecision(4) << getDielGrid(sliceNumber, j, k, 1) << " ";	
							if (j == (dim[1]-1))
								dumper << "\n\n";
								
								break;
						case 3:
							
							// dump potential - yet to come ...
							
							dumper << setprecision(4) << getPotential(sliceNumber, j, k) << " ";	
							if (j == (dim[1]-1))
								dumper << "\n\n";
								
							
							
							break;
						default:
							CXXException theException("ERROR CXXSpace::dumpSpaceSlice() - No such grid !");
							throw theException;
					}					
				}
			}	
				break;
			// y-dimension	
		case 1: // y-slice
			if (sliceNumber > dim[1]) {
				CXXException theException("ERROR CXXSpace::dumpSpaceSlice() 1 - outside of grid boundary ");
				throw theException;
			}
			dumper << "\nAtomGrid: Slicing through 3D grid at y= " << sliceNumber <<". Directions: (horizontal/vetical)  <-> (x/z): \n\n"; 
			// now dump y slice ...
			
			
			for (k = 0; k < dim[2]; k++) {
				for (i = 0; i < dim[0]; i++) {
					
					switch (gridType) {
						case 1: 		
							// dump charge grid
							
							dumper << setprecision(3) << getGridCharge(i, sliceNumber, k) << " ";	
							if (i == (dim[0]-1))
								dumper << "\n\n";
								
								break;
						case 2:
							// dump Boundary grid
							
							dumper << setprecision(3) << getBoundaryMap(i, sliceNumber, k) << " ";	
							if (i == (dim[0]-1))
								dumper << "\n\n";
								
								
								break;
						case 3:
							// dump potential - yet to come ...
							
							break;
						default:
							CXXException theException("ERROR CXXSpace::dumpSpaceSlice() - No such grid !");
							throw theException;
					}
				}
			}	
				
				break;
		case 2: // z-slice
			if (sliceNumber > dim[2]) {
				CXXException theException("ERROR CXXSpace::dumpSpaceSlice() 2 - outside of grid boundary ");
				throw theException;
			}
			
			// now dump z slice ...
			dumper << "\nAtomGrid: Slicing through 3D grid at z= " << sliceNumber <<". Directions: (horizontal/vetical)  <-> (x/y): \n\n"; 
			
			for (j = 0; j < dim[1]; j++) {
				for (i = 0; i < dim[0]; i++) {
					
					switch (gridType) {
						case 1: 		
							// dump charge grid
							
							dumper << setprecision(3) << getGridCharge(i, j, sliceNumber) << " ";	
							if (i == (dim[0]-1))
								dumper << "\n\n";
								
								break;
						case 2:
							// dump dielectric grid
							dumper << "BoundaryMap dump\n\n";	
							dumper << setprecision(3) << getBoundaryMap(i, j, sliceNumber) << " ";	
							if (i == (dim[0]-1))
								dumper << "\n\n";
								
								break;
						case 3:
							// dump potential - yet to come ...
							
							break;
						default:
							CXXException theException("ERROR CXXSpace::dumpSpaceSlice() - No such grid !");
							throw theException;
					}
				}
			}
				break;
			// end of outer switch clause - dimensiom
	}
	
	return 0;
	
}




