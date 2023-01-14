/*
 *  creator.cpp
 *  lpbSolver
 *
 *  Created by gruber on Fri Jul 02 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include "CXXCreator.h"
#include "CXXSpace.h"
#include <algorithm>

#if defined _OPENMP
#include <omp.h>
#else
#if __APPLE__
#include <dispatch/dispatch.h>
#endif
#endif

CXXCreator::CXXCreator (mmdb::pstr thePdb) {
    mmdb::InitMatType();
    
    int RC;
    
    theMMDBManager = new mmdb::Manager();
    theMMDBManager->SetFlag( mmdb::MMDBF_PrintCIFWarnings );
    
    
    RC = theMMDBManager->ReadCoorFile (thePdb);
    // if RC > 0 reading in file has failed - say so and quit
    if (RC) {
        CXXException theException = CXXException ("ERROR in: CXXCreator::CXXCreator( pstr thePdb) - could not read pdb file");
        throw theException;
    }
    
    init();
    int selHnd = selectAllAtoms();
    theMMDBManager->GetSelIndex(selHnd, SelAtom, nSelAtoms);
}

CXXCreator::CXXCreator (mmdb::Manager* theMMDBManager_in) {
    theMMDBManager = theMMDBManager_in;
    int selHnd = selectAllAtoms();
    theMMDBManager->GetSelIndex(selHnd, SelAtom, nSelAtoms);
}

CXXCreator::CXXCreator (mmdb::Manager* theMMDBManager, int selHnd, int context_selHnd ) {
    init();
    theMMDBManager->GetSelIndex(selHnd, SelAtom, nSelAtoms);
    
    int neighbour_selhnd;
    int nSelAtomsNeighbours;
    mmdb::Atom** SelAtomNeighbours;
    
    //theMMDBManager->Select(neighbour_selhnd,mmdb::STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","*","*","*",mmdb::SKEY_NEW);
    
    /* 35 angstroms seems overly generous, but I still get incorrect charges with ca 25 angstroms. */
    if (context_selHnd>0) {
        neighbour_selhnd = theMMDBManager->NewSelection();
        theMMDBManager->SelectNeighbours(neighbour_selhnd,mmdb::STYPE_ATOM,SelAtom,nSelAtoms,0.0,35.0,mmdb::SKEY_OR);
        theMMDBManager->Select(neighbour_selhnd,mmdb::STYPE_ATOM,context_selHnd,mmdb::SKEY_AND);
        theMMDBManager->GetSelIndex(neighbour_selhnd, SelAtomNeighbours, nSelAtomsNeighbours);
        //theMMDBManager->DeleteSelection(neighbour_selhnd);
        SelAtom = SelAtomNeighbours;
        nSelAtoms = nSelAtomsNeighbours;
    } else {
        theMMDBManager->GetSelIndex(selHnd, SelAtom, nSelAtoms);
    }
}

void CXXCreator::init(){
    // defaults can be changed using setGridSpacing and set probeRadius but only before space is made ..
    gridSpacing = 1.2;
    probeRadius = 1.6; // WARNING need to have these changable depending on sor options ect. ..
    zeroSpace = 6;
    
    antiAliasPara = 5;
    scale = 1.0/(antiAliasPara);
    
    
    // default values for this two parameters ...
    externalDielectric = 78;
    internalDielectric = 4;
    
    theChargeTable= CXXChargeTable();
}

int CXXCreator::selectAllAtoms(){
    //make new selection. First create new selection handle...
    int selHnd = theMMDBManager->NewSelection();
    theMMDBManager->SelectAtoms( selHnd, 0,"*",mmdb::ANY_RES,"*",mmdb::ANY_RES,"*","*","*","*","*" );
    
    // first select all atoms in the pdb
    // !! WARNING: This does not exclude HETATOM records - how do deal with this ??
    return selHnd;
}

int CXXCreator::setParameters( double IonicStrength, double Temperature, double gridConstant) {
    
    gridSpacing = gridConstant;
    temp = Temperature;
    ionicStrength = IonicStrength;
    
    return 0;
}



CXXCoord<CXXCoord_ftype>CXXCreator::getAtomCoord(int atomNr) {
    
    CXXCoord<CXXCoord_ftype>theCoord;
    
    if (atomNr >= nSelAtoms) {
        CXXException theException = CXXException("ERROR in: CXXCoord::getAtomCoord(atomNr) - atomNr out of range");
        throw theException;
    }
    if(SelAtom){
        mmdb::Atom* theAtom = SelAtom[atomNr];
        if(theAtom){
            theCoord.setX(theAtom->x);
            theCoord.setY(theAtom->y);
            theCoord.setZ(theAtom->z);
        }
    }
    return theCoord;
}

double CXXCreator::getAtomRadius(int atomNr) {
    
    double radius = 0;
    
    if (atomNr >= nSelAtoms) {
        CXXException theException = CXXException("ERROR in: CXXCoord::getAtomRadius(atomNr) - atomNr out of range");
        throw theException;
    }
    if(SelAtom){
        mmdb::Atom* theAtom = SelAtom[atomNr];
        if(theAtom){
            radius = mmdb::getVdWaalsRadius(theAtom->element);
        }
    }
    return radius;
}


string CXXCreator::getAtomElement(int atomNr) {
    
    string theElement;
    
    if (atomNr >= nSelAtoms) {
        CXXException theException = CXXException("ERROR in: CXXCoord::getAtomElement(atomNr) - atomNr out of range");
        throw theException;
    }
    mmdb::Atom* theAtom = SelAtom[atomNr];
    
    theElement = theAtom->element;
    
    return theElement;
}


string CXXCreator::getAtomName(int atomNr) {
    
    string theName;
    if (atomNr >= nSelAtoms) {
        CXXException theException = CXXException("ERROR in: CXXCoord::getAtomName(atomNr) - atomNr out of range");
        throw theException;
    }
    mmdb::Atom* theAtom = SelAtom[atomNr];
    theName = theAtom->name;
    
    
    return theName;
}


string CXXCreator::getAtomResidueName(int atomNr) {
    
    mmdb::pstr name;
    string theResidueName;
    
    if (atomNr >= nSelAtoms) {
        CXXException theException = CXXException("ERROR in: CXXCoord::getAtomResidue(atomNr) - atomNr out of range");
        throw theException;
    }
    mmdb::Atom* theAtom = SelAtom[atomNr];
    name = theAtom->GetResName();
    theResidueName = name;
    
    return theResidueName;
}



double CXXCreator::lookUpCharge(int atomNr) {
    
    string residueName, atomName, element, charge;
    double theCharge = 0.;
    
    if(SelAtom){
        mmdb::Atom* theAtom = SelAtom[atomNr];
        if(theAtom){
            theCharge = SelAtom[atomNr]->charge;
        }
    }
    return theCharge;
    
    /*
     residueName = getAtomResidueName(atomNr);
     atomName = getAtomName(atomNr);
     element = getAtomElement(atomNr);
     
     double theAtomCharge;
     if ((theAtomCharge = SelAtom[atomNr]->charge) == 0.)
     theAtomCharge = theChargeTable.getCharge(residueName, atomName);
     return theAtomCharge;
     */
}



int CXXCreator::distributeAtomCharge(CXXCoord<CXXCoord_ftype>gridOrigin, CXXCoord<CXXCoord_ftype>xGridVector, CXXCoord<CXXCoord_ftype>yGridVector, CXXCoord<CXXCoord_ftype>zGridVector, int atomNr) {
    
    // this is very much the same as for finding the atom grid volume - only that now each grid point
    // gets charge associated with it:
    
    double gridSpacing = xGridVector.get3DLength(); // WARNING need to make sure this is consitent throughout the process
    double atomRadius = getAtomRadius(atomNr);		// For now only have equistant orthogonal grids ..
    double currentAtomCharge = 0;
    
    currentAtomCharge = lookUpCharge(atomNr);
    
    double currentAtomGridVolume = getGridVolumeOfAtom(gridOrigin,xGridVector,yGridVector,zGridVector, atomNr);
    double currentAtomChargeDensity = currentAtomCharge/currentAtomGridVolume;
    
    // a) find all gridpoints that are close to the atom center
    
    // find lower, back, left corner for box containing current atom
    CXXCoord<CXXCoord_ftype>theAtomLocation;
    theAtomLocation = getAtomCoord(atomNr);
    
    CXXCoord<CXXCoord_ftype>distanceVector = theAtomLocation - gridOrigin;	// gridOrigin has to be choosen so that this is positive !
    
    if (distanceVector.x() < 0 | distanceVector.y() < 0 | distanceVector.z() < 0) {
        CXXException theException = CXXException(" ERROR: (CXXCreator::distributeAtomCharge) :GridOrigin is not the samlleset point ...!\n");
        throw theException;
    }
    // the distance of the box edge from the edge of the grid is the distance vector
    // component (center of atom minus gird origin) in that direction minus the atom
    // radius - in gird units ...
    // take the distance vector - atomRadius in units of grid distance and round down ...
    
    distanceVector.setX(distanceVector.x() - atomRadius);
    distanceVector.setY(distanceVector.y() - atomRadius);
    distanceVector.setZ(distanceVector.z() - atomRadius);
    
    // now convert into grid coordinates ....
    
    distanceVector = distanceVector.scalarMultiply(1/gridSpacing);
    
    int iMin = static_cast<int>(distanceVector.x());
    int jMin = static_cast<int>(distanceVector.y());
    int kMin = static_cast<int>(distanceVector.z());
    
    // this is the number of grid units from grid edge to box edge
    // subtract one to add safety boundary ...
    
    iMin = iMin - 1;
    jMin = jMin - 1;
    kMin = kMin - 1;
    
    if (iMin <= 0 | jMin <= 0 | kMin < 0) {
        CXXException theException = CXXException(" ERROR: (CXXCreator::distributeAtomCharge) : AtomBox extends to or beyond grid boundary ...!\n");
        throw theException;
    }
    
    
    // these are the grid coordinates for the lower, left, back point of the box covering the atom
    
    // determine extend of box containing current atom plus one gridSpacing in all directions
    // box needs to contain at least the atom and one grid spaceing worth of space on each direction
    
    // first calculate extend of box atomRadius in grid units rounded down
    int intBoxSize = (static_cast<int>(2*atomRadius/gridSpacing));
    
    // then add two on each side to make sure there is at least one empty girdpoint
    intBoxSize = intBoxSize + 4;
    
    // finally get upper, right, front gridPoint of this box by adding intBoxSize to lower, left, back grid point
    int iMax = iMin + intBoxSize;
    int jMax = jMin + intBoxSize;
    int kMax = kMin + intBoxSize;
    
    if (iMax >= space->getDimI() | jMin >= space->getDimJ() | kMin >= space->getDimK()) {
        CXXException theException = CXXException(" ERROR: (CXXCreator::distributeAtomCharge) : AtomBox extends to or beyond grid boundary ...!\n");
        throw theException;
    }
    
    
    // start calculating volumes
    CXXCoord<CXXCoord_ftype>currentGridPoint;
    double fractionCovered;		// fraction of current grid point that is covered		double atomGridVolume = 0;		// total sum of all fractional grid volumes containg the atom
    double currentGridPointCharge = 0;
    double cubeVolume = xGridVector.get3DLength()*yGridVector.get3DLength()*zGridVector.get3DLength(); // warning - this assumes orthogonal grid ...
    
    // now iterate over all grid point in the box and find fractional cube coverage for each grid point
    for (int i = iMin; i < iMax; i++) {
        for (int j = jMin; j < jMax; j++) {
            for (int k = kMin; k < kMax; k++) {
                
                // make each gridpoint in turn and get fraction of it that is covered by theAtom
                currentGridPoint = (gridOrigin + xGridVector.scalarMultiply(i) + yGridVector.scalarMultiply(j) + zGridVector.scalarMultiply(k));
                currentGridPointCharge = 0; // (make this the current charge of the point ant then add...
                try {
                    fractionCovered = fractionOfGridPointCoveredByAtom(currentGridPoint,xGridVector,yGridVector,zGridVector, atomNr);
                    
                    // b) for each of those gridpoints assigen fraction of charge based on the degree of coverage
                    currentGridPointCharge = currentGridPointCharge + fractionCovered*cubeVolume*currentAtomChargeDensity;
                    // c) add this fraction of charge to any charge already present
                    
                    space->addGridCharge(i, j, k, currentGridPointCharge);
                }
                catch (CXXException theException) {
                    theException.Report();
                    throw theException;
                }
            }
        }
    }
    
    // here end of loop over atomBox
    return 0;
}


double CXXCreator::getGridVolumeOfAtom(CXXCoord<CXXCoord_ftype>gridOrigin, CXXCoord<CXXCoord_ftype>xGridVector, CXXCoord<CXXCoord_ftype>yGridVector, CXXCoord<CXXCoord_ftype>zGridVector, int atomNr) {
    
    
    // WARNING - want to do this with setOptions() type accessor ...
    int fastVolume = 1;
    if (fastVolume == 1) {
        
        double atomRadius;
        atomRadius = getAtomRadius(atomNr);
        double volume = 4.0*mmdb::Pi*atomRadius*atomRadius*atomRadius/3.0;
        return volume;
    }
    
    
    double cubeVolume = xGridVector.get3DLength()*yGridVector.get3DLength()*zGridVector.get3DLength();
    double gridSpacing = xGridVector.get3DLength();  // WARNING if allow non equidistant grids need to do that more carefuly
    
    double atomRadius;
    atomRadius = getAtomRadius(atomNr);
    
    // find lower, back, left corner for box containing current atom
    
    CXXCoord<CXXCoord_ftype>theAtomLocation;
    theAtomLocation = getAtomCoord(atomNr);
    
    CXXCoord<CXXCoord_ftype>distanceVector = theAtomLocation - gridOrigin;	// gridOrigin has to be choosen so that this is positive !
    
    // the distance of the box edge from the edge of the grid is the distance vector
    // component (center of atom minus gird origin) in that direction minus the atom
    // radius - in gird units ...
    // take the distance vector - atomRadius in units of grid distance and round down ...
    
    distanceVector.setX(distanceVector.x() - atomRadius);
    distanceVector.setY(distanceVector.y() - atomRadius);
    distanceVector.setZ(distanceVector.z() - atomRadius);
    
    if (distanceVector.x() < 0 | distanceVector.y() < 0 | distanceVector.z() < 0) {
        CXXException theException = CXXException(" ERROR: (CXXCreator::distributeAtomCharge) :GridOrigin is not the samlleset point ...!\n");
        throw theException;
    }
    
    
    // now convert into grid coordinates ....
    
    distanceVector = distanceVector.scalarMultiply(1/gridSpacing);
    
    
    int iMin = static_cast<int>(distanceVector.x());
    int jMin = static_cast<int>(distanceVector.y());
    int kMin = static_cast<int>(distanceVector.z());
    
    
    //std::cout << "Atom location x: " << theAtomLocation.x() << " y: " << theAtomLocation.y() << " z: " << theAtomLocation.z() << "\n";
    // this is the number of grid units from grid edge to box edge
    // subtract one to add safety boundary ...
    
    iMin = iMin - 1;
    jMin = jMin - 1;
    kMin = kMin - 1;
    // these are the grid coordinates for the lower, left, back point of the box covering the atom
    
    if (iMin <= 0 | jMin <= 0 | kMin < 0) {
        CXXException theException = CXXException(" ERROR: (CXXCreator::distributeAtomCharge) : AtomBox extends to or beyond grid boundary ...!\n");
        throw theException;
    }
    
    // determine extend of box containing current atom plus one gridSpacing in all directions
    // box needs to contain at least the atom and one grid spaceing worth of space on each direction
    
    // first calculate extend of box atomRadius in grid units rounded down
    int intBoxSize = (static_cast<int>(2*atomRadius/gridSpacing));
    
    // then add two on each side to make sure there is at least one empty girdpoint
    intBoxSize = intBoxSize + 4;
    
    // finally get upper, right, front gridPoint of this box by adding intBoxSize to lower, left, back grid point
    int iMax = iMin + intBoxSize;
    int jMax = jMin + intBoxSize;
    int kMax = kMin + intBoxSize;
    
    if (iMax >= space->getDimI() | jMin >= space->getDimJ() | kMin >= space->getDimK()) {
        CXXException theException = CXXException(" ERROR: (CXXCreator::distributeAtomCharge) : AtomBox extends to or beyond grid boundary ...!\n");
        throw theException;
    }
    
    
    // start calculating volumes
    CXXCoord<CXXCoord_ftype>currentGridPoint;
    double fractionCovered;		// fraction of current grid point that is covered
    double atomGridVolume = 0;		// total sum of all fractional grid volumes containg the atom
    
    
    // now iterate over all grid point in the box and find fractional cube coverage for each grid point
    for (int i = iMin; i < iMax; i++) {
        for (int j = jMin; j < jMax; j++) {
            for (int k = kMin; k < kMax; k++) {
                
                // make each gridpoint in turn and get fraction ofit that is covered by theAtom
                currentGridPoint = (gridOrigin + xGridVector.scalarMultiply(i) + yGridVector.scalarMultiply(j) + zGridVector.scalarMultiply(k));
                //std::cout << "Current Grid Pointx: " << currentGridPoint.x() << " y: " << currentGridPoint.y() << " z: " << currentGridPoint.z() << "\n";
                try {
                    fractionCovered = fractionOfGridPointCoveredByAtom(currentGridPoint,xGridVector,yGridVector,zGridVector, atomNr);
                    // Volume of current gridPoint associated space covered by current atom is fractonCovered*cubeVolume
                    // the sum of all the subcubes multiplied by the subcube volume is the total cubeVolume of the atom
                    atomGridVolume = atomGridVolume + fractionCovered*cubeVolume;
                }
                catch (CXXException theException) {
                    theException.Report();
                    throw theException;
                }
            }
        }
    }
    
    return atomGridVolume;
}

double CXXCreator::fractionOfGridPointCoveredByAtom(CXXCoord<CXXCoord_ftype>theGridPoint, CXXCoord<CXXCoord_ftype>xGridVector, CXXCoord<CXXCoord_ftype>yGridVector, CXXCoord<CXXCoord_ftype>zGridVector, int atomNr) {
    
    
    // determines the how many subgrid points cover the neighbourhood of the central gird point
    
    
    double fraction;
    double atomRadius;
    
    atomRadius = getAtomRadius(atomNr);
    double atomRadiusSq = atomRadius * atomRadius;
    
    
    CXXCoord<CXXCoord_ftype>theAtomLocation;
    theAtomLocation = getAtomCoord(atomNr); // pass by reference or by value ?
    
    CXXCoord<CXXCoord_ftype>subGridPoint;
    CXXCoord<CXXCoord_ftype>distanceVector;
    
    // this is the - - - corner of the grid used for charge anti aliasing. The back, low, left corner of a
    // cube of grid spacing size centered around the gridpoint
    // The subgrid is defined to cover this cube- the amount charge assigned to the grid point depends on
    // the number of these subgrid points that are covered by the atom sphere.
    // Reference for this procedure: Bruccoleri et al J. Phys. Chem 1997
    
    CXXCoord<CXXCoord_ftype>lowerBackLeftSubgridCorner = theGridPoint - xGridVector.scalarMultiply(0.5)
    - yGridVector.scalarMultiply(0.5) - zGridVector.scalarMultiply(0.5);
    
    // grid step vectors in subgrid - resolution determined by scale that is antiAliasPara ...
    CXXCoord<CXXCoord_ftype>xStepVector = xGridVector.scalarMultiply(scale);
    CXXCoord<CXXCoord_ftype>yStepVector = yGridVector.scalarMultiply(scale);
    CXXCoord<CXXCoord_ftype>zStepVector = zGridVector.scalarMultiply(scale);
    
    // start at the CENTER of the first subcube ...
    CXXCoord<CXXCoord_ftype>centerOfFirstSubcube = lowerBackLeftSubgridCorner + xStepVector.scalarMultiply(0.5)
    + yStepVector.scalarMultiply(0.5) + zStepVector.scalarMultiply(0.5);
    
    
    // this counts the subgrid points inside the atom sphere
    int subGridCounter = 0;
    
    // now loop over whole subgrid around the gridpoint given
    
    for (int i = 0; i < antiAliasPara; i++) {
        for (int j = 0; j < antiAliasPara; j++) {
            for (int k = 0; k < antiAliasPara; k++) {
                
                // calculate each point in subgrid and look if it is in the interior of the atom sphere.
                subGridPoint = centerOfFirstSubcube + xStepVector.scalarMultiply(i)
                + yStepVector.scalarMultiply(j) + zStepVector.scalarMultiply(k);
                
                // calculate distance from subgrid point to atom location
                distanceVector = theAtomLocation - subGridPoint;
                
                // calculate square length of distance vector and compare to atomRadius
                // if subgrid point is inside the atom spere - count subGridCounter up by one
                // <=> more charge gets distributed to the gridPoint at the center of the subgrid
                if (distanceVector*distanceVector <= atomRadiusSq) {
                    subGridCounter++;
                }
            }
        }
    }
    
    
    // this the fraction of the grid size cube around the gridPoint that is covered by the volume of the atom
    
    fraction = double (subGridCounter*1.0/(antiAliasPara*antiAliasPara*antiAliasPara));
    
    return fraction;
    
}



int CXXCreator::createSpace() {
    
    // first find out how big space needs to be to contain the content of the pdb
    // this could allow optimisation by adding a rotation function to make box as small as possible !
    
    // loop over all atoms and find max and min coords in all directions
    
    
    CXXCoord<CXXCoord_ftype>currentAtomCoord;
    
    double xMin =1e10;
    double xMax =-1e10;
    double yMin =1e10;
    double yMax =-1e10;
    double zMin =1e10;
    double zMax =-1e10;
    double rMax =0;
    
    
    for (int atomNr = 0; atomNr < nSelAtoms; atomNr++) {
        currentAtomCoord = getAtomCoord(atomNr);
        rMax = max(rMax, getAtomRadius(atomNr));
        xMin = min(currentAtomCoord.x(),CXXCoord_ftype(xMin));
        yMin = min(currentAtomCoord.y(),CXXCoord_ftype(yMin));
        zMin = min(currentAtomCoord.z(),CXXCoord_ftype(zMin));
        xMax = max(currentAtomCoord.x(),CXXCoord_ftype(xMax));
        yMax = max(currentAtomCoord.y(),CXXCoord_ftype(yMax));
        zMax = max(currentAtomCoord.z(),CXXCoord_ftype(zMax));
    }
    
    // now choose conservative box size and add empty boundary of size zeroSpace on all sides (depends on
    // poisson boltzmann boundary condition option choosen)
    
    xMin = xMin - rMax - zeroSpace;
    yMin = yMin - rMax - zeroSpace;
    zMin = zMin - rMax - zeroSpace;
    
    xMax = xMax + rMax + zeroSpace;
    yMax = yMax + rMax + zeroSpace;
    zMax = zMax + rMax + zeroSpace;
    
    
    
    // call creator of space with the new coords to make space of the right size
    
    try {
        space = new CXXSpace(probeRadius, gridSpacing, xMin, xMax, yMin, yMax, zMin, zMax);
        space->defineBoundaryConditions(externalDielectric);
    }
    catch (CXXException theException) {
        theException.Report();
        return 1;
    }
    // If this went through we now have empty space to add charges, dielectric and potential to ...
    
    
    // synchronise parameters
    space->setSolventParameters(ionicStrength, temp);
    return 0;
    
}


int CXXCreator::introduceMatter(double dielectricOfSolvent, double dielectricOfProtein) {
    
    
    externalDielectric = dielectricOfSolvent;
    internalDielectric = dielectricOfProtein;
    
    std::cout << "\nNow inserting matter into space\n";
    
    int chargedAtomCount = 0;
    try {
        CXXCoord<CXXCoord_ftype>xGridVector = space->getSpaceSpanningVectorX();
        CXXCoord<CXXCoord_ftype>yGridVector = space->getSpaceSpanningVectorY();
        CXXCoord<CXXCoord_ftype>zGridVector = space->getSpaceSpanningVectorZ();
        
        
        CXXCoord<CXXCoord_ftype>gridOrigin = space->getOrigin();// grid origin is choosen so that all grid indices are positive ...
        
        // loop over all atoms and consider atoms with non zero charge
        
        double theAtomCharge;
        double theAtomRadius;
        CXXCoord<CXXCoord_ftype>theAtomCoord;
        std::cout << "Adding atoms\n";
        
        for (int atomNr = 0; atomNr < nSelAtoms; atomNr++) {
            
            theAtomCoord = getAtomCoord(atomNr);
            theAtomRadius = getAtomRadius(atomNr);
            
            theAtomCharge = lookUpCharge(atomNr);
            
            // adding each atom to space to generate fft based map to distinguish between inside and outside
            //std::cout << "Adding atom at: " << theAtomCoord.x() << " " << theAtomCoord.y() << " " << theAtomCoord.z() << "\n atomRadius: " << theAtomRadius << "\n " ;
            space->atomAt3f( theAtomCoord.x(), theAtomCoord.y(), theAtomCoord.z(),  theAtomRadius);
            
            // only run through this if atom adds charge to the system
            if (theAtomCharge != 0) {
                chargedAtomCount++;
                distributeAtomCharge(gridOrigin, xGridVector, yGridVector, zGridVector, atomNr);
            }
        }
        std::cout << nSelAtoms << " total atoms inserted into space.\n" << chargedAtomCount << " charged atoms inserted\n";
        std::cout << "Charge smoothing: charge anti aliasing\n";
        
        
        space->introduceMedium(externalDielectric, internalDielectric, probeRadius);
        
        //			space->dumpSpaceSlice(0,2,9);
        //			space->dumpSpaceSlice(0,1,9);
        
        
    }
    catch (CXXException theException) {
        std::cout << "Reporting exception terminating this run: \n";
        theException.Report();
        throw theException;
        return 1;
    }
    
    
    return 0;
}


int CXXCreator::evolve(int optionSOR, double convergenceCriterion) {
    
    /* when this is called the system is described by the four discrete properties in space:
     
     i)		Its charge distribution stored in the dim[0]xdim[1]xdim[2] grid: chargeGrid
     at each point i,j,k the charge density on that grid point can be accessed through
     
     getGridCharge(i,j,k)
     
     ii)		The dielectricGrid: dielGrid. At a given point i,j,k the dielGrid has three components:
     
     getDielGrid(i,j,k,0) - dielectric looking from i,j,k to i+1,j,k
     getDielGrid(i,j,k,1) - dielectric looking from i,j,k to i,j+1,k
     getDielGrid(i,j,k,2) - dielectric looking from i,j,k to i+1,j,k+1
     
     
     iii)	The electrostatic potential also is defined at each point i,j,k but is initialised to zero.
     
     accessor:			getPotential(i,j,k)
     
     
     
     The relationship between charge, medium (dielectric) and potential is described by the Poisson Botzmann
     equation. Linearised and expressed as finite differece equation in terms of i), ii) and iii) the PBE
     reads:
     
     charge(i,j,k)/h =   dielGrid(i-1,j,k,0)[potentialGrid(i,j,k) - potentialGrid(i-1,j,k)]
     +	dielGrid(i,j-1,k,0)[potentialGrid(i,j,k) - potentialGrid(i,j-1,k)]
     +	dielGrid(i,j  ,k,0)[potentialGrid(i,j,k) - potentialGrid(i,j+1,k)]
     +   dielGrid(i,j,k-1,0)[potentialGrid(i,j,k) - potentialGrid(i,j,k-1)]
     +	dielGrid(i,j,k 1,0)[potentialGrid(i,j,k) - potentialGrid(i,j,k-1)]
     +   gridSpacing*gridSpacing*dielectric*K(i,j,k)*K(i,j,k)*potentialGrid(i,j,k)
     
     where: K(i,j,k) is the Debey-Hueckel parameter at point i,j,k (dependent on salt concentration at this point)
     
     
     
     Instead of interating this directly =>O(N*N*N) we use the successive over-relaxation approach [???]
     
     Gauss-Seidel ...
     
     
     Analytially it can be shown:		w = 2/( 1 + sqrt( 1 - lambda))
     
     where lambda is the minimal spectral radius of the Gauss-Seidel matrix
     
     
     Ways to make this fast:
     
     1)  optimal omega
     
     The spectral radius can be approximated using the connected moment expansion (CMX - later) or
     approximative formula
     
     2)  chebyshev acceleration
     
     The optima over-relaxation paramber calculated in 1) is only optimal close to convergence. The process
     can be accelerated by having w converge from 1 to optimal value throughout the iteration:
     
     => w = w(n)
     
     with:			w(1)  = 1 / (1 - 0.5*lambda*lambda)
     for n > 1:		w(n+1) = 1/ (1 - 0.25*lambda*lambda*w(n) )
     
     3)  finally use "odd/even" pattern in indeces - that is only need to update even/odd i,j,k point every
     every other pass:
     
     update (i,j,k) only if (i+j+k)mod 2 == n mod 2
     
     */
    
    
    // need the numbers to be real so now introduce constants of significance ..
    
    // elementary charge in Coulomb
    double e = 1.6022e-19;
    // vacuum periability in Farad / m = C/(V*m)
    double epsilonZero=8.85e-12;
    
    
    // Use approximate formula to derive spectral radius of this problem based on dimensions of grid
    
    double lambda = 0;
    lambda = 1.0/3.0* (cos(mmdb::Pi/(space->getDimI())) + cos(mmdb::Pi/(space->getDimJ())) + cos(mmdb::Pi/(space->getDimK())) );
    std::cout << "\nSpecctral radius of problem approximated as: " << lambda << "\n";
    
    // caclculating initial over relaxation parameter
    double w = 1/(1-0.5*lambda*lambda);
    std::cout << "=> approximate optimal first step over relaxation parameter: " << w << "\n";
    
    // houskeeping variables for SOR algorithm
    //double currentPotential;
    //double newPotential;
    //double resid =0;
    int levelReported = 0;
    // keeps track of the norm for convergence criterion
    double largestChange = -1e30;
    double largestPotential = -1e30;
    double fractionalCorrection = 1;
    
    // components of current potential point 6 nearest epsilons + kappa factor)
    //double eFactor =0;
    // grid Factor puts together all the constants that need to be multiplied to the charge to give the source term
    double chargeFactor = e/(gridSpacing*1e-10*epsilonZero);
    // counts interation steps
    int n = 1;
    int maxIterations = 200;
    std::cout << "\nStarting SOR iteration \n\n";
    // now iterate over whole grid
    
    try {
        
        while (fabs(fractionalCorrection) > convergenceCriterion && n < maxIterations) { // while potential is still changeing...
            // reset largestChange
            largestChange = 0;
            
            // iterate over all grid points in the grid interior
            for (size_t i = 0; i < space->getDimI(); i++) {
                
                double largestChangeOfRow[space->getDimJ()];
                double *largestChangeOfRowPntr = &(largestChangeOfRow[0]);
                for (size_t j=0; j<space->getDimJ(); j++) largestChangeOfRow[j]=-1e30;
                double largestPotentialOfRow[space->getDimJ()];
                double *largestPotentialOfRowPntr = &(largestPotentialOfRow[0]);
                for (size_t j=0; j<space->getDimJ(); j++) largestPotentialOfRow[j]=-1e30;
                CXXSpace *passableSpace = space;
                
#if __APPLE__ && !defined _OPENMP
                dispatch_apply(space->getDimJ(), dispatch_get_global_queue(0, 0), ^(size_t j){
#elif defined _OPENMP
#pragma omp parallel for default(none) shared(i, passableSpace, w, largestPotentialOfRowPntr, largestChangeOfRowPntr, n, chargeFactor) schedule(dynamic, 10)
#warning Compiling for OMP
                for (size_t j = 0; j < passableSpace->getDimJ(); j++) {
#else
                for (size_t j = 0; j < passableSpace->getDimJ(); j++) {
#endif

                    for (size_t k = 0; k < passableSpace->getDimK(); k++) {
                        
                        // odd/even ordering - only update when required ..
                        if (((i+j+k)%2) == (n%2)) {
                            
                            // look up old potential at this grid point
                            double currentPotential = passableSpace->getPotential(i,j,k);
                            
                            // iteration gives new (n+1) potenial WARNING: performance- these can be combined into kappa grid
                            
                            // eFactor =(  space->getDielGrid(i,j,k,0) + space->getDielGrid(i-1,j,k,0) +
                            //			space->getDielGrid(i,j,k,1) + space->getDielGrid(i,j-1,k,1) +
                            //			space->getDielGrid(i,j,k,2) + space->getDielGrid(i,j,k-1,2)  + kapapaFactor);
                            
                            double eFactor = passableSpace->getEpsilonKappaSq(i,j,k);
                            
                            // eFactor*currentPotential = (sum(6 nearest Epsilons) + kappaFactor) sourceTerms: charge*chargeFactor
                            double resid = eFactor*currentPotential - passableSpace->getGridCharge(i,j,k)*chargeFactor
                            - passableSpace->getDielGrid(i-1,j,k,0)*passableSpace->getPotential(i-1,j,k) - passableSpace->getDielGrid(i,j,k,0)*space->getPotential(i+1,j,k)
                            - passableSpace->getDielGrid(i,j-1,k,1)*passableSpace->getPotential(i,j-1,k) - passableSpace->getDielGrid(i,j,k,1)*space->getPotential(i,j+1,k)
                            - passableSpace->getDielGrid(i,j,k-1,2)*passableSpace->getPotential(i,j,k-1) - passableSpace->getDielGrid(i,j,k,2)*passableSpace->getPotential(i,j,k+1);
                            
                            // use residual and over correct current potenttial by over relaxation to calculate new potential
                            double newPotential  = currentPotential - w*resid/eFactor;
                            
                            // update the potential at the current grid point
                            passableSpace->setPotential(i,j,k,newPotential);
                            
                            // see if the difference between the current and the new potential is bigger than the largest
                            // difference observed for this n (for convergence)
                            largestPotentialOfRowPntr[j] = max(largestPotentialOfRowPntr[j], double(fabs(currentPotential)));
                            largestChangeOfRowPntr[j] = max(largestChangeOfRowPntr[j], double(fabs(newPotential - currentPotential)));
                            
                        }
                    }
                    //}

#if __APPLE__ && !defined _OPENMP
                });
#else
                }
#endif
                // for (size_t j=0; j<space->getDimJ(); j++) largestChange = max(largestChange,largestChangeOfRow[j]);
                // for (size_t j=0; j<space->getDimJ(); j++) largestPotential = max(largestPotential,largestPotentialOfRow[j]);
            }
            
            
            
            
            // calculate largest fractional change
            fractionalCorrection = largestChange/largestPotential;
            // give some feedback ..
            if (fractionalCorrection < 0.25 && levelReported < 1) {
                std::cout << "25% convergence level reached \n";
                levelReported = 1;}
            
            if (fractionalCorrection < 0.1 && levelReported < 2) {
                std::cout << "10% convergence level reached \n";
                levelReported = 2;}
            
            if (fractionalCorrection < 0.05 && levelReported < 3) {
                std::cout << "5% convergence level reached \n";
                levelReported = 3;					}
            
            if (fractionalCorrection < 0.01 && levelReported < 4) {
                std::cout << "1% convergence level reached \n";
                levelReported = 4;					}
            
            if (fractionalCorrection < 0.005 && levelReported < 5) {
                std::cout << "0.5% convergence level reached \n";
                levelReported = 5;					}
            
            if (fractionalCorrection < 0.001 && levelReported < 6) {
                std::cout << "0.1% convergence level reached \n";
                levelReported = 6;					}
            
            // end of "while" loop - if the largestChange is still larger than the convergenceCriterion
            // then:	increment n
            //			iterate w
            //			reset largestChange (above) and start again ...
            
            
            w = 1/(1-0.25*lambda*lambda*w);
            n++;
        }
    }
    catch (CXXException theException) {
        theException.Report();
        return 1;
    }
    if (n >= maxIterations) {
        std::cout << "WARNING: max number of iterations reached - convergence to specified level not reached !\n";
        std::cout << fractionalCorrection*100 << " % convergence reached after " << n << " iterations. \n";
        CXXException theException("WARNING: max number of iterations reached - convergence to specified level not reached !\n");
        throw theException;
        return 1;
    }
    else {
        std::cout << "\n" << fractionalCorrection*100 << "% convergence reached after " << n << " iterations. \n";
        //		space->dumpSpaceSlice(0,3,9);
        return 0;
    }
}

double CXXCreator::getGridSpacing(){
    return gridSpacing;
}

int CXXCreator::calculate(){
    setParameters(0.15, 300, 1.2);
    createSpace();
    introduceMatter(78,4);
    evolve(0, 0.01);
    return 0;
}

clipper::NXmap<double> CXXCreator::coerceToClipperMap(clipper::Cell &cell){
    //Coerce it into a clipper NXmap:
    //First evaluate the limits, and corresponding cell dimensions etc
    double spacing = getGridSpacing();
    CXXCoord<CXXCoord_ftype>doubleOrigin = space->getOrigin();
    int iu1 = int(doubleOrigin.x()/spacing+(doubleOrigin.x()>=0?0.4:-0.4));
    int iv1 = int(doubleOrigin.y()/spacing+(doubleOrigin.y()>=0?0.4:-0.4));
    int iw1 = int(doubleOrigin.z()/spacing+(doubleOrigin.z()>=0?0.4:-0.4));
    int nu = int(space->getDimI());
    int nv = int(space->getDimJ());
    int nw = int(space->getDimK());
    int iu2 = iu1+nu-1;
    int iv2 = iv1+nv-1;
    int iw2 = iw1+nw-1;
    float a=(nu+1)*spacing;
    float b=(nv+1)*spacing;
    float c=(nw+1)*spacing;
    float alpha, beta, gamma;
    alpha = beta = gamma = 90.;
    
    //Now create the clipper objects needed to instantiate a map
    cell = clipper::Cell(clipper::Cell_descr(a, b, c, alpha, beta, gamma));
    clipper::Grid_sampling grid(nu+1, nv+1, nw+1);
    clipper::Grid_range grid_extent(clipper::Coord_grid(iu1, iv1, iw1), clipper::Coord_grid(iu2, iv2, iw2));
    //instantiate the NXmap
    clipper::NXmap<double> theClipperMap(cell, grid, grid_extent);
    
    //And fill it with data from our own internal data representation
    for (int i=0; i<nu; i++){
        for (int j=0; j<nv; j++){
            for (int k=0; k<nw; k++){
                theClipperMap.set_data(clipper::Coord_grid(i, j, k), space->getPotential(i, j, k));
            }
        }
    }
    return theClipperMap;
}
