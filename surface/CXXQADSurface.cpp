#include "CXXQADSurface.h"
#include <map>
#include <algorithm>
#include <vector>
#define fmax(A,B) ((A)>(B)?(A):(B));
#define fmin(A,B) ((A)<(B)?(A):(B));
#include <time.h>
using namespace clipper;
using namespace std;

CXX_mot::CXXQADSurface::~CXXQADSurface() {}

CXX_mot::CXXQADSurface::CXXQADSurface(mmdb::PManager theMMDBManager_in, int selHndl_in, 
							 double probeRadius_in, double sample_in){
	
	theMMDBManager = theMMDBManager_in;
	selHndl = selHndl_in;
	probeRadius = probeRadius_in;
	sample = sample_in;
	clock_t time = clock();
	clock_t newTime;
	
	//isolate the array of selected atoms
	theMMDBManager->GetSelIndex(selHndl, selectedAtoms, nSelectedAtoms);
	
	//determine the maximum atomic radius
	maxAtomRadius = 0.;
	for (int i=0; i< nSelectedAtoms; i++){
		atomRadii.push_back(getAtomRadius(selectedAtoms[i]));
		maxAtomRadius = max(maxAtomRadius, getAtomRadius(selectedAtoms[i]));
	}
	
	//Pre identify all contacts 
	mmdb::Contact *contacts = NULL;
	contacts = 0;
	int nContacts = 0;
	cout << "Off to precalculate contacts..."; cout.flush();
	theMMDBManager->SeekContacts( selectedAtoms, nSelectedAtoms, 
								  selectedAtoms, nSelectedAtoms, 
								  0., 2.*maxAtomRadius+2.*probeRadius, 
								  0, contacts, nContacts, 0, 0);
	neighbourhoods.resize(nSelectedAtoms);
	for (int i=0; i< nContacts; i++){
		neighbourhoods[contacts[i].id1].push_back(contacts[i].id2);
	}
	newTime = clock(); double deltaTime = (newTime - time) / CLOCKS_PER_SEC; time = newTime;
	cout << "...Done " << deltaTime << nContacts << " found\n"; cout.flush();	
	if (contacts) delete (contacts);
	
	//determine size of grids to be used in surface calculation and make the maps
	cout << "Off to prepare grids..."; cout.flush();
	prepareGrids();
	newTime = clock();
	newTime = clock();  deltaTime = double(newTime - time) / CLOCKS_PER_SEC; time = newTime;
	cout << "...Done " << deltaTime << "\n"; cout.flush();
	
	//make first map which reflects at each point the probebal distance to
	//a nearby solvent probe position
	cout << "Off to make a distance map..."; cout.flush();
	makeDistanceSqMap();
	newTime = clock();  deltaTime = double(newTime - time) / CLOCKS_PER_SEC; time = newTime;
	cout << "...Done " << deltaTime <<"\n"; cout.flush();
	
	//Set the distanceSq for edible pixels to (probeRadius + sample)^2
	cout << "Off to modify a distance map..."; cout.flush();
	setInaccessibleDistanceSq();
	newTime = clock();  deltaTime = double(newTime - time) / CLOCKS_PER_SEC; time = newTime;
	cout << "...Done " << deltaTime << "\n"; cout.flush(); 
	
	//make some candidate probe positions that will exactly define the extrapolated
	//probePositions from the VDW surface element
	cout << "Off to deal with vdW surface..."; cout.flush();
	addProbesFromVdwSurface();
	newTime = clock();  deltaTime = double(newTime - time) / CLOCKS_PER_SEC; time = newTime;
	cout << "...Done " << deltaTime << " "<< probePositions.size() << " probe positions now\n"; cout.flush();

	//Now identify torus sections and allow them also to eat into the volume
	cout << "Off to deal with tori..."; cout.flush();
	toruses();
	newTime = clock();  deltaTime = double(newTime - time) / CLOCKS_PER_SEC; time = newTime;
	cout << "...Done " << deltaTime <<" "<< probePositions.size() << " probe positions now\n"; cout.flush();

	//Now allow the probes to eat into the solvent surface
	cout << "Off to allow probes to eat into the molecular volume..."; cout.flush();
	allowProbesToEat();
	newTime = clock();  deltaTime = double(newTime - time) / CLOCKS_PER_SEC; time = newTime;
	cout << "...Done " << deltaTime << "\n"; cout.flush();
	
	//Now convert from distSq map to dist map
	cout << "Converting distSq to dist..."; cout.flush();
	sqrtDistanceSq();	newTime = clock();  deltaTime = double(newTime - time) / CLOCKS_PER_SEC; time = newTime;
	cout << "...Done " << deltaTime << "\n"; cout.flush();
	
	//Now contour the volume
	cout << "Contour the map..."; cout.flush();
	contourMap(probeRadius-1e-06);
	newTime = clock();  deltaTime = double(newTime - time) / CLOCKS_PER_SEC; time = newTime;
	cout << "...Done " << deltaTime << "\n"; cout.flush();
	
	//Now assign normals to each vertex based on the triangles that share it
	cout << "Off to calculate averaged normals..."; cout.flush();
	calculateAveragedNormals();	
	newTime = clock();  deltaTime = double(newTime - time) / CLOCKS_PER_SEC; time = newTime;
	cout << "...Done " << deltaTime << "\n"; cout.flush();

	//assignNormalsFromProbes();
	//Make prettier by using normals derived from actual atom positions where permissible
	//cout << "Off to apply atom based normals..."; cout.flush();
	//applyAtomBasedNormals();
	//cout << "...Done\n"; cout.flush();
}

int CXX_mot::CXXQADSurface::prepareGrids (){
	//Pass through finding the limits of the coordinates
	
	double xyzMin[3], xyzMax[3];
	for (int i=0; i<3; i++){
		xyzMin[i] = 1e30;
		xyzMax[i] = -1e30;
	}
	for (int i=0; i< nSelectedAtoms; i++){
		xyzMin[0] = fmin(xyzMin[0], selectedAtoms[i]->x);
		xyzMin[1] = fmin(xyzMin[1], selectedAtoms[i]->y);
		xyzMin[2] = fmin(xyzMin[2], selectedAtoms[i]->z);
		xyzMax[0] = fmax(xyzMax[0], selectedAtoms[i]->x);
		xyzMax[1] = fmax(xyzMax[1], selectedAtoms[i]->y);
		xyzMax[2] = fmax(xyzMax[2], selectedAtoms[i]->z);
	}
	
	//Expand max and min coordinates by maximum estimated atom radius + probe radius + 1 grid point
	for (int i=0; i<3; i++){
		xyzMin[i] -= (maxAtomRadius + 2.*probeRadius + 2.*sample);
		xyzMax[i] += (maxAtomRadius + 2.*probeRadius + 2.*sample);
	}
	
	//Convert to grid coordinates	
	double uvwMin[3], uvwMax[3];
	for (int i=0; i<3; i++){
		uvwMin[i] = xyzMin[i] / sample;
		uvwMax[i] = xyzMax[i] / sample;
	}
	
	//Convert to grid unit integers such that the range entirely contains the coordinate limits
	int iUvwMin[3], iUvwMax[3], nUvw[3];
	for (int i=0; i<3; i++){
		iUvwMin[i] = int(uvwMin[i] - (uvwMin[i]<=0.?1.:0.));
		iUvwMax[i] = int(uvwMax[i] + (uvwMax[i]>=0.?1.:0.));
		nUvw[i] = (iUvwMax[i] - iUvwMin[i]) + 1;
	}

	//Generate a pseudo cell that is slightly larger than this box
	double a = (nUvw[0]+1) * sample;
	double b = (nUvw[1]+1) * sample;
	double c = (nUvw[2]+1) * sample;
	
	//Now create the clipper objects needed to instantiate a map
	clipperSpacegroup = Spacegroup(Spacegroup::P1);
	clipperCell = clipper::Cell(clipper::Cell_descr(a, b, c, 90., 90., 90.));
	clipperGridSampling = clipper::Grid_sampling(nUvw[0]+1, nUvw[1]+1, nUvw[2]+1);
	clipperGridRange = clipper::Grid_range(Coord_grid(iUvwMin[0], iUvwMin[1], iUvwMin[2]),
										   Coord_grid(iUvwMax[0], iUvwMax[1], iUvwMax[2]));			
	
	//dump(clipperGridSampling);
	//dump(clipperGridRange);
	
	//instantiate the Xmaps
	
	theDoubleMap = clipper::Xmap<double>(clipperSpacegroup, clipperCell, clipperGridSampling);	
	theFlagMap = clipper::Xmap<int>(clipperSpacegroup, clipperCell, clipperGridSampling);	
	clipper::Xmap_base::Map_reference_index ix;
	
	// for ( ix = theDoubleMap.first(); !ix.last(); ix.next() ) {
	//	 theFlagMap[ix] = Solvent;
	//	 theDoubleMap[ix] = 0.;
	// }

	return 0;
}

int CXX_mot::CXXQADSurface::makeDistanceSqMap(){
	Xmap<double>::Map_reference_coord i0, iu, iv, iw;
	
	for (int iAtom = 0; iAtom < nSelectedAtoms; iAtom++){
		double vdWRadius = fastGetAtomRadius(iAtom) ;
		double accessibleRadius = vdWRadius + probeRadius;
		double accessibleRadiusSq = accessibleRadius*accessibleRadius;

		Grid_range gd (clipperCell, clipperGridSampling, accessibleRadius);
		
		Coord_orth atomCoordOrth(selectedAtoms[iAtom]->x, 
								 selectedAtoms[iAtom]->y, 
								 selectedAtoms[iAtom]->z);
		
		Coord_frac uvw = atomCoordOrth.coord_frac( clipperCell);
		Coord_grid g0 = uvw.coord_grid(clipperGridSampling) + gd.min();
		Coord_grid g1 = uvw.coord_grid(clipperGridSampling) + gd.max();
		
		i0 = Xmap<double>::Map_reference_coord( theDoubleMap, g0 );
		
		for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() ){
			for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() ){
				for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ){
					Coord_orth gridCoord = iw.coord_orth();
					double dx = gridCoord[0] - atomCoordOrth[0];
					if (fabs(dx) < accessibleRadius){
						double dy = gridCoord[1] - atomCoordOrth[1];
						if (fabs(dy) < accessibleRadius){
							double dxdysq = dx*dx + dy*dy;
							if (dxdysq < accessibleRadiusSq){
								double dz = gridCoord[2] - atomCoordOrth[2];
								if (fabs(dz) < accessibleRadius){
									double dxdydzsq = dz*dz+dxdysq;
									if (dxdydzsq < accessibleRadiusSq){
										double distanceFromAtom = sqrt(dxdydzsq);
										double distanceFromProbe = accessibleRadius - distanceFromAtom;
										double distanceFromProbeSq = distanceFromProbe * distanceFromProbe;
										if (distanceFromProbeSq > theDoubleMap[iw]){
											theDoubleMap[iw] = distanceFromProbeSq;
											switch (theFlagMap[iw]) {
//												if (distanceFromProbeSq<vdWRadiusSq) theFlagMap[iw] = vdW;
//												else theFlagMap[iw] = Accessible;
//												break;
												case Solvent:
												case vdW:
												case Inaccessible:
												case Accessible:
													theFlagMap[iw] = Inaccessible;
													break;
											}
										}	
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return 0;
}

int CXX_mot::CXXQADSurface::addProbesFromVdwSurface(){
	Xmap<double>::Map_reference_coord i0, iu, iv, iw;
	
	for (int iAtom = 0; iAtom < nSelectedAtoms; iAtom++){
		double vdwRadius = fastGetAtomRadius(iAtom);

		Coord_orth atomCoordOrth(selectedAtoms[iAtom]->x, 
								 selectedAtoms[iAtom]->y, 
								 selectedAtoms[iAtom]->z);
		Coord_frac uvw = atomCoordOrth.coord_frac( clipperCell);
		
		Grid_range gd (clipperCell, clipperGridSampling, vdwRadius);
		Coord_grid g0 = uvw.coord_grid(clipperGridSampling) + gd.min();
		Coord_grid g1 = uvw.coord_grid(clipperGridSampling) + gd.max();

		Grid_range gd1 (clipperCell, clipperGridSampling, vdwRadius+probeRadius);
		Coord_grid g2 = uvw.coord_grid(clipperGridSampling) + gd1.min();
		Coord_grid g3 = uvw.coord_grid(clipperGridSampling) + gd1.max();
		Grid_range atomRange = Grid_range(g2, g3);
		
		i0 = Xmap<double>::Map_reference_coord( theDoubleMap, g0 );
		double vdwRadiusSq = vdwRadius*vdwRadius;
		double extrapolate = ((probeRadius+vdwRadius)/vdwRadius);
		
		//Now loop over u and v dimensions, calculating w coordinate which will be on the
		//vdW surface of the atom
		nVdwProbePositions = 0;
		for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() ){
			for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() ){
				Coord_orth gridCoord = iv.coord_orth();
				double dx = gridCoord[0] - atomCoordOrth[0];
				if (fabs(dx) <= vdwRadius){
					double dy = gridCoord[1] - atomCoordOrth[1];
					if (fabs(dy) <= vdwRadius){
						double dxdysq = dx*dx + dy*dy;
						if (dxdysq <= vdwRadiusSq){
							double dz = sqrt (vdwRadiusSq - dxdysq);
							
							Coord_orth candidate = gridCoord;
							candidate[2] = atomCoordOrth[2] + dz;
							Coord_orth delta = candidate - atomCoordOrth;
							for (int i=0; i<3; i++) delta[i] = delta[i] * extrapolate;
							candidate = atomCoordOrth + delta;
							int buried = coordIsBuriedByNeighbours(candidate, iAtom);
							if (!buried) {
								allowProbeToEatWithinGridRange(candidate, atomRange);
							}
							
							candidate = gridCoord;
							candidate[2] = atomCoordOrth[2] - dz;
							delta = candidate - atomCoordOrth;
							for (int i=0; i<3; i++) delta[i] = delta[i] * extrapolate;
							candidate = atomCoordOrth + delta;
							buried = coordIsBuriedByNeighbours(candidate, iAtom);
							if (!buried) {
								allowProbeToEatWithinGridRange(candidate, atomRange);
							}
						}
					}
				}
			}
		}
		for ( iv = i0; iv.coord().v() <= g1.v(); iv.next_v() ){
			for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ){
				Coord_orth gridCoord = iw.coord_orth();
				double dy = gridCoord[1] - atomCoordOrth[1];
				if (fabs(dy) <= vdwRadius){
					double dz = gridCoord[2] - atomCoordOrth[2];
					if (fabs(dz) <= vdwRadius){
						double dydzsq = dy*dy + dz*dz;
						if (dydzsq <= vdwRadiusSq){
							double dx = sqrt (vdwRadiusSq - dydzsq);
							
							Coord_orth candidate = gridCoord;
							candidate[0] = atomCoordOrth[0] + dx;
							Coord_orth delta = candidate - atomCoordOrth;
							for (int i=0; i<3; i++) delta[i] = delta[i] * extrapolate;
							candidate = atomCoordOrth + delta;
							int buried = coordIsBuriedByNeighbours(candidate, iAtom);
							if (!buried) {
								allowProbeToEatWithinGridRange(candidate, atomRange);
							}
							
							candidate = gridCoord;
							candidate[0] = atomCoordOrth[0] - dx;
							delta = candidate - atomCoordOrth;
							for (int i=0; i<3; i++) delta[i] = delta[i] * extrapolate;
							candidate = atomCoordOrth + delta;
							buried = coordIsBuriedByNeighbours(candidate, iAtom);
							if (!buried) {
								allowProbeToEatWithinGridRange(candidate, atomRange);
							}
						}
					}
				}
			}
			for ( iw = i0; iw.coord().w() <= g1.w(); iw.next_w() ){
				for ( iu = iw; iu.coord().u() <= g1.u(); iu.next_u() ){
					Coord_orth gridCoord = iu.coord_orth();
					double dz = gridCoord[2] - atomCoordOrth[2];
					if (fabs(dz) <= vdwRadius){
						double dx = gridCoord[0] - atomCoordOrth[0];
						if (fabs(dx) <= vdwRadius){
							double dzdxsq = dz*dz + dx*dx;
							if (dzdxsq <= vdwRadiusSq){
								double dy = sqrt (vdwRadiusSq - dzdxsq);
								
								Coord_orth candidate = gridCoord;
								candidate[1] = atomCoordOrth[1] + dy;
								Coord_orth delta = candidate - atomCoordOrth;
								for (int i=0; i<3; i++) delta[i] = delta[i] * extrapolate;
								candidate = atomCoordOrth + delta;
								int buried = coordIsBuriedByNeighbours(candidate, iAtom);
								if (!buried) {
									allowProbeToEatWithinGridRange(candidate, atomRange);
								}
								
								candidate = gridCoord;
								candidate[1] = atomCoordOrth[1] - dy;
								delta = candidate - atomCoordOrth;
								for (int i=0; i<3; i++) delta[i] = delta[i] * extrapolate;
								candidate = atomCoordOrth + delta;
								buried = coordIsBuriedByNeighbours(candidate, iAtom);
								if (!buried) {
									allowProbeToEatWithinGridRange(candidate, atomRange);
								}
							}
						}
					}
				}
			}
		}
	}
	return 0;
}

clipper::Xmap<double> &CXX_mot::CXXQADSurface::getDoubleMap(){
	return theDoubleMap;
}

clipper::Cell &CXX_mot::CXXQADSurface::getCell(){
	return clipperCell;
}

void CXX_mot::CXXQADSurface::dump(Grid_sampling theObject){
	cout << "Grid sampling :" << endl <<
	theObject.nu() << " " <<
	theObject.nv() << " " <<
	theObject.nw() << endl;
}

void CXX_mot::CXXQADSurface::dump(Coord_grid theObject){
	cout << "Coord grid :" << endl <<
	theObject.u() << " " <<
	theObject.v() << " " <<
	theObject.w() << endl;
}

void CXX_mot::CXXQADSurface::dump(Grid_range theObject){
	cout << "Grid limits :"  << endl <<
	theObject.min().u() << " " << theObject.max().u() << " " <<
	theObject.min().v() << " " << theObject.max().v() << " " <<
	theObject.min().w() << " " << theObject.max().w() << endl;
}

int CXX_mot::CXXQADSurface::allowProbesToEat(){
	Xmap<double>::Map_reference_coord i0, iu, iv, iw;
	
	for (unsigned iProbe = 0; iProbe < probePositions.size(); iProbe++){
		double eatingRadius = probeRadius+sample;
		double eatingRadiusSq = eatingRadius*eatingRadius;
		Grid_range gd (clipperCell, clipperGridSampling, eatingRadius);
		
		Coord_orth probeCoordOrth = probePositions[iProbe];
		
		Coord_frac uvw = probeCoordOrth.coord_frac(clipperCell);
		Coord_grid g0 = uvw.coord_grid(clipperGridSampling) + gd.min();
		Coord_grid g1 = uvw.coord_grid(clipperGridSampling) + gd.max();
		
		i0 = Xmap<double>::Map_reference_coord( theDoubleMap, g0 );
		
		for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() ){
			for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() ){
				for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ){
					if ( theFlagMap[iw] == Inaccessible){
						Coord_orth gridCoord = iw.coord_orth();
						double dx = gridCoord[0] - probeCoordOrth[0];
						if (fabs(dx) <= eatingRadius){
							double dy = gridCoord[1] - probeCoordOrth[1];
							if (fabs(dy) <= eatingRadius){
								double dxdysq = dx*dx + dy*dy;
								if (dxdysq <= eatingRadiusSq){
									double dz = gridCoord[2] - probeCoordOrth[2];
									if (fabs(dz) <= eatingRadius){
										double dxdydzsq = dz*dz+dxdysq;
										if (dxdydzsq <= eatingRadiusSq){
											double distanceFromProbeSq = dxdydzsq;
											if (distanceFromProbeSq <= theDoubleMap[iw]){
												theDoubleMap[iw] = distanceFromProbeSq;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return 0;
}

int CXX_mot::CXXQADSurface::allowProbeToEatWithinGridRange(Coord_orth probeCoordOrth, Grid_range theRange){
	Xmap<double>::Map_reference_coord i0, iu, iv, iw;
	
	double eatingRadius = probeRadius+sample;
	double eatingRadiusSq = eatingRadius*eatingRadius;
	Grid_range gd (clipperCell, clipperGridSampling, eatingRadius);
	
	Coord_frac uvw = probeCoordOrth.coord_frac(clipperCell);
	Coord_grid g0 = uvw.coord_grid(clipperGridSampling) + gd.min();
	Coord_grid g1 = uvw.coord_grid(clipperGridSampling) + gd.max();
	Grid_range probeRange(g0, g1);
	Grid_range theIntersection;
	
	gdIntersection(probeRange, theRange, theIntersection);
	
	g0 = theIntersection.min();
	g1 = theIntersection.max();
	
	i0 = Xmap<double>::Map_reference_coord( theDoubleMap, g0 );
	
	for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() ){
		for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() ){
			for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ){
				if ( theFlagMap[iw] == Inaccessible){
					Coord_orth gridCoord = iw.coord_orth();
					double dx = gridCoord[0] - probeCoordOrth[0];
					if (fabs(dx) <= eatingRadius){
						double dy = gridCoord[1] - probeCoordOrth[1];
						if (fabs(dy) <= eatingRadius){
							double dxdysq = dx*dx + dy*dy;
							if (dxdysq <= eatingRadiusSq){
								double dz = gridCoord[2] - probeCoordOrth[2];
								if (fabs(dz) <= eatingRadius){
									double dxdydzsq = dz*dz+dxdysq;
									if (dxdydzsq < eatingRadiusSq){
										double distanceFromProbeSq = dxdydzsq;
										if (distanceFromProbeSq < theDoubleMap[iw]){
											theDoubleMap[iw] = distanceFromProbeSq;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return 0;
}

int CXX_mot::CXXQADSurface::coordIsBuriedByNeighbours(Coord_orth &point,int iAtom1){
	mmdb::PAtom Atom1 = selectedAtoms[iAtom1];
	vector<int> &theNeighbourhood = neighbourhoods[iAtom1];
	int buried = 0;
	
	for (unsigned iAtom2 = 0; iAtom2<theNeighbourhood.size() && !buried; iAtom2++){
		mmdb::PAtom Atom2 = selectedAtoms[theNeighbourhood[iAtom2]];
		if (Atom2 != Atom1){
			Coord_orth atomCoordOrth(Atom2->x, Atom2->y, Atom2->z);
			double accessibleRadius  = fastGetAtomRadius(iAtom2)+probeRadius;
			double accessibleRadiusSq = accessibleRadius * accessibleRadius;
			double dx = point[0] - atomCoordOrth[0];
			if (fabs(dx) <= accessibleRadius){
				double dy = point[1] - atomCoordOrth[1];
				if (fabs(dy) <= accessibleRadius){
					double dxdysq = dx*dx + dy*dy;
					if (dxdysq <= accessibleRadiusSq){
						double dz = point[2] - atomCoordOrth[2];
						if (fabs(dz) <= accessibleRadius){
							double dxdydzsq = dz*dz+dxdysq;
							if (dxdydzsq < accessibleRadiusSq){
								buried = 1;
							}
						}
					}
				}
			}
		}
	}
	return buried;
}

void CXX_mot::CXXQADSurface::addProbe(Coord_orth newProbePosition){
	probePositions.push_back(newProbePosition);
}

void CXX_mot::CXXQADSurface::copyFlagToDouble(){
	clipper::Xmap_base::Map_reference_index ix;
	for ( ix = theFlagMap.first(); !ix.last(); ix.next() ) {
		theDoubleMap[ix] = theFlagMap[ix];
	}
}

int CXX_mot::CXXQADSurface::contourMap(double isoLevel){
	int perPixelTriangles[24][3];
	Coord_orth perPixelVertices[12];
	vector<vector<int> > referencingTriangles;
	
	clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw;
	i0 = clipper::Xmap_base::Map_reference_coord( theDoubleMap, clipperGridRange.min() );
	
	int compareStart = 0;
	int oldCompareStart  = 0;
	vector<clipper::Coord_orth> vertVector;
	vertVector.reserve(12);
	for ( iu = i0; iu.coord().u() <= clipperGridRange.max().u(); iu.next_u() ){
		for ( iv = iu; iv.coord().v() <= clipperGridRange.max().v(); iv.next_v() ){
			for ( iw = iv; iw.coord().w() <= clipperGridRange.max().w(); iw.next_w() ) {
				vertVector.resize(0);
				int ntriangs = contourPixel(iw, isoLevel, perPixelVertices, 
											vertVector, perPixelTriangles);
				for (int iTriang = 0; iTriang< ntriangs; iTriang++){
					for (int i=0; i<3; i++){
						//Compare this vertex with all other ones on this and previous layer:
						//Note that we can probably make this faster by searching backwards !
						int iEquivalent = 0;
						if (vertices.size() >=1){
							for (int iCompare = vertices.size()-1; 
								 iCompare >=compareStart && !iEquivalent; 
								 iCompare--){
								if (compareVertices(vertVector[perPixelTriangles[iTriang][i]],
													vertices[iCompare])){
									iEquivalent = iCompare;
								}
							}
						}
						//Check to see if this is a brand new vertex
						if (!iEquivalent) {
							iEquivalent = vertices.size();
							vertices.push_back(vertVector[perPixelTriangles[iTriang][i]]);
						}
						triangles.push_back(iEquivalent);
					}
				}
			}
		}
		compareStart = oldCompareStart;
		oldCompareStart = vertices.size();
	}
	return 0;
}

int CXX_mot::CXXQADSurface::contourPixel(Xmap_base::Map_reference_coord index, double isoLevel,
								Coord_orth vertlist[], vector<clipper::Coord_orth> &vertVector, 
								int triangles[][3]){
	//cubeindex is the index into thelookup table of what edges and what
	//triangles are implicit for a given voxel
	GRIDCELL grid;
	int cubeindex = 0;
	
	clipper::Xmap_base::Map_reference_coord localIndex = index;
	if (theDoubleMap[localIndex]<isoLevel) cubeindex |=1;
	grid.p[0] = localIndex.coord_orth();
	grid.val[0] = theDoubleMap[localIndex];
	
	localIndex = index; localIndex.next_u();
	if (theDoubleMap[localIndex]<isoLevel) cubeindex |=2;
	grid.p[1] = localIndex.coord_orth();
	grid.val[1] = theDoubleMap[localIndex];
	
	localIndex = index; localIndex.next_u(); localIndex.next_v();
	if (theDoubleMap[localIndex]<isoLevel) cubeindex |=4;
	grid.p[2] = localIndex.coord_orth();
	grid.val[2] = theDoubleMap[localIndex];
	
	localIndex = index; localIndex.next_v();
	if (theDoubleMap[localIndex]<isoLevel) cubeindex |=8;
	grid.p[3] = localIndex.coord_orth();
	grid.val[3] = theDoubleMap[localIndex];
	
	localIndex = index; localIndex.next_w();
	if (theDoubleMap[localIndex]<isoLevel) cubeindex |=16;
	grid.p[4] = localIndex.coord_orth();
	grid.val[4] = theDoubleMap[localIndex];
	
	localIndex = index; localIndex.next_u(); localIndex.next_w();
	if (theDoubleMap[localIndex]<isoLevel) cubeindex |=32;
	grid.p[5] = localIndex.coord_orth();
	grid.val[5] = theDoubleMap[localIndex];
	
	localIndex = index; localIndex.next_u(); localIndex.next_v(); localIndex.next_w();
	if (theDoubleMap[localIndex]<isoLevel) cubeindex |=64;
	grid.p[6] = localIndex.coord_orth();
	grid.val[6] = theDoubleMap[localIndex];
	
	localIndex = index; localIndex.next_v(); localIndex.next_w();
	if (theDoubleMap[localIndex]<isoLevel) cubeindex |=128;
	grid.p[7] = localIndex.coord_orth();
	grid.val[7] = theDoubleMap[localIndex];
	
	//cubeindex of zero implies all corners of the cube are below the isoLevel
	if (cubeindex == 0 || cubeindex == 255) return 0;
	
	if (edgeTable[cubeindex] & 1){
		vertlist[0] =
		coordOrthInterp(isoLevel,grid.p[0],grid.p[1],grid.val[0],grid.val[1]);
		vertVector.push_back(vertlist[0]);
	}
	
	if (edgeTable[cubeindex] & 2){
		vertlist[1] =
		coordOrthInterp(isoLevel,grid.p[1],grid.p[2],grid.val[1],grid.val[2]);
		vertVector.push_back(vertlist[1]);
	}
	
	if (edgeTable[cubeindex] & 4){
		vertlist[2] =
		coordOrthInterp(isoLevel,grid.p[2],grid.p[3],grid.val[2],grid.val[3]);
		vertVector.push_back(vertlist[2]);
	}
	
	if (edgeTable[cubeindex] & 8){
		vertlist[3] =
		coordOrthInterp(isoLevel,grid.p[3],grid.p[0],grid.val[3],grid.val[0]);
		vertVector.push_back(vertlist[3]);
	}
	
	if (edgeTable[cubeindex] & 16){
		vertlist[4] =
		coordOrthInterp(isoLevel,grid.p[4],grid.p[5],grid.val[4],grid.val[5]);
		vertVector.push_back(vertlist[4]);
	}
	
	if (edgeTable[cubeindex] & 32){
		vertlist[5] =
		coordOrthInterp(isoLevel,grid.p[5],grid.p[6],grid.val[5],grid.val[6]);
		vertVector.push_back(vertlist[5]);
	}
	
	if (edgeTable[cubeindex] & 64){
		vertlist[6] =
		coordOrthInterp(isoLevel,grid.p[6],grid.p[7],grid.val[6],grid.val[7]);
		vertVector.push_back(vertlist[6]);
	}
	
	if (edgeTable[cubeindex] & 128){
		vertlist[7] =
		coordOrthInterp(isoLevel,grid.p[7],grid.p[4],grid.val[7],grid.val[4]);
		vertVector.push_back(vertlist[7]);
	}
	
	if (edgeTable[cubeindex] & 256){
		vertlist[8] =
		coordOrthInterp(isoLevel,grid.p[0],grid.p[4],grid.val[0],grid.val[4]);
		vertVector.push_back(vertlist[8]);
	}
	
	if (edgeTable[cubeindex] & 512){
		vertlist[9] =
		coordOrthInterp(isoLevel,grid.p[1],grid.p[5],grid.val[1],grid.val[5]);
		vertVector.push_back(vertlist[9]);
	}
	
	if (edgeTable[cubeindex] & 1024){
		vertlist[10] =
		coordOrthInterp(isoLevel,grid.p[2],grid.p[6],grid.val[2],grid.val[6]);
		vertVector.push_back(vertlist[10]);
	}
	
	if (edgeTable[cubeindex] & 2048){
		vertlist[11] =
		coordOrthInterp(isoLevel,grid.p[3],grid.p[7],grid.val[3],grid.val[7]);
		vertVector.push_back(vertlist[11]);
	}
	
	int ntriang = 0;
	for (int i=0;triTable[cubeindex][i]!=-1;i+=3) {
		triangles[ntriang][0] = triTablePrime[cubeindex][i  ];
		triangles[ntriang][1] = triTablePrime[cubeindex][i+1];
		triangles[ntriang][2] = triTablePrime[cubeindex][i+2];
		ntriang++;
	}
	return ntriang;
}	

Coord_orth CXX_mot::CXXQADSurface::coordOrthInterp(double isoLevel, Coord_orth p1, Coord_orth p2, double v1, double v2){
	Coord_orth result;
	if (fabs(v1-isoLevel) < 1e-08) return p1;
	if (fabs(v2-isoLevel) < 1e-08) return p2;
	if (fabs(v2-v1) < 1e-08) return p1;
	double fraction = (isoLevel-v1) / (v2-v1);
	Coord_orth diffvec = p2 - p1;
	for (int i=0; i<3; i++){
		result[i] = p1[i] + fraction * diffvec[i];
	}
	return result;
}

int CXX_mot::CXXQADSurface::compareVertices(Coord_orth &v1, Coord_orth &v2){
	return (v1.x() == v2.x() && v1.y() == v2.y() && v1.z() == v2.z());
}

vector<clipper::Coord_orth> &CXX_mot::CXXQADSurface::getVertices(){
	return vertices;
}
vector<clipper::Coord_orth> &CXX_mot::CXXQADSurface::getNormals(){
	return vertexNormals;
}
vector<int> &CXX_mot::CXXQADSurface::getTriangles(){
	return triangles;
}

int CXX_mot::CXXQADSurface::calculateAveragedNormals(){
	//Now calculate normals
	vertexNormals.resize(vertices.size());
	for (unsigned i=0; i<vertices.size(); i++){
		vertexNormals[i] = Coord_orth(0., 0., 0.);
	}
	//Calculate a per triangle normal that gets added into the vertex normals of
	//each of the participant vertices
	for (unsigned iTriang=0; iTriang < triangles.size()/3; iTriang++){
		int iv0 = triangles[3*iTriang];
		int iv1 = triangles[3*iTriang+1];
		int iv2 = triangles[3*iTriang+2];
		Coord_orth v01 = vertices[iv1] - vertices[iv0];
		Coord_orth v02 = vertices[iv2] - vertices[iv0];
		Coord_orth normal (  v01[1]*v02[2]-v02[1]*v01[2],
							 -(v01[0]*v02[2]-v02[0]*v01[2]),
							 v01[0]*v02[1]-v02[0]*v01[1] );
		double nLengthSq = normal.lengthsq();
		
		//The following excludes zero area triangles from the normal averaging procedure
		//I remain unclear how I can end up with zero area triangles however 
		if (nLengthSq > 1e-14){
			for (int i=0; i<3; i++) normal[i] /= nLengthSq;
			vertexNormals[iv0] = vertexNormals[iv0] + normal;
			vertexNormals[iv1] = vertexNormals[iv1] + normal;
			vertexNormals[iv2] = vertexNormals[iv2] + normal;
		}
	}
	//Then convert to per vertex normal
	for (unsigned iVertex = 0; iVertex < vertices.size(); iVertex++){
		double nLength = sqrt(vertexNormals[iVertex].lengthsq());
		if (nLength  < 1e-12) cout << "On dear \n";
		for (int i=0; i<3; i++) vertexNormals[iVertex][i] = vertexNormals[iVertex][i] / nLength;
	}
	return 0;
	
}

double CXX_mot::CXXQADSurface::getAtomRadius(mmdb::PAtom theAtom){
	//Here get handle of a radius data type from MMDB if such has been stored
	int iRadiusHandle = theMMDBManager->GetUDDHandle(mmdb::UDR_ATOM, "PerAtomRadius");
	double theRadius;
	if (iRadiusHandle>0){
		int success = theAtom->GetUDData (iRadiusHandle, theRadius);
		if (success != mmdb::UDDATA_Ok) theRadius = mmdb::getVdWaalsRadius(theAtom->element);
	}
	else theRadius = mmdb::getVdWaalsRadius(theAtom->element);
	return theRadius;
}

int CXX_mot::CXXQADSurface::transformTriTable(){
	cout << "{ ";
	for (int iTriang = 0; iTriang< 256; iTriang++){
		cout << "{ ";
		vector<int> oldIDs;
		map<int, int> correspondence;
		for (int iVertex = 0; iVertex < 16 && triTable[iTriang][iVertex] != -1; iVertex++){
			oldIDs.push_back(triTable[iTriang][iVertex]);
		}
		sort(oldIDs.begin(), oldIDs.end());
		for (unsigned iID=0; iID<oldIDs.size(); iID++){
			if ( correspondence.find(oldIDs[iID]) == correspondence.end()){
				correspondence[oldIDs[iID]] = correspondence.size()-1;
			}
		}
		for (int iVertex = 0; iVertex < 16; iVertex++){
			if (triTable[iTriang][iVertex] == -1) {
				cout << " -1,";
			}
			else {
				cout << " " << correspondence[triTable[iTriang][iVertex]] << ",";
			}
		}
		cout << "}" << endl;
	}
	cout << "}";
	return 0;
}

int CXX_mot::CXXQADSurface::setInaccessibleDistanceSq(){
	//Here adopt Solvent accessible surface contour vertices into the probe list
	clipper::Xmap_base::Map_reference_index ix;
	for (ix = theDoubleMap.first(); !ix.last(); ix.next() ) {
		if (theFlagMap[ix] == Inaccessible || theFlagMap[ix] == Accessible){
			theDoubleMap[ix] = (probeRadius+sample) * (probeRadius+sample);
			theFlagMap[ix] = Inaccessible;
		}
	}
	return 0;
}

int CXX_mot::CXXQADSurface::sqrtDistanceSq(){
	//Here convert from DistanceSqmap to Distance Map for contouring
	clipper::Xmap_base::Map_reference_index ix;
	for (ix = theDoubleMap.first(); !ix.last(); ix.next() ) {
		if (theDoubleMap[ix]!= Solvent){
			theDoubleMap[ix] = sqrt(theDoubleMap[ix]);
		}
	}
	return 0;
}

int CXX_mot::CXXQADSurface::toruses()
{
	int nTorusSegments = 0;
	Xmap<double>::Map_reference_coord i0, iu, iv, iw;
	CXXCoord xAxis(1.,0.,0.);
	CXXCoord yAxis(0.,1.,0.);
	CXXCoord zAxis(0.,0.,1.);
			
	//Identify a big box around this torus, and check if all the pixels inside
	//the box are within the allowed part of the torus
	double eatingRadius = probeRadius + sample;
	double eatingRadiusSq = eatingRadius * eatingRadius;

	//loop over all atoms in selection
	for (int atomNr  = 0;atomNr < nSelectedAtoms; atomNr++) { 
		mmdb::PAtom centralAtom = selectedAtoms[atomNr];

		Coord_orth atomCoordOrth(centralAtom->x, centralAtom->y, centralAtom->z);
		Coord_frac uvw = atomCoordOrth.coord_frac(clipperCell);
		Grid_range gd (clipperCell, clipperGridSampling, getAtomRadius(centralAtom)+probeRadius);
		Coord_grid g0 = uvw.coord_grid(clipperGridSampling) + gd.min();
		Coord_grid g1 = uvw.coord_grid(clipperGridSampling) + gd.max();
		Grid_range gdAtom(g0, g1);
		
		// now add central atom to the hood
		double radiusOfAtom1 = fastGetAtomRadius(atomNr);
		CXXNewHood theNewHood(centralAtom, radiusOfAtom1, probeRadius);
		
		for (unsigned sphereAtomNr = 0; sphereAtomNr < neighbourhoods[atomNr].size(); sphereAtomNr++) {
			mmdb::PAtom sphereAtom = selectedAtoms[neighbourhoods[atomNr][sphereAtomNr]];
			double radiusOfAtom2 = fastGetAtomRadius(sphereAtomNr);
			theNewHood.addAtom(sphereAtom, radiusOfAtom2);
		}
		
		//Find the non-hidden segments of the circles
		theNewHood.findSegments();
		
		//Loop over the segments that are left, using them to "eat" into the molecular volume
		std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> > &theCircles(theNewHood.getCircles());
		std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::iterator circlesEnd = theCircles.end();
		for (std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::iterator circleIter = theCircles.begin();
			 circleIter != circlesEnd;
			 ++circleIter){
			CXXCircle &theCircle(*circleIter);
			
			if (centralAtom->serNum < theCircle.getAtomJ()->serNum){
				if (theCircle.nSegments()>0){
					const CXXCoord &torusCentre(theCircle.getCentreOfCircle());
					const CXXCoord &torusAxis(theCircle.getNormal());
					double rTraj = theCircle.getRadiusOfCircle();	
					const CXXCoord &v1unit(theCircle.getReferenceUnitRadius());
					
					double rMax = rTraj + probeRadius;
					double rMin = rTraj - probeRadius;
					double rMaxSq = rMax*rMax;
					double rMinSq = rMin*rMin;

					double xDotAxis = xAxis*torusAxis;
					double yDotAxis = yAxis*torusAxis;
					double zDotAxis = zAxis*torusAxis;
					
					double xCoeff = sqrt (1.0 - xDotAxis*xDotAxis);
					double yCoeff = sqrt (1.0 - yDotAxis*yDotAxis);
					double zCoeff = sqrt (1.0 - zDotAxis*zDotAxis);
					
					Coord_orth minOrth(torusCentre[0] - (xCoeff*rTraj + probeRadius),
									   torusCentre[1] - (yCoeff*rTraj + probeRadius) ,
									   torusCentre[2] - (zCoeff*rTraj + probeRadius));
					for (int i=0; i<3; i++) if (minOrth[i]<0.) minOrth[i] -=sample;
					
					Coord_orth maxOrth(torusCentre[0] + (xCoeff*rTraj + probeRadius),
									   torusCentre[1] + (yCoeff*rTraj + probeRadius) ,
									   torusCentre[2] + (zCoeff*rTraj + probeRadius));
					for (int i=0; i<3; i++) if (maxOrth[i]>0.) maxOrth[i] +=sample;
					
					g0 = minOrth.coord_frac(clipperCell).coord_grid(clipperGridSampling);
					g1 = maxOrth.coord_frac(clipperCell).coord_grid(clipperGridSampling);				

					Grid_range gdTorus (g0, g1);
					Grid_range grCombined; 
					
					gdIntersection(gdTorus, gdAtom, grCombined);
					g0 = grCombined.min();
					g1 = grCombined.max();
					
					i0 = Xmap<double>::Map_reference_coord( theDoubleMap, g0 );
					
					for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() ){
						for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() ){
							for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ){
								if ( theFlagMap[iw] == Inaccessible){
									Coord_orth gridCoord = iw.coord_orth();
									CXXCoord gridCXXCoord(gridCoord[0], gridCoord[1], gridCoord[2]);
									
									CXXCoord centreToGrid = gridCXXCoord - torusCentre;
									double rParallel= centreToGrid*torusAxis;
									if (fabs(rParallel) < eatingRadius){
										//Project the grid point bac into the plane of the circle
										CXXCoord parallelOffset = torusAxis;
										parallelOffset *= rParallel;
										CXXCoord projected = gridCXXCoord - parallelOffset;
										CXXCoord radiusProjected = projected - torusCentre;
										double rPerpSq = radiusProjected.get3DLengthSq();
										if (rPerpSq < rMaxSq && rPerpSq > rMinSq){
											double rSection = sqrt(eatingRadiusSq - rParallel*rParallel);
											double rPerpMax = rTraj + rSection;
											double rPerpMaxSq = rPerpMax*rPerpMax;
											double rPerpMin = rTraj - rSection;
											double rPerpMinSq = rPerpMin * rPerpMin;
											if (rPerpSq > rPerpMinSq &&
												rPerpSq < rPerpMaxSq){
												double rPerp = sqrt(rPerpSq);
												radiusProjected /= rPerp;
												//The last thing we have to evaluate is how far around the
												//torus trajectory we are, and whether this is within the allowed range
												radiusProjected *= rTraj;
												CXXCoord probePosition = torusCentre + radiusProjected;
												CXXCoord offset = gridCXXCoord - probePosition;
												double distanceFromProbeSq = offset.get3DLengthSq();
												if (distanceFromProbeSq < theDoubleMap[iw]){
													int inRange = 0;
													for (unsigned iEdge = 0; 
														 iEdge < theCircle.nSegments() && !inRange; 
														 iEdge++){
														
														double thetaMin = theCircle.start(iEdge)->getAngle();
														double thetaMax = theCircle.stop(iEdge)->getAngle();
														if (thetaMax == 2.*M_PI){
															inRange = 1;
														}
														else {
															double theta = torusAxis.angleBetween(v1unit, radiusProjected);
															if (theta>=thetaMin && theta <= thetaMax) inRange = 1;
														}
													}
													if (inRange){
														//Here if we have established we are inside the donut !
														theDoubleMap[iw] = distanceFromProbeSq;
													}
												}
											}
										}
									}
								}
							}
						}
					}							
				}
			}
		}
		//Now collect a list of points where the probes are in contact with three atoms
		circlesEnd = theCircles.end();
		for (std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::iterator circleIter = theCircles.begin();
			 circleIter != circlesEnd;
			 ++circleIter){
			CXXCircle &theCircle(*circleIter);
			if (!theCircle.getEaten()){
                const list<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> > &nodes = theCircle.getNodes();
                list<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> >::const_iterator endNode = nodes.end();
                for (list<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> >::const_iterator nodeIter = nodes.begin();
                     nodeIter != endNode;
                     ++nodeIter){
					const CXXCircleNode &aNode(*nodeIter);
					if (!aNode.isDeleted()){
						if (aNode.getAtomK()->serNum > aNode.getAtomJ()->serNum &&
							aNode.getAtomJ()->serNum > aNode.getAtomI()->serNum){
							CXXCoord nodeCXXCoord(aNode.getCoord());
							Coord_orth newProbe(nodeCXXCoord[0], nodeCXXCoord[1], nodeCXXCoord[2]);
							allowProbeToEatWithinGridRange(newProbe, gdAtom);
						}
					}
				}
			}
		}			
	}
	cout << "Number of Torus segments was " << nTorusSegments << " ";
	return 0;	
}

Grid_range &CXX_mot::CXXQADSurface::gdIntersection (Grid_range &g0, Grid_range &g1, Grid_range &theResult){
	Coord_grid newMin(max(g0.min()[0], g1.min()[0]),
					  max(g0.min()[1], g1.min()[1]),
					  max(g0.min()[2], g1.min()[2]));
	Coord_grid newMax(min(g0.max()[0], g1.max()[0]),
					  min(g0.max()[1], g1.max()[1]),
					  min(g0.max()[2], g1.max()[2]));
	
	theResult=Grid_range(newMin, newMax);

	return theResult;
}


