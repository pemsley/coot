/*
 *  CXXBall.cpp
 *  CXXSurface
 *
 *  Created by Martin Noble on 10/06/2009.
 *  Copyright 2009 LMB, Oxford University. All rights reserved.
 *
 */

#include "CXXBall.h"
#include "CXXNewHood.h"
#include "CXXSurface.h"

CXX_old::CXXAlloc<CXX_mot::CXXReentrantProbeBall> CXX_mot::CXXReentrantProbeBall::allocator = CXX_old::CXXAlloc<CXX_mot::CXXReentrantProbeBall>();

int CXX_mot::CXXBall::triangulateBalls(vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall*> > &ballPntrs,
                              vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall*> > &contextBallPntrs,
							  double delta, CXXSurface *aSurface, int insideOrOutside)
{
    std::map<const CXXBall *, std::vector<const CXXBall *, CXX_old::CXXAlloc<const CXXBall *> > >contactMap;
    ballContacts(ballPntrs, contextBallPntrs, contactMap);	
    std::cout << "Established contact map\n";
    
    //Now pass through our reentrant Probe list, triangulating them
    int nBalls = ballPntrs.size();
	
	//It goes like this...  
	// Each sphere as it intersects with other spheres will generate loose ends...
	// We will store those loose ends and attribute them
	// in an stl::map to the spheres with which the intersection occurs
	std::vector<std::map<const CXXBall *, std::vector<CXXCoord, CXX_old::CXXAlloc<CXXCoord> > > > raggedEdges;
	raggedEdges.resize(nBalls);
#ifndef NEED_OPENMP_PRAGMA_HACK	
#pragma omp parallel for default(none) shared(delta, nBalls, contactMap, aSurface, raggedEdges, insideOrOutside, ballPntrs) schedule(dynamic, 10) //num_threads(2)
#else	
#pragma omp parallel for default(none) shared(delta, nBalls, contactMap, aSurface, raggedEdges, insideOrOutside) schedule(dynamic, 10) //num_threads(2)
#endif	
    for (int i=0; i< nBalls; i++){
        const CXXBall &ball(*ballPntrs[i]);
        CXXNewHood ballHood;
        ballHood.initWith(&ball);
        std::map<const CXXBall *, std::vector<const CXXBall *, CXX_old::CXXAlloc<const CXXBall *> > >::iterator contacts;
        contacts = contactMap.find(ballPntrs[i]);
        if (contacts!=contactMap.end()){
            std::vector<const CXXBall *, CXX_old::CXXAlloc<const CXXBall *> >::iterator neighboursEnd(contacts->second.end());
            for (std::vector<const CXXBall *, CXX_old::CXXAlloc<const CXXBall *> >::iterator neighbour = contacts->second.begin();
                 neighbour!=neighboursEnd;
                 ++neighbour){
                ballHood.addBall(**neighbour);
            }
        }	
        ballHood.triangulateAsBallHoodInto(aSurface, delta, raggedEdges[i], false, insideOrOutside);
    }
	
	//Reformat the ragged Edges, so that we have them separated such that all of the ragged edges
	//that need to be added to a particular probe are in an appropriate map associated with that probe	
	std::map<const CXXBall*, std::map<const CXXBall *, std::vector<CXXCoord, CXX_old::CXXAlloc<CXXCoord> > > >reformattedEdges;
	for (unsigned int i=0; i<raggedEdges.size(); i++){
		std::map<const CXXBall *, std::vector<CXXCoord, CXX_old::CXXAlloc<CXXCoord> > >::iterator raggedEdgesEnd = raggedEdges[i].end();
		for (std::map<const CXXBall *, std::vector<CXXCoord, CXX_old::CXXAlloc<CXXCoord> > >::iterator raggedEdgeSet = raggedEdges[i].begin();
			 raggedEdgeSet != raggedEdgesEnd; 
			 ++raggedEdgeSet){
			if (contactMap.find(raggedEdgeSet->first)!=contactMap.end()){
				std::map<const CXXBall *, std::vector<CXXCoord, CXX_old::CXXAlloc<CXXCoord> > > &ballsEdges = reformattedEdges[raggedEdgeSet->first];
				ballsEdges[ballPntrs[i]] = raggedEdgeSet->second;
			}
		}
	}
	raggedEdges.clear();
    
	int nToRedraw = reformattedEdges.size();
	std::cout << "There are " << nToRedraw << " trimmed balls\n";
    
	std::vector<const CXXBall *, CXX_old::CXXAlloc<const CXXBall*> >trimmedBalls;
	//Copy this list of balls into a vector to allow subsequent OpenMP parallelisation
	std::map<const CXXBall*, std::map<const CXXBall*, std::vector<CXXCoord, CXX_old::CXXAlloc<CXXCoord> > > >::iterator reformattedEdgeEnd = reformattedEdges.end();
	for (std::map<const CXXBall*, std::map<const CXXBall *, std::vector<CXXCoord, CXX_old::CXXAlloc<CXXCoord> > > >::iterator reformattedEdge = reformattedEdges.begin();
		 reformattedEdge != reformattedEdgeEnd; 
		 reformattedEdge++){
		trimmedBalls.push_back(reformattedEdge->first);
	};
	
#pragma omp parallel for default(none) shared(trimmedBalls, contactMap, insideOrOutside, reformattedEdges, delta, aSurface) schedule(dynamic, 10) //num_threads(2)
	for (unsigned int i=0; i<trimmedBalls.size(); i++){
		const CXXBall &ball(*(trimmedBalls[i]));
		std::map<const CXXBall *, std::vector<const CXXBall *, CXX_old::CXXAlloc<const CXXBall *> > >::iterator contacts;
		contacts = contactMap.find(&ball);
		CXXNewHood ballHood;
		ballHood.initWith(&ball);
		if (contacts!=contactMap.end()){
			std::vector<const CXXBall *, CXX_old::CXXAlloc<const CXXBall *> >::iterator neighboursEnd(contacts->second.end());
			for (std::vector<const CXXBall *, CXX_old::CXXAlloc<const CXXBall *> >::iterator neighbour = contacts->second.begin();
				 neighbour!=neighboursEnd;
				 ++neighbour){
				ballHood.addBall(**neighbour);
			}
		}
		ballHood.triangulateAsBallHoodInto(aSurface, delta, (reformattedEdges[&ball]), 
										   true, insideOrOutside);
	}
	return 0;
}

int CXX_mot::CXXBall::ballContacts(std::vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall*> > &balls, 
                          std::vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall*> > &contextBalls, 
						  std::map<const CXXBall*, std::vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall*> > > &contactMap) 
{
	int maxNBins = 20;
	if (balls.size() == 0) return 1;
	
	double maxBallRadius = -1.e30;
	std::vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall*> >::iterator ballsEnd = balls.end();
	for (std::vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall*> >::iterator ball = balls.begin();
		 ball != ballsEnd;
		 ++ball){
		maxBallRadius = (maxBallRadius > (*ball)->getRadius() ? maxBallRadius : (*ball)->getRadius() );
	}
	std::vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall*> >::iterator contextBallsEnd = contextBalls.end();
	for (std::vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall*> >::iterator ball = contextBalls.begin();
		 ball != contextBallsEnd;
		 ++ball){
		maxBallRadius = (maxBallRadius > (*ball)->getRadius() ? maxBallRadius : (*ball)->getRadius() );
	}
	std::cout << "Maximum radius was " << maxBallRadius << std::endl;
	
	double ballRadiusX2 = 2.*maxBallRadius;
	
	double limits[3][2];
	for (int i=0; i<3; i++){
		limits[i][0] = 1e30;
		limits[i][1] = -1e30;
	}
    //Establish limits of volume contaiing context Balls
	for (std::vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall*> >::iterator ballIter = contextBalls.begin(); 
		 ballIter!=contextBallsEnd; 
		 ++ballIter){
		for (int i=0; i<3; i++){
			limits[i][0] = (limits[i][0]<(**ballIter)[i]?limits[i][0]:(**ballIter)[i]);
			limits[i][1] = (limits[i][1]>(**ballIter)[i]?limits[i][1]:(**ballIter)[i]);
		}
	}
	for (int i=0; i<3; i++){
		limits[i][0] -= 0.0001;
		limits[i][1] += 0.0001;
	}
	
	double minimumBinSize = 2.*maxBallRadius;
	int nBins[3];
	double binWidth[3];
	
	for (int i=0; i<3; i++){
		int nMinimumBins = (int)ceil((limits[i][1]-limits[i][0])/minimumBinSize);
		if (nMinimumBins<maxNBins){
			nBins[i] = nMinimumBins;
			binWidth[i] = minimumBinSize;
		}
		else {
			nBins[i] = maxNBins;
			binWidth[i] = (limits[i][1]-limits[i][0])/double(maxNBins);
		}
	}
    std::cout << "Bins in x,y,z: "<<nBins[0]<<" "<<nBins[1]<<" "<<nBins[2]<<std::endl;
    std::cout << "Size in x,y,z: "<<binWidth[0]<<" "<<binWidth[1]<<" "<<binWidth[2]<<std::endl;
	//Prepare the bin vectors
	std::vector<std::vector<std::vector<std::vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall *> > > > >binnedballs;
	binnedballs.resize(nBins[0]);
	for (int i=0; i<nBins[0]; i++){
		binnedballs[i].resize(nBins[1]);
		for (int j=0; j<nBins[1]; j++){
			binnedballs[i][j].resize(nBins[2]);
		}
	}
	
	//Distribute balls among bins
	for (std::vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall*> >::iterator ballIter = contextBalls.begin(); 
         ballIter!=contextBallsEnd; 
         ++ballIter){
		int iBin[3];
		for (int i=0; i<3; i++){
			iBin[i] = (int)floor(((**ballIter)[i]-limits[i][0]) / binWidth[i]);
			iBin[i]  = max(0,iBin[i]);
			iBin[i] = min(nBins[i]-1,iBin[i]);
		}
		binnedballs[iBin[0]][iBin[1]][iBin[2]].push_back(*ballIter);
	}
    
    //To allow subsequent parallelization, create map entries in contactMap for each ball
	for (std::vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall*> >::iterator ballIter = balls.begin(); 
         ballIter!=ballsEnd; 
         ++ballIter){
		const CXXBall *binBall = *ballIter;
        contactMap[binBall].resize(0);
    }
	std::cout << "Ready for parallel contact search \n";
	std::cout.flush();
#ifndef NEED_OPENMP_PRAGMA_HACK	
#pragma omp parallel for default(none) shared (binnedballs, ballRadiusX2, limits, binWidth, nBins, balls, contactMap) schedule(dynamic,10) //num_threads(2)
#else	
#pragma omp parallel for default(none) shared (binnedballs, ballRadiusX2, limits, binWidth, nBins) schedule(dynamic,10) //num_threads(2)
#endif	
    for (unsigned int iBall=0; iBall< balls.size(); iBall++){
		int iBin[3];
		for (int i=0; i<3; i++){
			iBin[i] = (int)floor(((*balls[iBall])[i]-limits[i][0]) / binWidth[i]);
		}
		int startBinX = max(0,iBin[0]-1);
		int endBinX = min(nBins[0],iBin[0]+2);
		int startBinY = max(0,iBin[1]-1);
		int endBinY = min(nBins[1],iBin[1]+2);
		int startBinZ = max(0,iBin[2]-1);
		int endBinZ = min(nBins[2],iBin[2]+2);
        
        const CXXBall &centralballRef(*balls[iBall]);
        
        //Compare this in turn with all balls in this or neighbouring bins
        for (int searchBinX = startBinX; searchBinX<endBinX; searchBinX++){
            for (int searchBinY = startBinY; searchBinY<endBinY; searchBinY++){
                for (int searchBinZ = startBinZ; searchBinZ<endBinZ; searchBinZ++){
                    std::vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall*> > &searchBin(binnedballs[searchBinX][searchBinY][searchBinZ]);
                    std::vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall*> >::iterator binSearchEnd = searchBin.end();
                    for (std::vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall*> >::iterator otherball = searchBin.begin(); 
                         otherball!=binSearchEnd; 
                         ++otherball){
                        const CXXBall &otherballRef(**otherball);
                        if (&centralballRef != &otherballRef){
                            if (!centralballRef.getCoord().isNearly(otherballRef.getCoord(), 0.00001)){
                                if (centralballRef.getCoord().isNearly(otherballRef.getCoord(), ballRadiusX2)){
                                    CXXCoord diff(centralballRef.getCoord() - otherballRef.getCoord());
                                    if (diff.get3DLengthSq()<pow(centralballRef.getRadius() + otherballRef.getRadius(), 2.0)){
                                        contactMap[&centralballRef].push_back(&otherballRef);
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

double CXX_mot::CXXAtomBall::mostRecentDelta = 1.;

CXX_mot::CXXSphereElement CXX_mot::CXXAtomBall::unitSphereAtOrigin =
   CXX_mot::CXXSphereElement(CXX_mot::CXXCoord(0.,0.,0.), 1., 1.);
