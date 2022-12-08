/*
 *  CXXBall.cpp
 *  CXXSurface
 *
 *  Created by Martin Noble on 10/06/2009.
 *  Copyright 2009 LMB, Oxford University. All rights reserved.
 *
 */

#include "CXXSurfaceVertex.h"
#include "CXXBall.h"
#include "CXXNewHood.h"
#include "CXXSurface.h"

#if defined _OPENMP
#include <omp.h>
#else
#if __APPLE__
#include <dispatch/dispatch.h>
#endif
#endif

void CXXBall::initSphereElement(CXXSphereElement &theSphere, const double &delta, const CXXSphereElement &unitCellAtOriginForDelta) const{
    theSphere = unitCellAtOriginForDelta;
    theSphere.scaleBy(theRadius);
    theSphere.translateBy (theCoord);
}

void CXXAtomBall::initSphereElement(CXXSphereElement &theSphere, const double &delta, const CXXSphereElement &unitSphereAtOriginForDelta) const{
    theSphere = unitSphereAtOriginForDelta;
    theSphere.scaleBy(theRadius);
    theSphere.translateBy (theCoord);
    theSphere.setAtom(theAtom);
}


void CXXReentrantProbeBall::initSphereElement(CXXSphereElement &theSphere, const double &delta, const CXXSphereElement &unitCellAtOriginForDelta) const {
    theSphere.initWith(theCoord, theAtomI, theAtomJ, theAtomK,
                       delta, theRadius, includeAtoms);
    
}


int CXXBall::triangulateBalls(vector<const CXXBall*  > &ballPntrs,
                              vector<const CXXBall*  > &contextBallPntrs,
                              double delta, CXXSurfaceMaker *aSurface, int insideOrOutside)
{
    std::map<const CXXBall *, std::vector<const CXXBall *  > >contactMap;
    CXXBall::ballContacts(ballPntrs, contextBallPntrs, contactMap);
    //std::cout << "Established contact map\n";
    
    //Now pass through our reentrant Probe list, triangulating them
    size_t nBalls = ballPntrs.size();
    
    //It goes like this...
    // Each sphere as it intersects with other spheres will generate loose ends...
    // We will store those loose ends and attribute them
    // in an stl::map to the spheres with which the intersection occurs
    std::map<const CXXBall *, std::vector<CXXCoord<CXXCoord_ftype> > > * raggedEdges = new std::map<const CXXBall *, std::vector<CXXCoord<CXXCoord_ftype> > >[nBalls];
    //Prepare a unit sphere at origin for this delta
    const CXXSphereElement unitCellAtOriginForDelta = CXXSphereElement(CXXCoord<CXXCoord_ftype>(0.,0.,0.), 1., delta);
    
    size_t oldSize = aSurface->getChildSurfaces().size();
    aSurface->getChildSurfaces().resize(oldSize+nBalls);
    CXXSurface *ballSurfacesArray = &(aSurface->getChildSurfaces()[oldSize]);
    
#if __APPLE__ && !defined _OPENMP
    dispatch_apply(nBalls, dispatch_get_global_queue(0, 0), ^(size_t i){
#elif defined _OPENMP
          // #pragma omp parallel for default(none) shared(delta, nBalls, aSurface, raggedEdges, insideOrOutside, contactMap, ballSurfacesArray, ballPntrs, cout) schedule(dynamic, 10)
#pragma omp parallel for default(none) shared(delta, nBalls, aSurface, raggedEdges, insideOrOutside, contactMap, ballSurfacesArray, ballPntrs, unitCellAtOriginForDelta, cout) schedule(dynamic, 10)
#warning Compiling for OMP
        for (int i=0; i< nBalls; i++){
#else
        for (int i=0; i< nBalls; i++){
#endif
        CXXSurface &ballSurface = ballSurfacesArray[i];
        if (!(i%100) || i==nBalls-1) {
#if __APPLE__ && !defined _OPENMP
            dispatch_sync(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_LOW, NULL),^(void){
                cout << "Dealing with ball number " <<i <<endl;
            });
#elif defined _OPENMP
#pragma omp critical (cout)
	    cout << "Dealing with ball number " <<i << endl;
#endif
        }
        
        
        const CXXBall &ball(*ballPntrs[i]);
        CXXNewHood ballHood;
        ballHood.initWith(&ball);
        std::map<const CXXBall *, std::vector<const CXXBall *  > >::const_iterator contacts;
        contacts = contactMap.find(ballPntrs[i]);
        if (contacts!=contactMap.end()){
            std::vector<const CXXBall *  >::const_iterator neighboursEnd(contacts->second.end());
            for (std::vector<const CXXBall *  >::const_iterator neighbour = contacts->second.begin();
                 neighbour!=neighboursEnd;
                 ++neighbour){
                ballHood.addBall(**neighbour);
            }
        }
        std::map<const CXXBall *, std::vector<CXXCoord<CXXCoord_ftype> > >&raggedEdgesOfBall = raggedEdges[i];
        ballHood.triangulateAsBallHoodInto(ballSurface, delta, raggedEdgesOfBall, false, insideOrOutside, unitCellAtOriginForDelta);
#if defined TARGET_OS_MAC && !defined _OPENMP
    });
#else
    }
#endif
    
    //Reformat the ragged Edges, so that we have them separated such that all of the ragged edges
    //that need to be added to a particular probe are in an appropriate map associated with that probe
    std::map<const CXXBall*, std::map<const CXXBall *, std::vector<CXXCoord<CXXCoord_ftype> > > >reformattedEdges;
    for (int i=0; i<nBalls; i++){
        std::map<const CXXBall *, std::vector<CXXCoord<CXXCoord_ftype> > >::iterator raggedEdgesEnd = raggedEdges[i].end();
        for (std::map<const CXXBall *, std::vector<CXXCoord<CXXCoord_ftype> > >::iterator raggedEdgeSet = raggedEdges[i].begin();
             raggedEdgeSet != raggedEdgesEnd;
             ++raggedEdgeSet){
            if (contactMap.find(raggedEdgeSet->first)!=contactMap.end()){
                std::map<const CXXBall *, std::vector<CXXCoord<CXXCoord_ftype> > > &ballsEdges = reformattedEdges[raggedEdgeSet->first];
                ballsEdges[ballPntrs[i]] = raggedEdgeSet->second;
            }
        }
    }
    delete [] raggedEdges;
    
    std::vector<const CXXBall *  >trimmedBalls;
    //Copy this list of balls into a vector to allow subsequent OpenMP parallelisation
    std::map<const CXXBall*, std::map<const CXXBall*, std::vector<CXXCoord<CXXCoord_ftype> > > >::iterator reformattedEdgeEnd = reformattedEdges.end();
      for (std::map<const CXXBall*, std::map<const CXXBall *, std::vector<CXXCoord<CXXCoord_ftype> > > >::iterator reformattedEdge = reformattedEdges.begin();
         reformattedEdge != reformattedEdgeEnd;
         reformattedEdge++){
        trimmedBalls.push_back(reformattedEdge->first);
    };
    
    oldSize = aSurface->getChildSurfaces().size();
    aSurface->getChildSurfaces().resize(oldSize+trimmedBalls.size());
    CXXSurface *probeBallSurfacesArray = &(aSurface->getChildSurfaces()[oldSize]);
#if __APPLE__ && !defined _OPENMP
    dispatch_apply(trimmedBalls.size(), dispatch_get_global_queue(0, 0), ^(size_t i){
#elif defined _OPENMP
#warning Compiling for OMP
#pragma omp parallel for default(none) shared(trimmedBalls, contactMap, insideOrOutside, reformattedEdges, delta, aSurface, probeBallSurfacesArray, unitCellAtOriginForDelta) schedule(dynamic, 10)
#else
    for (int i=0; i<trimmedBalls.size(); i++){
#endif
        CXXSurface &ballSurface = probeBallSurfacesArray[i];
        const CXXBall &ball(*(trimmedBalls[i]));
        std::map<const CXXBall *, std::vector<const CXXBall *  > >::const_iterator contacts;
        contacts = contactMap.find(&ball);
        CXXNewHood ballHood;
        ballHood.initWith(&ball);
        if (contacts!=contactMap.end()){
            std::vector<const CXXBall *  >::const_iterator neighboursEnd(contacts->second.end());
            for (std::vector<const CXXBall *  >::const_iterator neighbour = contacts->second.begin();
                 neighbour!=neighboursEnd;
                 ++neighbour){
                ballHood.addBall(**neighbour);
            }
        }
        std::map<const CXXBall*, std::map<const CXXBall *, std::vector<CXXCoord<CXXCoord_ftype> > > >::const_iterator edgesOfBallIterator = reformattedEdges.find(&ball);
        if (edgesOfBallIterator != reformattedEdges.end()){
            std::map<const CXXBall *, std::vector<CXXCoord<CXXCoord_ftype> > > &edgesOfBallMap = const_cast< std::map<const CXXBall *, std::vector<CXXCoord<CXXCoord_ftype> > > & >(edgesOfBallIterator->second);
            ballHood.triangulateAsBallHoodInto(ballSurface, delta, edgesOfBallMap,true, insideOrOutside, unitCellAtOriginForDelta);
        }
#if defined TARGET_OS_MAC && !defined _OPENMP
    });
#else
    }
#endif
    
    
    return 0;
}

int CXXBall::ballContacts(std::vector<const CXXBall*  > &balls,
                          std::vector<const CXXBall*  > &contextBalls,
                          std::map<const CXXBall*, std::vector<const CXXBall*  > > &contactMap)
{
    int maxNBins = 20;
    if (balls.size() == 0) return 1;
    
    double maxBallRadius = -1.e30;
    std::vector<const CXXBall*  >::iterator ballsEnd = balls.end();
    for (std::vector<const CXXBall*  >::iterator ball = balls.begin();
         ball != ballsEnd;
         ++ball){
        maxBallRadius = (maxBallRadius > (*ball)->getRadius() ? maxBallRadius : (*ball)->getRadius() );
    }
    std::vector<const CXXBall*  >::iterator contextBallsEnd = contextBalls.end();
    for (std::vector<const CXXBall*  >::iterator ball = contextBalls.begin();
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
    for (std::vector<const CXXBall*  >::iterator ballIter = contextBalls.begin();
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
        int nMinimumBins = ceil((limits[i][1]-limits[i][0])/minimumBinSize);
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
    std::vector<std::vector<std::vector<std::vector<const CXXBall*  > > > >binnedballs;
    binnedballs.resize(nBins[0]);
    for (int i=0; i<nBins[0]; i++){
        binnedballs[i].resize(nBins[1]);
        for (int j=0; j<nBins[1]; j++){
            binnedballs[i][j].resize(nBins[2]);
        }
    }
    
    //Distribute balls among bins
    for (std::vector<const CXXBall*  >::iterator ballIter = contextBalls.begin();
         ballIter!=contextBallsEnd;
         ++ballIter){
        int iBin[3];
        for (int i=0; i<3; i++){
            iBin[i] = floor(((**ballIter)[i]-limits[i][0]) / binWidth[i]);
            iBin[i]  = max(0,iBin[i]);
            iBin[i] = min(nBins[i]-1,iBin[i]);
        }
        binnedballs[iBin[0]][iBin[1]][iBin[2]].push_back(*ballIter);
    }
    
    //To allow subsequent parallelization, create map entries in contactMap for each ball
    for (std::vector<const CXXBall*  >::iterator ballIter = balls.begin();
         ballIter!=ballsEnd;
         ++ballIter){
        const CXXBall *binBall = *ballIter;
        contactMap[binBall].resize(0);
    }
    std::cout << "Ready for parallel contact search \n";
    std::cout.flush();
    
    double *limitsPntr = &(limits[0][0]);
    double *binWidthPntr = binWidth;
    int *nBinsPntr = nBins;
    
    
#if __APPLE__ && !defined _OPENMP
    dispatch_apply(balls.size(), dispatch_get_global_queue(0, 0), ^(size_t iBall){
#elif defined _OPENMP
#pragma omp parallel for default(none) shared (binnedballs, ballRadiusX2, limits, binWidth, nBins, contactMap, limitsPntr, nBinsPntr, binWidthPntr, balls) schedule(dynamic,10)
    for (int iBall=0; iBall< balls.size(); iBall++){
#else
    for (int iBall=0; iBall< balls.size(); iBall++){
#endif
        int iBin[3];
        int *iBinPntr = iBin;
        for (int i=0; i<3; i++){
            iBinPntr[i] = floor(((*balls[iBall])[i]-limitsPntr[2*i]) / binWidthPntr[i]);
        }
        int startBinX = max(0,iBin[0]-1);
        int endBinX = min(nBinsPntr[0],iBinPntr[0]+2);
        int startBinY = max(0,iBinPntr[1]-1);
        int endBinY = min(nBinsPntr[1],iBinPntr[1]+2);
        int startBinZ = max(0,iBinPntr[2]-1);
        int endBinZ = min(nBinsPntr[2],iBinPntr[2]+2);
        
        const CXXBall &centralballRef(*balls[iBall]);
        std::vector<const CXXBall*  > &contactors = contactMap.find(&centralballRef)->second;
        
        //Compare this in turn with all balls in this or neighbouring bins
        for (int searchBinX = startBinX; searchBinX<endBinX; searchBinX++){
            for (int searchBinY = startBinY; searchBinY<endBinY; searchBinY++){
                for (int searchBinZ = startBinZ; searchBinZ<endBinZ; searchBinZ++){
                    const std::vector<const CXXBall*  > &searchBin(binnedballs[searchBinX][searchBinY][searchBinZ]);
                    std::vector<const CXXBall*  >::const_iterator binSearchEnd = searchBin.end();
                    for (std::vector<const CXXBall*  >::const_iterator otherball = searchBin.begin();
                         otherball!=binSearchEnd;
                         ++otherball){
                        const CXXBall &otherballRef(**otherball);
                        if (&centralballRef != &otherballRef){
                            if (!centralballRef.getCoord().isNearly(otherballRef.getCoord(), 0.00001)){
                                if (centralballRef.getCoord().isNearly(otherballRef.getCoord(), ballRadiusX2)){
                                    CXXCoord<CXXCoord_ftype>diff(centralballRef.getCoord() - otherballRef.getCoord());
                                    if (diff.get3DLengthSq()<pow(centralballRef.getRadius() + otherballRef.getRadius(), 2.0)){
                                        
                                        contactors.push_back(&otherballRef);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
#if defined TARGET_OS_MAC && !defined _OPENMP
    });
#else
    }
#endif
    
    return 0;
}

double CXXAtomBall::mostRecentDelta = 1.;
CXXSphereElement CXXAtomBall::unitSphereAtOrigin = CXXSphereElement(CXXCoord<CXXCoord_ftype>(0.,0.,0.), 1., 1.);
