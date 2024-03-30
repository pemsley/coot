//
//  CXXSurfaceMaker.hpp
//  MoleculesToTriangles
//
//  Created by Martin Noble on 9/5/16.
//  Copyright © 2016 MartinNobleSoftware. All rights reserved.
//

#ifndef CXXSurfaceMaker_h
#define CXXSurfaceMaker_h

#include "CXXSurface.h"

#include <stdio.h>
#include <vector>
#include <string>
#include <map>

#ifndef __MMDB_Manager__
#include "mmdb2/mmdb_manager.h"
#endif

class CXXVCN;
class CXXBall;
class CXXCircleNode;

class CXXSurfaceMaker
{
private:
    mmdb::Manager *allAtomsManager;
    std::vector<CXXSurface> childSurfaces;

public:
    CXXSurfaceMaker() { allAtomsManager = 0; };
    ~CXXSurfaceMaker(){};
    CXXSurfaceMaker(std::string path);
    CXXSurfaceMaker(mmdb::Manager *, const std::string selectionString);
    CXXSurfaceMaker(mmdb::Manager *, const std::string selectionString, const std::string contextString);
    CXXSurfaceMaker(mmdb::Manager *, const std::string selectionString, const std::string contextString, const double delta, const double probeRadius, const bool blend_edges);
    CXXSurfaceMaker(mmdb::Manager *, const int);
    CXXSurfaceMaker(mmdb::Manager *, const int, const double, const double, const bool);
    CXXSurfaceMaker(mmdb::Manager *, const int, const int, const double, const double, const bool);
    CXXSurfaceMaker(int);
    std::string report();
    int writeAsGrasp(std::string path);
    int readGraspFile(std::string fileName);

    std::vector<CXXSurface> &getChildSurfaces()
    {
        return childSurfaces;
    };
    mmdb::Manager *getMMDBManager() const;
    int calculateFromAtoms(mmdb::Manager *, const std::string, const std::string, const double, const double, const bool);
    int calculateFromAtoms(mmdb::Manager *, const int, const int, const double, const double, const double, const bool);
    int calculateFromAtoms(mmdb::Manager *, const int, const int, const double, const double, const bool);
    int calculateVDWFromAtoms(mmdb::Manager *, const int, const int, const double, const double, const double, const bool);
    int calculateAccessibleFromAtoms(mmdb::Manager *, const int, const int, const double, const double, const double, const bool);
    int calculateFromAtoms(mmdb::Manager *, const std::string, const double, const double, const bool);
    int calculateFromAtoms(mmdb::Manager *, const int, const double, const double, const bool);

    SurfaceParameters measuredProperties();

    int selectionStringToSelHnd(mmdb::Manager *allAtomsManager_in, std::string selectionString);
    double getAtomRadius(mmdb::Atom *);
    int assignAtom(mmdb::Manager *, int);
    int colorByAssignedAtom();
    int colorByColourArray(const std::vector<double *> &colours, mmdb::Manager *molHnd, int selHnd);

    void generateArrays(std::vector<CXXVCN> &vcns, std::vector<CXXCoord<float>> &accessibles, std::vector<mmdb::Atom *> &atoms);
    void memberHandleCentralAtoms(
        const int atomNr,
        const vector<const CXXBall *> *vdwBallPntrs,
        CXXSurface *elementSurfacesArray,
        const float radiusMultiplier,
        const float probeRadius,
        const float delta,
        const CXXSphereElement *unitSphereAtOrigin,
        const std::map<const CXXBall *, std::vector<const CXXBall *>> *contactMap,
        std::vector<CXXCircleNode> *splitReentrantProbesArray,
        const int selHnd
        );
};

#endif /* CXXSurfaceMaker_h */
