/*
 * MoleculesToTriangles/CXXSurface/CXXSurfaceMaker.h
 *
 * Copyright 2016 by Martin Noble, University of Oxford
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
