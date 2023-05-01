//
//  CXXSurfaceMaker.cpp
//  MoleculesToTriangles
//
//  Created by Martin Noble on 9/5/16.
//  Copyright © 2016 MartinNobleSoftware. All rights reserved.
//

#include "CXXSurfaceVertex.h"
#include "CXXSurfaceMaker.h"
#include "CXXVCN.h"
#include "mmdb2/mmdb_manager.h"
#include "mmdb2/mmdb_tables.h"
#include "CXXBall.h"
#include "CXXNewHood.h"
#include "CXXSurface.h"
#include <sstream>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

// can return -1 if name is invalid
//
long get_number_of_threads_by_system_call()
{

#ifdef _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
#else
    return sysconf(_SC_NPROCESSORS_CONF);
#endif
}

#if defined _OPENMP
#include <omp.h>
#else
#if __APPLE__
#include <dispatch/dispatch.h>
#else
#include <thread>
#include "ctpl_boost.h"
#endif
#endif

void CXXSurfaceMaker::generateArrays(std::vector<CXXVCN> &vcns, std::vector<CXXCoord<GLfloat>> &accessibles, std::vector<mmdb::Atom *> &atoms)
{
}

double CXXSurfaceMaker::getAtomRadius(mmdb::Atom *theAtom)
{
    // Here get handle of a radius data type from MMDB if such has been stored
    int iRadiusHandle = allAtomsManager->GetUDDHandle(mmdb::UDR_ATOM, "PerAtomRadius");
    double theRadius;
    if (iRadiusHandle > 0)
    {
        int success = theAtom->GetUDData(iRadiusHandle, theRadius);
        if (success != mmdb::UDDATA_Ok)
            theRadius = 1.8;
    }
    else
        theRadius = mmdb::getVdWaalsRadius(theAtom->element);
    return theRadius;
}

CXXSurfaceMaker::CXXSurfaceMaker(mmdb::Manager *allAtomsManager_in, const std::string selectionString,
                                 const std::string contextString)
{
    double delta = 30. * 2. * M_PI / 360.;
    double probeRadius = 1.4;
    bool blend_edges = false;
    calculateFromAtoms(allAtomsManager_in, selectionString, contextString, probeRadius, delta, blend_edges);
}

CXXSurfaceMaker::CXXSurfaceMaker(mmdb::Manager *allAtomsManager_in, const std::string selectionString,
                                 const std::string contextString, const double delta, const double probeRadius, const bool blend_edges)
{
    calculateFromAtoms(allAtomsManager_in, selectionString, contextString, probeRadius, delta, blend_edges);
}

CXXSurfaceMaker::CXXSurfaceMaker(mmdb::Manager *allAtomsManager_in, const std::string selectionString)
{
    double delta = 30. * 2. * M_PI / 360.;
    double probeRadius = 1.4;
    bool blend_edges = false;
    calculateFromAtoms(allAtomsManager_in, selectionString, probeRadius, delta, blend_edges);
}

CXXSurfaceMaker::CXXSurfaceMaker(mmdb::Manager *allAtomsManager_in, const int selHnd)
{
    double delta = 30. * 2. * M_PI / 360.;
    double probeRadius = 1.4;
    bool blend_edges = false;
    double radiusMultiplier = 1.;
    calculateFromAtoms(allAtomsManager_in, selHnd, selHnd, probeRadius, delta, radiusMultiplier, blend_edges);
}

CXXSurfaceMaker::CXXSurfaceMaker(mmdb::Manager *allAtomsManager_in, const int selHnd, const double delta, const double probeRadius, const bool blend_edges)
{
    double radiusMultiplier = 1.;
    calculateFromAtoms(allAtomsManager_in, selHnd, selHnd, probeRadius, delta, radiusMultiplier, blend_edges);
}

CXXSurfaceMaker::CXXSurfaceMaker(mmdb::Manager *allAtomsManager_in, const int selHnd, const int contextSelHnd,
                                 const double delta, const double probeRadius, const bool blend_edges)
{
    double radiusMultiplier = 1.;
    calculateFromAtoms(allAtomsManager_in, selHnd, contextSelHnd, probeRadius, delta, radiusMultiplier, blend_edges);
}

int CXXSurfaceMaker::calculateFromAtoms(mmdb::Manager *allAtomsManager_in, const std::string selectionString, const double probeRadius, const double delta, const bool blend_edges)
{
    int selHnd = selectionStringToSelHnd(allAtomsManager_in, selectionString);
    return calculateFromAtoms(allAtomsManager_in, selHnd, selHnd, probeRadius, delta, blend_edges);
}

int CXXSurfaceMaker::calculateFromAtoms(mmdb::Manager *allAtomsManager_in, const std::string selectionString, const std::string contextString, const double probeRadius, const double delta, const bool blend_edges)
{
    int selHnd = selectionStringToSelHnd(allAtomsManager_in, selectionString);
    int contextHnd = selectionStringToSelHnd(allAtomsManager_in, contextString);
    return calculateFromAtoms(allAtomsManager_in, selHnd, contextHnd, probeRadius, delta, blend_edges);
}

int CXXSurfaceMaker::calculateFromAtoms(mmdb::Manager *allAtomsManager_in, const int selHnd, const double probeRadius, const double delta, const bool blend_edges)
{
    return calculateFromAtoms(allAtomsManager_in, selHnd, selHnd, probeRadius, delta, blend_edges);
}

int CXXSurfaceMaker::calculateVDWFromAtoms(mmdb::Manager *allAtomsManager_in, const int selHnd, const int contextSelHnd, const double probeRadius, const double delta, const double radiusMultiplier, const bool blend_edges)
{
    allAtomsManager = allAtomsManager_in;

    int nSelAtoms;
    mmdb::Atom **SelAtom;
    allAtomsManager->GetSelIndex(selHnd, SelAtom, nSelAtoms);
    std::cout << "Surface selection includes " << nSelAtoms << "atoms" << endl;
    vector<const CXXBall *> vdwBallPntrs;
    for (int atomNr = 0; atomNr < nSelAtoms; atomNr++)
    {
        vdwBallPntrs.push_back(new CXXAtomBall(SelAtom[atomNr], radiusMultiplier * getAtomRadius(SelAtom[atomNr])));
    }
    map<mmdb::Atom *, const CXXBall *> mainAtoms;
    for (int atomNr = 0; atomNr < nSelAtoms; atomNr++)
    {
        mainAtoms[vdwBallPntrs[atomNr]->getAtomI()] = vdwBallPntrs[atomNr];
    }

    int nContextSelAtoms;
    mmdb::Atom **ContextSelAtom;
    allAtomsManager->GetSelIndex(contextSelHnd, ContextSelAtom, nContextSelAtoms);
    cout << "Context selection includes " << nSelAtoms << "atoms" << endl;

    vector<const CXXBall *> contextBallPntrs;
    int nUniqueContextAtoms = 0;
    int nSharedContextAtoms = 0;
    for (int atomNr = 0; atomNr < nContextSelAtoms; atomNr++)
    {
        map<mmdb::Atom *, const CXXBall *>::iterator equivalentMainAtom = mainAtoms.find(ContextSelAtom[atomNr]);
        if (equivalentMainAtom != mainAtoms.end())
        {
            contextBallPntrs.push_back(equivalentMainAtom->second);
            nSharedContextAtoms++;
        }
        else
        {
            contextBallPntrs.push_back(new CXXAtomBall(ContextSelAtom[atomNr], radiusMultiplier * getAtomRadius(ContextSelAtom[atomNr])));
            nUniqueContextAtoms++;
        }
    }
    std::cout << "nUniqueContextAtoms " << nUniqueContextAtoms << " nSharedContextAtoms " << nSharedContextAtoms << std::endl;

    CXXBall::triangulateBalls(vdwBallPntrs, contextBallPntrs, delta, this, CXXSphereElement::VDW);
    for (int i = 0; i < vdwBallPntrs.size(); i++)
    {
        if (vdwBallPntrs[i])
            delete static_cast<const CXXAtomBall *>(vdwBallPntrs[i]);
    }
    for (int i = 0; i < contextBallPntrs.size(); i++)
    {
        map<mmdb::Atom *, const CXXBall *>::iterator equivalentMainAtom = mainAtoms.find(ContextSelAtom[i]);
        if (equivalentMainAtom == mainAtoms.end())
        {
            if (contextBallPntrs[i])
                delete static_cast<const CXXAtomBall *>(contextBallPntrs[i]);
        }
    }
    report();
    return 0;
}

int CXXSurfaceMaker::calculateAccessibleFromAtoms(mmdb::Manager *allAtomsManager_in, const int selHnd, const int contextSelHnd, const double probeRadius, const double delta, const double radiusMultiplier, const bool blend_edges)
{
    allAtomsManager = allAtomsManager_in;

    int nSelAtoms;
    mmdb::Atom **SelAtom;
    allAtomsManager->GetSelIndex(selHnd, SelAtom, nSelAtoms);
    cout << "Surface selection includes " << nSelAtoms << "atoms" << endl;
    vector<const CXXBall *> vdwBallPntrs;
    for (int atomNr = 0; atomNr < nSelAtoms; atomNr++)
    {
        vdwBallPntrs.push_back(new CXXAtomBall(SelAtom[atomNr], probeRadius + (radiusMultiplier * getAtomRadius(SelAtom[atomNr]))));
    }
    map<mmdb::Atom *, const CXXBall *> mainAtoms;
    for (int atomNr = 0; atomNr < nSelAtoms; atomNr++)
    {
        mainAtoms[vdwBallPntrs[atomNr]->getAtomI()] = vdwBallPntrs[atomNr];
    }

    int nContextSelAtoms;
    mmdb::Atom **ContextSelAtom;
    allAtomsManager->GetSelIndex(contextSelHnd, ContextSelAtom, nContextSelAtoms);
    cout << "Context selection includes " << nSelAtoms << "atoms" << endl;

    vector<const CXXBall *> contextBallPntrs;
    int nUniqueContextAtoms = 0;
    int nSharedContextAtoms = 0;
    for (int atomNr = 0; atomNr < nContextSelAtoms; atomNr++)
    {
        map<mmdb::Atom *, const CXXBall *>::iterator equivalentMainAtom = mainAtoms.find(ContextSelAtom[atomNr]);
        if (equivalentMainAtom != mainAtoms.end())
        {
            contextBallPntrs.push_back(equivalentMainAtom->second);
            nSharedContextAtoms++;
        }
        else
        {
            contextBallPntrs.push_back(new CXXAtomBall(ContextSelAtom[atomNr], probeRadius + (radiusMultiplier * getAtomRadius(ContextSelAtom[atomNr]))));
            nUniqueContextAtoms++;
        }
    }
    std::cout << "nUniqueContextAtoms " << nUniqueContextAtoms << " nSharedContextAtoms " << nSharedContextAtoms << std::endl;

    CXXBall::triangulateBalls(vdwBallPntrs, contextBallPntrs, delta, this, CXXSphereElement::Accessible);
    for (int i = 0; i < vdwBallPntrs.size(); i++)
    {
        if (vdwBallPntrs[i])
            delete static_cast<const CXXAtomBall *>(vdwBallPntrs[i]);
    }
    for (int i = 0; i < contextBallPntrs.size(); i++)
    {
        map<mmdb::Atom *, const CXXBall *>::iterator equivalentMainAtom = mainAtoms.find(ContextSelAtom[i]);
        if (equivalentMainAtom == mainAtoms.end())
        {
            if (contextBallPntrs[i])
                delete static_cast<const CXXAtomBall *>(contextBallPntrs[i]);
        }
    }
    report();
    return 0;
}

int CXXSurfaceMaker::calculateFromAtoms(mmdb::Manager *allAtomsManager_in, const int selHnd, const int contextSelHnd, const double probeRadius, const double delta, const bool blend_edges)
{
    double radiusMultiplier = 1.0;
    return calculateFromAtoms(allAtomsManager_in, selHnd, contextSelHnd, probeRadius, delta, radiusMultiplier, blend_edges);
}

void handleCentralAtom(
    const int atomNr,
    CXXSurfaceMaker *parentSurfaceMaker,
    const vector<const CXXBall *> *vdwBallPntrs,
    CXXSurface *elementSurfacesArray,
    const float radiusMultiplier,
    const float probeRadius,
    const float delta,
    const CXXSphereElement *unitSphereAtOrigin,
    const std::map<const CXXBall *, std::vector<const CXXBall *>> *contactMap,
    std::vector<CXXCircleNode> *splitReentrantProbesArray,
    const int selHnd)
{
    /*
          if (!(atomNr % 100) || atomNr == nSelAtoms - 1)
          {
    #if __APPLE__ && !defined _OPENMP
              dispatch_sync(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_LOW, NULL), ^(void) {
                cout << "Dealing with atom number " << atomNr << endl;
              });
    #elif defined _OPENMP
    #pragma omp critical(cout)
                cout << "Dealing with atom number " << atomNr << " on thread " << omp_get_thread_num() << endl;
    #endif
  }
    */

    mmdb::Atom *centralAtom = static_cast<const CXXAtomBall *>((*vdwBallPntrs)[atomNr])->getAtomI();

    CXXSurface &elementSurface = elementSurfacesArray[atomNr];

    double radiusOfAtom1 = radiusMultiplier * parentSurfaceMaker->getAtomRadius(centralAtom);
    CXXNewHood theNewHood;
    theNewHood.initWith(centralAtom, radiusOfAtom1, probeRadius);

    // We have precalculated neighbours of the central atom, and now can use that
    // to our advantage
    std::map<const CXXBall *, std::vector<const CXXBall *>>::const_iterator contactMapIter = contactMap->find((*vdwBallPntrs)[atomNr]);

    const std::vector<const CXXBall *> &neighbours(contactMapIter->second);
    for (unsigned int sphereAtomNr = 0; sphereAtomNr < neighbours.size(); sphereAtomNr++)
    {
        theNewHood.addBall(*neighbours[sphereAtomNr]);
    }

    // Find the non-hidden segments of the circles
    theNewHood.findSegments();
    if (CXXNewHood::containsDrawable(theNewHood))
    {
        theNewHood.triangulateAsRegularHoodInto(elementSurface, delta, unitSphereAtOrigin);
        theNewHood.identifyUniqueNodes(splitReentrantProbesArray[atomNr], selHnd);
        elementSurface.compress(0.00001);
    }
}

void ctplHandleCentralAtom(int id,
                           const int atomNr,
                           CXXSurfaceMaker *parentSurfaceMaker,
                           const vector<const CXXBall *> *vdwBallPntrs,
                           CXXSurface *elementSurfacesArray,
                           const float radiusMultiplier,
                           const float probeRadius,
                           const float delta,
                           const CXXSphereElement *unitSphereAtOrigin,
                           const std::map<const CXXBall *, std::vector<const CXXBall *>> *contactMap,
                           std::vector<CXXCircleNode> *splitReentrantProbesArray,
                           const int selHnd)
{
    handleCentralAtom(atomNr, parentSurfaceMaker, vdwBallPntrs, elementSurfacesArray, radiusMultiplier, probeRadius, delta, unitSphereAtOrigin, contactMap, splitReentrantProbesArray, selHnd);
}

int CXXSurfaceMaker::calculateFromAtoms(mmdb::Manager *allAtomsManager_in, const int selHnd, const int contextSelHnd, const double probeRadius, const double delta, const double radiusMultiplier, const bool blend_edges)
{
    allAtomsManager = allAtomsManager_in;

    int nSelAtoms;
    mmdb::Atom **SelAtom;
    allAtomsManager->GetSelIndex(selHnd, SelAtom, nSelAtoms);
    cout << "Surface selection includes " << nSelAtoms << "atoms" << endl;
    vector<const CXXBall *> vdwBallPntrs;
    for (int atomNr = 0; atomNr < nSelAtoms; atomNr++)
    {
        vdwBallPntrs.push_back(new CXXAtomBall(SelAtom[atomNr], probeRadius + (radiusMultiplier * getAtomRadius(SelAtom[atomNr]))));
    }
    int nContextSelAtoms;
    mmdb::Atom **ContextSelAtom;
    allAtomsManager->GetSelIndex(contextSelHnd, ContextSelAtom, nContextSelAtoms);
    cout << "Context selection includes " << nSelAtoms << "atoms" << endl;
    vector<const CXXBall *> contextBallPntrs;
    for (int atomNr = 0; atomNr < nContextSelAtoms; atomNr++)
    {
        contextBallPntrs.push_back(new CXXAtomBall(ContextSelAtom[atomNr], probeRadius + (radiusMultiplier * getAtomRadius(ContextSelAtom[atomNr]))));
    }

    // Precalculate contacts
    std::map<const CXXBall *, std::vector<const CXXBall *>> contactMap;
    CXXBall::ballContacts(vdwBallPntrs, contextBallPntrs, contactMap);
    std::cout << "Established contact map\n";

    // Start up with reentrant Prbes separated per-atom...removes one locking state
    // for multi-threadig
    std::vector<std::vector<CXXCircleNode>> splitReentrantProbes(nSelAtoms);
    std::vector<CXXCircleNode> *splitReentrantProbesArray = &(splitReentrantProbes[0]);

    CXXSphereElement unitSphereAtOrigin(CXXCoord<CXXCoord_ftype>(0., 0., 0.), 1., delta);

    size_t oldSize = childSurfaces.size();
    childSurfaces.resize(oldSize + nSelAtoms);
    CXXSurface *elementSurfacesArray = &(childSurfaces[oldSize]);

#if !defined __APPLE__ && !defined _OPENMP
    // vector container stores threads
    unsigned int n_threads = 2;//get_number_of_threads_by_system_call();
    ctpl::thread_pool the_thread_pool(n_threads);
    std::cout << "n_threads" << n_threads << "\n";
#endif

#if __APPLE__ && !defined _OPENMP
    dispatch_apply(nSelAtoms, dispatch_get_global_queue(0, 0), ^(size_t atomNr) {
      handleCentralAtom(atomNr, this, &vdwBallPntrs, elementSurfacesArray, radiusMultiplier,
                        probeRadius, delta, &unitSphereAtOrigin, &contactMap, splitReentrantProbesArray,
                        selHnd);
#elif defined _OPENMP
#warning Compiling for OMP
      // compilation failure 9.2.1
      // #pragma omp parallel for default(none) shared(vdwBallPntrs,contactMap, splitReentrantProbesArray, nSelAtoms, cout, unitSphereAtOrigin, elementSurfacesArray) schedule(dynamic, 100)
      // 20230120-PE ... and again
      // #pragma omp parallel for default(none) shared(vdwBallPntrs, &contactMap, splitReentrantProbesArray, nSelAtoms, cout, unitSphereAtOrigin, elementSurfacesArray, radiusMultiplier, probeRadius, delta, selHnd) schedule(dynamic, 100)
#pragma omp parallel for default(none) shared(vdwBallPntrs, contactMap, splitReentrantProbesArray, nSelAtoms, cout, unitSphereAtOrigin, elementSurfacesArray, radiusMultiplier, probeRadius, delta, selHnd) schedule(dynamic, 100)
    for (int atomNr = 0; atomNr < nSelAtoms; atomNr++)
    {
        handleCentralAtom(atomNr, this, &vdwBallPntrs, elementSurfacesArray, radiusMultiplier,
                          probeRadius, delta, &unitSphereAtOrigin, &contactMap, splitReentrantProbesArray,
                          selHnd);
#else
    for (int atomNr = 0; atomNr < nSelAtoms; atomNr++)
    {
        if (true)
        {
            the_thread_pool.push(ctplHandleCentralAtom, atomNr, this, &vdwBallPntrs, elementSurfacesArray, radiusMultiplier,
                                 probeRadius, delta, &unitSphereAtOrigin, &contactMap, splitReentrantProbesArray,
                                 selHnd);
        }
        else
        {
            handleCentralAtom(atomNr, this, &vdwBallPntrs, elementSurfacesArray, radiusMultiplier,
                              probeRadius, delta, &unitSphereAtOrigin, &contactMap, splitReentrantProbesArray,
                              selHnd);
        }
#endif
        /*
        CXXSurface &elementSurface = elementSurfacesArray[atomNr];

        double radiusOfAtom1 = radiusMultiplier * getAtomRadius(centralAtom);
        CXXNewHood theNewHood;
        theNewHood.initWith(centralAtom, radiusOfAtom1, probeRadius);

        //We have precalculated neighbours of the central atom, and now can use that
        //to our advantage
        std::map<const CXXBall *, std::vector<const CXXBall *  > >::const_iterator contactMapIter = contactMap.find(vdwBallPntrs[atomNr]);

        const std::vector<const CXXBall *  > &neighbours(contactMapIter->second);
        for (unsigned int sphereAtomNr = 0; sphereAtomNr < neighbours.size(); sphereAtomNr++) {
            theNewHood.addBall(*neighbours[sphereAtomNr]);
        }

        //Find the non-hidden segments of the circles
        theNewHood.findSegments();
        if (CXXNewHood::containsDrawable(theNewHood)){
            theNewHood.triangulateAsRegularHoodInto(elementSurface, delta, &unitSphereAtOrigin);
            theNewHood.identifyUniqueNodes(splitReentrantProbesArray[atomNr], selHnd);
            elementSurface.compress(0.00001);
        }
    */
#if defined TARGET_OS_MAC && !defined _OPENMP
    });
#else
    }
#if __APPLE__ && !defined _OPENMP
    the_thread_pool.stop(true);
#endif
#endif

    vector<const CXXBall *> reentrantProbes;
    for (int i = 0; i < nSelAtoms; i++)
    {
        std::vector<CXXCircleNode>::iterator reentrantProbesEnd(splitReentrantProbes[i].end());
        for (std::vector<CXXCircleNode>::iterator reentrantProbe = splitReentrantProbes[i].begin();
             reentrantProbe != reentrantProbesEnd;
             ++reentrantProbe)
        {
            reentrantProbes.push_back(new CXXReentrantProbeBall(*reentrantProbe, selHnd, probeRadius));
        }
    }
    splitReentrantProbes.clear();
    for (int i = 0; i < vdwBallPntrs.size(); i++)
    {
        if (vdwBallPntrs[i])
            delete static_cast<const CXXAtomBall *>(vdwBallPntrs[i]);
    }
    for (int i = 0; i < contextBallPntrs.size(); i++)
    {
        if (contextBallPntrs[i])
            delete static_cast<const CXXAtomBall *>(contextBallPntrs[i]);
    }
    CXXBall::triangulateBalls(reentrantProbes, reentrantProbes, delta, this, CXXSphereElement::Reentrant);
    for (int i = 0; i < reentrantProbes.size(); i++)
    {
        delete static_cast<const CXXReentrantProbeBall *>(reentrantProbes[i]);
    }

    if (blend_edges)
    {
        cout << "Starting to blend edges" << endl;
        assignAtom(allAtomsManager, selHnd);
    }

    report();
    // cout << "Starting default colour surface" <<endl;
    // colorByAssignedAtom();
    // cout << "Finished default colour surface" <<endl;

    return 0;
}

int CXXSurfaceMaker::selectionStringToSelHnd(mmdb::Manager *allAtomsManager_in, std::string selectionString)
{
    int selHnd = allAtomsManager_in->NewSelection();
    char *pstring = (char *)malloc(sizeof(selectionString.c_str()) + 1);
    strcpy(pstring, selectionString.c_str());
    allAtomsManager_in->Select(selHnd, mmdb::STYPE_ATOM, pstring, mmdb::SKEY_NEW);
    free(pstring);
    return selHnd;
}

int CXXSurfaceMaker::assignAtom(mmdb::Manager *manager, int intValue)
{
    for (std::vector<CXXSurface>::iterator surfaceIter = childSurfaces.begin();
         surfaceIter != childSurfaces.end();
         surfaceIter++)
    {
        surfaceIter->assignAtom(manager, intValue);
    }
    return 0;
}

SurfaceParameters CXXSurfaceMaker::measuredProperties()
{
    std::vector<SurfaceParameters> surfaceParametersVector(childSurfaces.size());
    SurfaceParameters *surfaceParametersArray = &(surfaceParametersVector[0]);
    CXXSurface *childSurfacesArray = &(childSurfaces[0]);

#if __APPLE__ && !defined _OPENMP
    dispatch_apply(childSurfaces.size(), dispatch_get_global_queue(0, 0), ^(size_t iChildSurface) {
#elif defined _OPENMP
#warning Compiling for OMP
#pragma omp parallel for default(none) shared(surfaceParametersArray, childSurfacesArray) schedule(dynamic, 100)
    for (int iChildSurface = 0; iChildSurface < childSurfaces.size(); iChildSurface++)
    {
#else
        for (int iChildSurface = 0; iChildSurface < childSurfaces.size(); iChildSurface++)
        {
#endif
      surfaceParametersArray[iChildSurface] = childSurfacesArray[iChildSurface].measuredProperties();
#if defined TARGET_OS_MAC && !defined _OPENMP
    });
#else
    }
#endif

    SurfaceParameters surfaceParameters;
    for (vector<SurfaceParameters>::iterator childSurfaceParametersIter = surfaceParametersVector.begin();
         childSurfaceParametersIter != surfaceParametersVector.end();
         childSurfaceParametersIter++)
    {
        surfaceParameters += *childSurfaceParametersIter;
    }
    return surfaceParameters;
}

std::string CXXSurfaceMaker::report()
{
    std::ostringstream output;

    SurfaceParameters surfaceParameters = measuredProperties();

    output << "This surface has " << surfaceParameters.nVertices << " vertices and " << surfaceParameters.nTriangles
           << " triangles\n";
    output << "The Molecular surface area is " << surfaceParameters.MSA << endl;
    output << "The Accessible surface area is " << surfaceParameters.ASA << endl;
    for (std::map<std::string, double>::iterator prop = surfaceParameters.pMins.begin();
         prop != surfaceParameters.pMins.end();
         prop++)
    {
        output << "Property " << prop->first << " Range " << surfaceParameters.pMins[prop->first] << " to " << surfaceParameters.pMaxes[prop->first] << " mean " << surfaceParameters.pMeans[prop->first] << endl;
    }

    return output.str();
}

int CXXSurfaceMaker::writeAsGrasp(std::string path)
{
    CXXSurface collectedSurface;
    for (vector<CXXSurface>::iterator subSurfIter = childSurfaces.begin();
         subSurfIter != childSurfaces.end();
         subSurfIter++)
    {
        collectedSurface.appendSurface(*subSurfIter);
    }
    return collectedSurface.writeAsGrasp(path);
}

int CXXSurfaceMaker::readGraspFile(std::string fileName)
{
    childSurfaces.push_back(CXXSurface());
    return childSurfaces.back().readGraspFile(fileName);
}
