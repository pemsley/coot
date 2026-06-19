/*
 * MoleculesToTriangles/CXXSurface/CXXSurfaceMaker-gemmi.cc
 *
 * gemmi-native twin of CXXSurfaceMaker.cpp. Driven by a gemmi::Structure/Model and
 * selection sets (std::set<const gemmi::Atom*>). Atom radii come from CXXUtils-gemmi
 * get_atom_radius(atom, residue) at ball-construction time (residue context in scope);
 * downstream the per-atom radius is recovered from the ball (ball->getRadius()-probe).
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#define _USE_MATH_DEFINES
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif
#include "CXXSurfaceVertex.h"
#include "CXXSurfaceMaker-gemmi.hh"
#include "CXXVCN.h"
#include "CXXBall-gemmi.hh"
#include "CXXNewHood-gemmi.hh"
#include "CXXSurface-gemmi.hh"
#include "CXXSphereElement-gemmi.hh"
#include "CXXCircleNode-gemmi.hh"
#include "CXXUtils-gemmi.hh"
#include <sstream>

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

using namespace std;
using namespace coot::m2t;

void CXXSurfaceMaker::generateArrays(std::vector<CXXVCN> &vcns, std::vector<CXXCoord<float>> &accessibles, std::vector<const gemmi::Atom *> &atoms)
{
}

double CXXSurfaceMaker::getAtomRadius(const gemmi::Atom *atom, const gemmi::Residue &residue)
{
   return coot::m2t::get_atom_radius(*atom, residue);
}

CXXSurfaceMaker::CXXSurfaceMaker(gemmi::Structure *structure, gemmi::Model &model,
                                 const std::set<const gemmi::Atom*> &selSet,
                                 const std::set<const gemmi::Atom*> &contextSet,
                                 double delta, double probeRadius, bool blend_edges)
{
   allAtomsManager = 0;
   calculateFromAtoms(structure, model, selSet, contextSet, delta, probeRadius, blend_edges);
}

void CXXSurfaceMaker::memberHandleCentralAtoms(
   const int atomNr,
   const std::vector<const CXXBall *> *vdwBallPntrs,
   CXXSurface *elementSurfacesArray,
   const float radiusMultiplier,
   const float probeRadius,
   const float delta,
   const CXXSphereElement *unitSphereAtOrigin,
   const std::map<const CXXBall *, std::vector<const CXXBall *>> *contactMap,
   std::vector<CXXCircleNode> *splitReentrantProbesArray,
   const std::set<const gemmi::Atom*> &selSet)
{
   const gemmi::Atom *centralAtom = static_cast<const CXXAtomBall *>((*vdwBallPntrs)[atomNr])->getAtomI();

   CXXSurface &elementSurface = elementSurfacesArray[atomNr];

   // radiusMultiplier*getAtomRadius(centralAtom) == ballRadius - probeRadius (no residue context needed here)
   double radiusOfAtom1 = (*vdwBallPntrs)[atomNr]->getRadius() - probeRadius;
   (void) radiusMultiplier;
   CXXNewHood theNewHood;
   theNewHood.initWith(centralAtom, radiusOfAtom1, probeRadius);

   std::map<const CXXBall *, std::vector<const CXXBall *>>::const_iterator contactMapIter = contactMap->find((*vdwBallPntrs)[atomNr]);

   const std::vector<const CXXBall *> &neighbours(contactMapIter->second);
   for (unsigned int sphereAtomNr = 0; sphereAtomNr < neighbours.size(); sphereAtomNr++)
   {
      theNewHood.addBall(*neighbours[sphereAtomNr]);
   }

   theNewHood.findSegments();
   if (CXXNewHood::containsDrawable(theNewHood))
   {
      theNewHood.triangulateAsRegularHoodInto(elementSurface, delta, unitSphereAtOrigin);
      theNewHood.identifyUniqueNodes(splitReentrantProbesArray[atomNr], selSet);
      elementSurface.compress(0.00001);
   }
}

namespace {
   void handleCentralAtom(
      const int atomNr,
      CXXSurfaceMaker *parentSurfaceMaker,
      const std::vector<const CXXBall *> *vdwBallPntrs,
      CXXSurface *elementSurfacesArray,
      const float radiusMultiplier,
      const float probeRadius,
      const float delta,
      const CXXSphereElement *unitSphereAtOrigin,
      const std::map<const CXXBall *, std::vector<const CXXBall *>> *contactMap,
      std::vector<CXXCircleNode> *splitReentrantProbesArray,
      const std::set<const gemmi::Atom*> &selSet)
   {
      parentSurfaceMaker->memberHandleCentralAtoms(atomNr, vdwBallPntrs, elementSurfacesArray, radiusMultiplier, probeRadius,
                                                   delta, unitSphereAtOrigin, contactMap, splitReentrantProbesArray, selSet);
   }
}

int CXXSurfaceMaker::calculateFromAtoms(gemmi::Structure *structure, gemmi::Model &model,
                                        const std::set<const gemmi::Atom*> &selSet,
                                        const std::set<const gemmi::Atom*> &contextSet,
                                        double delta, double probeRadius, bool blend_edges)
{
   allAtomsManager = structure;
   double radiusMultiplier = 1.0;

   // Build balls by iterating the model (residue context needed for the radius lookup)
   vector<const CXXBall *> vdwBallPntrs;
   vector<const CXXBall *> contextBallPntrs;
   for (gemmi::Chain &chain : model.chains)
   {
      for (gemmi::Residue &res : chain.residues)
      {
         for (gemmi::Atom &atom : res.atoms)
         {
            const gemmi::Atom *ap = &atom;
            double r = probeRadius + (radiusMultiplier * getAtomRadius(ap, res));
            if (selSet.count(ap)) vdwBallPntrs.push_back(new CXXAtomBall(ap, r));
            if (contextSet.count(ap)) contextBallPntrs.push_back(new CXXAtomBall(ap, r));
         }
      }
   }
   int nSelAtoms = int(vdwBallPntrs.size());
   cout << "Surface selection includes " << nSelAtoms << " atoms" << endl;
   cout << "Context selection includes " << contextBallPntrs.size() << " atoms" << endl;

   // Precalculate contacts
   std::map<const CXXBall *, std::vector<const CXXBall *>> contactMap;
   CXXBall::ballContacts(vdwBallPntrs, contextBallPntrs, contactMap);
   std::cout << "Established contact map\n";

   std::vector<std::vector<CXXCircleNode>> splitReentrantProbes(nSelAtoms);
   std::vector<CXXCircleNode> *splitReentrantProbesArray = &(splitReentrantProbes[0]);

   CXXSphereElement unitSphereAtOrigin(CXXCoord<CXXCoord_ftype>(0., 0., 0.), 1., delta);

   size_t oldSize = childSurfaces.size();
   childSurfaces.resize(oldSize + nSelAtoms);
   CXXSurface *elementSurfacesArray = &(childSurfaces[oldSize]);

#if !defined __APPLE__ && !defined _OPENMP
   unsigned int n_threads = 2;
   ctpl::thread_pool the_thread_pool(n_threads);
#endif

#if __APPLE__ && !defined _OPENMP
   dispatch_apply(nSelAtoms, dispatch_get_global_queue(0, 0), ^(size_t atomNr) {
      handleCentralAtom(atomNr, this, &vdwBallPntrs, elementSurfacesArray, radiusMultiplier,
                        probeRadius, delta, &unitSphereAtOrigin, &contactMap, splitReentrantProbesArray,
                        selSet);
#elif defined _OPENMP
#pragma omp parallel for default(none) shared(vdwBallPntrs, contactMap, splitReentrantProbesArray, nSelAtoms, cout, unitSphereAtOrigin, elementSurfacesArray, radiusMultiplier, probeRadius, delta, selSet) schedule(dynamic, 100)
   for (int atomNr = 0; atomNr < nSelAtoms; atomNr++)
   {
      handleCentralAtom(atomNr, this, &vdwBallPntrs, elementSurfacesArray, radiusMultiplier,
                        probeRadius, delta, &unitSphereAtOrigin, &contactMap, splitReentrantProbesArray,
                        selSet);
#else
   for (int atomNr = 0; atomNr < nSelAtoms; atomNr++)
   {
      handleCentralAtom(atomNr, this, &vdwBallPntrs, elementSurfacesArray, radiusMultiplier,
                        probeRadius, delta, &unitSphereAtOrigin, &contactMap, splitReentrantProbesArray,
                        selSet);
#endif
#if defined TARGET_OS_MAC && !defined _OPENMP
   });
#else
   }
#if !defined __APPLE__ && !defined _OPENMP
   the_thread_pool.stop(true);
#endif
#endif

   vector<const CXXBall *> reentrantProbes;
   for (int i = 0; i < nSelAtoms; i++)
   {
      std::vector<CXXCircleNode>::iterator reentrantProbesEnd(splitReentrantProbes[i].end());
      for (std::vector<CXXCircleNode>::iterator reentrantProbe = splitReentrantProbes[i].begin();
           reentrantProbe != reentrantProbesEnd; ++reentrantProbe)
      {
         reentrantProbes.push_back(new CXXReentrantProbeBall(*reentrantProbe, selSet, probeRadius));
      }
   }
   splitReentrantProbes.clear();
   for (int i = 0; i < (int)vdwBallPntrs.size(); i++)
   {
      if (vdwBallPntrs[i]) delete static_cast<const CXXAtomBall *>(vdwBallPntrs[i]);
   }
   for (int i = 0; i < (int)contextBallPntrs.size(); i++)
   {
      if (contextBallPntrs[i]) delete static_cast<const CXXAtomBall *>(contextBallPntrs[i]);
   }
   CXXBall::triangulateBalls(reentrantProbes, reentrantProbes, delta, this, CXXSphereElement::Reentrant);
   for (int i = 0; i < (int)reentrantProbes.size(); i++)
   {
      delete static_cast<const CXXReentrantProbeBall *>(reentrantProbes[i]);
   }

   if (blend_edges)
   {
      cout << "Starting to blend edges" << endl;
      assignAtom(structure, selSet);
   }

   return 0;
}

namespace {
   // Ball-only surface (VdW / Accessible): radius = radiusOffset + radiusMultiplier*atomRadius.
   int ballSurface(CXXSurfaceMaker *maker, gemmi::Model &model,
                   const std::set<const gemmi::Atom*> &selSet,
                   const std::set<const gemmi::Atom*> &contextSet,
                   double delta, double radiusOffset, double radiusMultiplier, int sense)
   {
      std::vector<const CXXBall *> vdwBallPntrs;
      std::map<const gemmi::Atom *, const CXXBall *> mainAtoms;
      for (gemmi::Chain &chain : model.chains)
         for (gemmi::Residue &res : chain.residues)
            for (gemmi::Atom &atom : res.atoms) {
               const gemmi::Atom *ap = &atom;
               if (selSet.count(ap)) {
                  CXXAtomBall *b = new CXXAtomBall(ap, radiusOffset + (radiusMultiplier * maker->getAtomRadius(ap, res)));
                  vdwBallPntrs.push_back(b);
                  mainAtoms[ap] = b;
               }
            }
      std::vector<const CXXBall *> contextBallPntrs;
      std::vector<bool> contextUnique;
      for (gemmi::Chain &chain : model.chains)
         for (gemmi::Residue &res : chain.residues)
            for (gemmi::Atom &atom : res.atoms) {
               const gemmi::Atom *ap = &atom;
               if (contextSet.count(ap)) {
                  std::map<const gemmi::Atom *, const CXXBall *>::iterator it = mainAtoms.find(ap);
                  if (it != mainAtoms.end()) { contextBallPntrs.push_back(it->second); contextUnique.push_back(false); }
                  else { contextBallPntrs.push_back(new CXXAtomBall(ap, radiusOffset + (radiusMultiplier * maker->getAtomRadius(ap, res)))); contextUnique.push_back(true); }
               }
            }
      CXXBall::triangulateBalls(vdwBallPntrs, contextBallPntrs, delta, maker, sense);
      for (size_t i = 0; i < vdwBallPntrs.size(); i++) delete static_cast<const CXXAtomBall *>(vdwBallPntrs[i]);
      for (size_t i = 0; i < contextBallPntrs.size(); i++)
         if (contextUnique[i]) delete static_cast<const CXXAtomBall *>(contextBallPntrs[i]);
      return 0;
   }
}

int CXXSurfaceMaker::calculateVDWFromAtoms(gemmi::Structure *structure, gemmi::Model &model,
                                           const std::set<const gemmi::Atom*> &selSet,
                                           const std::set<const gemmi::Atom*> &contextSet,
                                           double delta, double probeRadius, double radiusMultiplier)
{
   allAtomsManager = structure;
   return ballSurface(this, model, selSet, contextSet, delta, 0.0, radiusMultiplier, CXXSphereElement::VDW);
}

int CXXSurfaceMaker::calculateAccessibleFromAtoms(gemmi::Structure *structure, gemmi::Model &model,
                                                  const std::set<const gemmi::Atom*> &selSet,
                                                  const std::set<const gemmi::Atom*> &contextSet,
                                                  double delta, double probeRadius, double radiusMultiplier)
{
   allAtomsManager = structure;
   return ballSurface(this, model, selSet, contextSet, delta, probeRadius, radiusMultiplier, CXXSphereElement::Accessible);
}

int CXXSurfaceMaker::assignAtom(gemmi::Structure *structure, const std::set<const gemmi::Atom*> &selSet)
{
   for (std::vector<CXXSurface>::iterator surfaceIter = childSurfaces.begin();
        surfaceIter != childSurfaces.end(); surfaceIter++)
   {
      surfaceIter->assignAtom(structure, selSet);
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
#pragma omp parallel for default(none) shared(surfaceParametersArray, childSurfacesArray) schedule(dynamic, 100)
   for (int iChildSurface = 0; iChildSurface < (int)childSurfaces.size(); iChildSurface++)
   {
#else
   for (int iChildSurface = 0; iChildSurface < (int)childSurfaces.size(); iChildSurface++)
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
        childSurfaceParametersIter != surfaceParametersVector.end(); childSurfaceParametersIter++)
   {
      surfaceParameters += *childSurfaceParametersIter;
   }
   return surfaceParameters;
}
