/*
 * MoleculesToTriangles/CXXSurface/CXXSurface-gemmi.hh
 *
 * gemmi-native twin of CXXSurface.h. allAtomsManager (mmdb::Manager*) becomes a
 * gemmi::Structure*; the per-vertex atom is stored via the opaque void* pointer
 * property (it carries a const gemmi::Atom*). Electrostatics colour-array helpers
 * are deferred to Phase 6.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */
#ifndef CXXSurface_gemmi_included
#define CXXSurface_gemmi_included
#include <string>
#include <vector>
#include <map>
#include <set>
#include "CXXTriangle-gemmi.hh"
#include "CXXSurfaceVertex.h"
#include "CXXCoord.h"
#include <gemmi/model.hpp>

namespace coot {
   namespace m2t {

      class CXXSphereElement;
      class CXXTorusElement;

      class SurfaceParameters {
      public:
         size_t nVertices;
         size_t nTriangles;
         double MSA;
         double ASA;
         std::map<std::string, double> pMins;
         std::map<std::string, double> pMaxes;
         std::map<std::string, double> pMeans;
         SurfaceParameters() : nVertices(0), nTriangles(0), MSA(0.), ASA(0.) {}
         SurfaceParameters& operator += (const SurfaceParameters &otherParameters) {
            MSA += otherParameters.MSA;
            ASA += otherParameters.ASA;
            for (std::map<std::string,double>::const_iterator prop=otherParameters.pMins.begin(); prop!=otherParameters.pMins.end(); prop++) {
               if (pMins.find(prop->first)!=pMins.end()) pMins[prop->first] = std::min(pMins[prop->first], prop->second);
               else pMins[prop->first] = prop->second;
            }
            for (std::map<std::string,double>::const_iterator prop=otherParameters.pMaxes.begin(); prop!=otherParameters.pMaxes.end(); prop++) {
               if (pMaxes.find(prop->first)!=pMaxes.end()) pMaxes[prop->first] = std::max(pMaxes[prop->first], prop->second);
               else pMaxes[prop->first] = prop->second;
            }
            for (std::map<std::string,double>::const_iterator prop=otherParameters.pMeans.begin(); prop!=otherParameters.pMeans.end(); prop++) {
               if (pMeans.find(prop->first)!=pMeans.end()) pMeans[prop->first] = (pMeans[prop->first] + prop->second) / (nVertices + otherParameters.nVertices);
               else pMeans[prop->first] = prop->second;
            }
            nTriangles += otherParameters.nTriangles;
            nVertices += otherParameters.nVertices;
            return *this;
         }
      };

      typedef std::map<std::string, size_t> StringIntMap;

      class CXXSurface {
      private:
         std::string name;
         StringIntMap vectors;
         StringIntMap scalars;
         StringIntMap pointers;
         std::vector<CXXTriangle> triangles;
         std::vector<CXXSurfaceVertex> vertices;
         gemmi::Structure* allAtomsManager;
         int init();
         size_t nTriangles;
         char fileName[512];
      public:
         std::string report();

         const std::vector<CXXTriangle> &getTriangles() const { return triangles; }
         const std::vector<CXXSurfaceVertex> &getVertices() const { return vertices; }
         const StringIntMap &vectorNames() const { return vectors; }
         const StringIntMap &scalarNames() const { return scalars; }
         const StringIntMap &pointerNames() const { return pointers; }

         const CXXCoord<CXXCoord_ftype>& coordRef(size_t coordType, size_t iTriangle, size_t corner) const;
         const CXXCoord<CXXCoord_ftype>& coordRef(size_t coordType, size_t iVertex) const;
         int getCoord(const std::string &type, const size_t iTriangle, const size_t corner, double *buffer);
         int getCoord(const std::string &type, const size_t iVertex, double *buffer);
         int getCoord(const size_t handle, const size_t iVertex, double *buffer);
         int getPointer(const std::string &type, size_t iVertex, void **return_p);
         gemmi::Structure* getStructure() const;
         size_t setCoord(const std::string &type, size_t iVertex, const CXXCoord<CXXCoord_ftype>&crd);
         void setScalar(size_t scalarHandle, size_t iVertex, double &value);
         void setScalar(const std::string name, size_t iVertex, double &value);

         int addPerVertexVector(const std::string name, double *vectorBuffer);
         int addPerVertexScalar(const std::string name, double *scalarBuffer);
         int addPerVertexPointer(const std::string name, void **pointerBuffer);

         size_t extendWithVectorData(size_t count, const std::string name, double *vectorBuffer);
         size_t updateWithVectorData(size_t count, const std::string name, size_t start, double *data);
         size_t updateWithPointerData(size_t count, const std::string name, size_t start, void **data);
         size_t extendTriangles(int *triangleBuffer, int count);

         size_t numberOfTriangles() const;
         size_t numberOfVertices() const;
         size_t vertex(size_t iTriangle, size_t iCorner) const;

         int upLoadSphere(CXXSphereElement &theSphere, double probeRadius, const int sense);
         int uploadTorus(CXXTorusElement &theTorus);
         int operator == (const CXXSurface &comparator) const { return (this == &comparator); }
         size_t getVectorHandle(const std::string name);
         size_t getReadVectorHandle(const std::string name);
         size_t getScalarHandle(const std::string name);
         size_t getReadScalarHandle(const std::string name);
         size_t getPointerHandle(const std::string name);
         int getScalar(int handle, int iVertex, double &result);
         int getVector(int handle, int iVertex, double *result);
         int addTriangle(const CXXTriangle &aTriangle);
         void compress(double tolerance);
         SurfaceParameters measuredProperties();

         // gemmi: assign each vertex its nearest selected atom (stored as a void* pointer property)
         int assignAtom(gemmi::Structure *structure, const std::set<const gemmi::Atom*> &selSet);
         int colorByAssignedAtom();
         void appendSurface(const CXXSurface &otherSurface);
      };
   }
}

#endif
