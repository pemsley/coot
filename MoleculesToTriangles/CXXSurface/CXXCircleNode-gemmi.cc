/*
 * MoleculesToTriangles/CXXSurface/CXXCircleNode-gemmi.cc
 *
 * gemmi-native twin of CXXCircleNode.cpp.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */
#include "CXXSurfaceVertex.h"
#include "CXXCircleNode-gemmi.hh"
#include "CXXCoord.h"
#include "CXXCircle-gemmi.hh"
#include "CXXNewHood-gemmi.hh"
#include <cmath>
#include <iostream>

using namespace coot::m2t;

CXXCircleNode::CXXCircleNode() :
   theParent(0), theOtherCircle(0),
   theCoord(CXXCoord<CXXCoord_ftype>(0.,0.,0.)),
   unitRadius(CXXCoord<CXXCoord_ftype>(0.,0.,0.)),
   theAngle(0.), theFlag(-1), thisIsDeleted(0),
   atomI(0), atomJ(0), atomK(0) {}

CXXCircleNode::CXXCircleNode(const CXXCircle *aParent, const CXXCircle *anOtherCircle, const CXXCoord<CXXCoord_ftype>&crd, int aFlag) :
   theParent(aParent), theOtherCircle(anOtherCircle),
   theCoord(crd), unitRadius(crd - aParent->getCentreOfCircle()),
   theAngle(0.), theFlag(aFlag), thisIsDeleted(0),
   atomI(0), atomJ(0), atomK(0)
{
   unitRadius /= aParent->getRadiusOfCircle();
   if (theParent) {
      atomJ = theParent->getAtomJ();
      if (theParent->getParent()) {
         atomI = theParent->getParent()->getAtomI();
      }
   }
   if (theOtherCircle) atomK = theOtherCircle->getAtomJ();
}

bool CXXCircleNode::shouldDelete(const CXXCircleNode &aNode) { return (aNode.isDeleted() == 1); }

void CXXCircleNode::setParent(CXXCircle *parent) {
   theParent = parent;
   if (theParent) {
      atomJ = theParent->getAtomJ();
      if (theParent->getParent()) {
         atomI = theParent->getParent()->getAtomI();
      }
   }
}

void CXXCircleNode::setOtherCircle(CXXCircle *parent) {
   theOtherCircle = parent;
   if (theOtherCircle) atomK = theOtherCircle->getAtomJ();
}

int CXXCircleNode::setReference(const CXXCoord<CXXCoord_ftype>&referenceVector) {
   theAngle = theParent->getNormal().angleBetween(getUnitRadius(), referenceVector);
   return 0;
}

void CXXCircleNode::setCoord(const CXXCoord<CXXCoord_ftype>&coord) {
   theCoord = coord;
   unitRadius = theCoord - theParent->getCentreOfCircle();
   unitRadius /= theParent->getRadiusOfCircle();
}

int CXXCircleNode::probeContacts(std::vector<CXXCircleNode> &probes,
                                 double probeRadius,
                                 std::map<CXXCircleNode *, std::vector<CXXCircleNode *> > &contactMap)
{
   if (probes.size() == 0) return 1;
   double probeRadiusX2 = 2.*probeRadius;
   double probeRadiusX2Sq = probeRadiusX2*probeRadiusX2;

   double limits[3][2];
   for (int i=0; i<3; i++) { limits[i][0] = 1e30; limits[i][1] = -1e30; }
   std::vector<CXXCircleNode>::iterator end = probes.end();
   for (std::vector<CXXCircleNode>::iterator probeIter = probes.begin(); probeIter!=end; ++probeIter) {
      for (int i=0; i<3; i++) {
         limits[i][0] = (limits[i][0]<(*probeIter)[i]?limits[i][0]:(*probeIter)[i]);
         limits[i][1] = (limits[i][1]>(*probeIter)[i]?limits[i][1]:(*probeIter)[i]);
      }
   }
   for (int i=0; i<3; i++) { limits[i][0] -= 0.00000001; limits[i][1] += 0.00000001; }

   double minimumBinSize = probeRadius;
   int nBins[3];
   double binWidth[3];
   for (int i=0; i<3; i++) {
      int nMinimumBins = ceil((limits[i][1]-limits[i][0])/minimumBinSize);
      if (nMinimumBins<10) { nBins[i] = nMinimumBins; binWidth[i] = minimumBinSize; }
      else { nBins[i] = 10; binWidth[i] = (limits[i][1]-limits[i][0])/10.; }
   }

   std::vector<std::vector<std::vector<std::vector<std::vector<CXXCircleNode>::iterator> > > >binnedProbes;
   binnedProbes.resize(nBins[0]);
   for (int i=0; i<nBins[0]; i++) {
      binnedProbes[i].resize(nBins[1]);
      for (int j=0; j<nBins[1]; j++) binnedProbes[i][j].resize(nBins[2]);
   }

   for (std::vector<CXXCircleNode>::iterator probeIter = probes.begin(); probeIter!=end; ++probeIter) {
      int iBin[3];
      for (int i=0; i<3; i++) iBin[i] = floor((*probeIter)[i]-limits[i][0]) / binWidth[i];
      binnedProbes[iBin[0]][iBin[1]][iBin[2]].push_back(probeIter);
      contactMap[&(*probeIter)].resize(0);
   }

   for (int iBinX = 0; iBinX < nBins[0]; iBinX++) {
      int startBinX = std::max(0,iBinX-1);
      int endBinX = std::min(nBins[0],iBinX+2);
      for (int iBinY = 0; iBinY < nBins[1]; iBinY++) {
         int startBinY = std::max(0,iBinY-1);
         int endBinY = std::min(nBins[1],iBinY+2);
         for (int iBinZ = 0; iBinZ < nBins[2]; iBinZ++) {
            int startBinZ = std::max(0,iBinZ-1);
            int endBinZ = std::min(nBins[2],iBinZ+2);
            std::vector<std::vector<CXXCircleNode>::iterator> &centralBin(binnedProbes[iBinX][iBinY][iBinZ]);
            for (int iCentralProbe = 0; iCentralProbe < (int)centralBin.size(); iCentralProbe++) {
               CXXCircleNode &centralProbeRef(*centralBin[iCentralProbe]);
               std::vector<const gemmi::Atom*> ijkCentral(3);
               std::vector<const gemmi::Atom*> ijkOther(3);
               ijkCentral[0] = centralProbeRef.getAtomI();
               ijkCentral[1] = centralProbeRef.getAtomJ();
               ijkCentral[2] = centralProbeRef.getAtomK();
               std::sort(ijkCentral.begin(), ijkCentral.end());
               for (int searchBinX = startBinX; searchBinX<endBinX; searchBinX++) {
                  for (int searchBinY = startBinY; searchBinY<endBinY; searchBinY++) {
                     for (int searchBinZ = startBinZ; searchBinZ<endBinZ; searchBinZ++) {
                        std::vector<std::vector<CXXCircleNode>::iterator> &searchBin(binnedProbes[searchBinX][searchBinY][searchBinZ]);
                        std::vector<std::vector<CXXCircleNode>::iterator>::iterator binSearchEnd = searchBin.end();
                        for (std::vector<std::vector<CXXCircleNode>::iterator>::iterator otherProbe = searchBin.begin();
                             otherProbe!=binSearchEnd; ++otherProbe) {
                           CXXCircleNode &otherProbeRef(**otherProbe);
                           if (&centralProbeRef != &otherProbeRef) {
                              ijkOther[0] = otherProbeRef.getAtomI();
                              ijkOther[1] = otherProbeRef.getAtomJ();
                              ijkOther[2] = otherProbeRef.getAtomK();
                              std::sort(ijkOther.begin(), ijkOther.end());
                              if (ijkCentral[0] != ijkOther[0] || ijkCentral[1] != ijkOther[1] || ijkCentral[2] != ijkOther[2]
                                  || !centralProbeRef.getCoord().isNearly(otherProbeRef.getCoord(), 0.00001)) {
                                 if (centralProbeRef.getCoord().isNearly(otherProbeRef.getCoord(), probeRadiusX2)) {
                                    CXXCoord<CXXCoord_ftype>diff(centralProbeRef.getCoord() - otherProbeRef.getCoord());
                                    if (diff.get3DLengthSq()<probeRadiusX2Sq) {
                                       contactMap[&centralProbeRef].push_back(&otherProbeRef);
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
   std::map<CXXCircleNode *, std::vector<CXXCircleNode *> >::iterator endIter = contactMap.end();
   int nContacts = 0;
   for (std::map<CXXCircleNode *, std::vector<CXXCircleNode *> >::iterator iter = contactMap.begin(); iter !=endIter; ++iter)
      nContacts += iter->second.size();
   std::cout << "Total of " << nContacts << " contacts among the reentrant Probes\n";
   return 0;
}

bool CXXCircleNode::shouldDeletePointer(CXXCircleNode* &aNodePointer) { return (aNodePointer->isDeleted()==1); }

bool CXXCircleNode::equals(CXXCircleNode &node1, CXXCircleNode &node2) {
   std::vector<int>ijkCentral(3);
   std::vector<int>ijkOther(3);
   ijkCentral[0] = node1.getAtomI()->serial;
   ijkCentral[1] = node1.getAtomJ()->serial;
   ijkCentral[2] = node1.getAtomK()->serial;
   std::sort(ijkCentral.begin(), ijkCentral.end());
   ijkOther[0] = node2.getAtomI()->serial;
   ijkOther[1] = node2.getAtomJ()->serial;
   ijkOther[2] = node2.getAtomK()->serial;
   std::sort(ijkOther.begin(), ijkOther.end());
   if (ijkCentral[0] != ijkOther[0]) return false;
   if (ijkCentral[1] != ijkOther[1]) return false;
   if (ijkCentral[2] != ijkOther[2]) return false;
   if (!node1.getCoord().isNearly(node2.getCoord(), 0.0001)) return false;
   return true;
}

bool CXXCircleNode::equalsPntr(CXXCircleNode* &node1, CXXCircleNode* &node2) {
   std::vector<const gemmi::Atom*>ijkCentral(3);
   std::vector<const gemmi::Atom*>ijkOther(3);
   ijkCentral[0] = node1->getAtomI();
   ijkCentral[1] = node1->getAtomJ();
   ijkCentral[2] = node1->getAtomK();
   std::sort(ijkCentral.begin(), ijkCentral.end());
   ijkOther[0] = node2->getAtomI();
   ijkOther[1] = node2->getAtomJ();
   ijkOther[2] = node2->getAtomK();
   std::sort(ijkOther.begin(), ijkOther.end());
   if (ijkCentral[0] != ijkOther[0]) return false;
   if (ijkCentral[1] != ijkOther[1]) return false;
   if (ijkCentral[2] != ijkOther[2]) return false;
   if (!node1->getCoord().isNearly(node2->getCoord(), 0.00001)) return false;
   return true;
}

void CXXCircleNode::filterContacts(std::map<CXXCircleNode *, std::vector<CXXCircleNode *> > &contactMap) {
   std::map< CXXCircleNode *, std::vector<CXXCircleNode *> >::iterator contactMapEnd(contactMap.end());
   for (std::map< CXXCircleNode *, std::vector<CXXCircleNode *> >::iterator contactMapIter = contactMap.begin();
        contactMapIter != contactMapEnd; ++contactMapIter) {
      std::vector<CXXCircleNode *> &neighbours(contactMapIter->second);
      std::vector<CXXCircleNode *>::iterator neighboursEnd(neighbours.end());
      std::vector<CXXCircleNode *>::iterator neighboursBegin(neighbours.begin());
      neighbours.erase(std::unique(neighboursBegin, neighboursEnd, CXXCircleNode::equalsPntr), neighboursEnd);
   }
}

bool CXXCircleNode::angleLessThan(const CXXCircleNode &node1, const CXXCircleNode &node2) {
   return node1.getAngle() < node2.getAngle();
}
