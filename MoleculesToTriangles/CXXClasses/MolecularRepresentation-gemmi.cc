/*
 * MoleculesToTriangles/CXXClasses/MolecularRepresentation-gemmi.cc
 *
 * gemmi-native twin of MolecularRepresentation.cpp (mmdb->gemmi migration).
 *
 * Ported so far: drawSpheres (BallsPrimitive is mmdb-free).
 * The bond/cylinder/stick/fan draw methods are stubbed pending gemmi twins of the
 * atom-storing primitives (BondsPrimitive / SticksPrimitive / CylindersPrimitive /
 * CylinderPoint / FlatFanPrimitive / BoxSectionPrimitive), which currently key on
 * mmdb::Atom*. Surfaces are Phase 5; colorByPotential is Phase 6.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include "MolecularRepresentation-gemmi.hh"
#include "BallsPrimitive.h"
#include "CylindersPrimitive-gemmi.hh"
#include "BoxSectionPrimitive-gemmi.hh"
#include "SticksPrimitive-gemmi.hh"
#include "FlatFanPrimitive-gemmi.hh"
#include "DiscreteSegment-gemmi.hh"
#include "DishyBase-gemmi.hh"
#include "ColorScheme-gemmi.hh"
#include "AtomPropertyRampColorRule-gemmi.hh"
#include "MoleculesToTriangles/CXXSurface/CXXUtils-gemmi.hh"
#include "SurfacePrimitive-gemmi.hh"
#include <gemmi/neighbor.hpp>
#include <gemmi/resinfo.hpp>

namespace {
   // map a stable atom serial (assigned by gemmi-bonds make_bonds over model 0)
   // back to its CRA, so a bond's two serials can be resolved to atoms.
   std::vector<gemmi::CRA> build_cra_index(gemmi::Model &model) {
      int n = 0;
      for (gemmi::Chain &c : model.chains)
         for (gemmi::Residue &r : c.residues)
            n += (int)r.atoms.size();
      std::vector<gemmi::CRA> idx(n, gemmi::CRA{nullptr, nullptr, nullptr});
      for (gemmi::Chain &c : model.chains)
         for (gemmi::Residue &r : c.residues)
            for (gemmi::Atom &a : r.atoms)
               if (a.serial >= 0 && a.serial < n)
                  idx[a.serial] = gemmi::CRA{&c, &r, &a};
      return idx;
   }
   inline double atom_distance(const gemmi::Atom *a, const gemmi::Atom *b) {
      double dx = a->pos.x - b->pos.x, dy = a->pos.y - b->pos.y, dz = a->pos.z - b->pos.z;
      return std::sqrt(dx*dx + dy*dy + dz*dz);
   }
}

std::string coot::m2t::MolecularRepresentation::renderStyles[] = {
   "Ribbon", "Sticks", "Calpha", "Spheres", "Cylinders",
   "MolecularSurface", "VdWSurface", "AccessibleSurface", "HydrogenBonds"
};

int coot::m2t::MolecularRepresentation::nStyles() {
   return sizeof(renderStyles) / sizeof(std::string);
}

// gemmi equivalent of the original (commented-out) prepareForSelectionInMMDB:
// for each ResidueNumber ramp rule, set its start/end value to the min/max
// residue seq-number of the residues its selection matches, so the colour ramp
// spans the selected range rather than the rule's default 1..1000.
void coot::m2t::MolecularRepresentation::prepareRampRules() {
   if (!colorScheme || !myMolecule) return;
   gemmi::Model *model = myMolecule->getModel();
   if (!model) return;
   for (auto &rule : colorScheme->getColorRules()) {
      auto ramp = std::dynamic_pointer_cast<AtomPropertyRampColorRule>(rule);
      if (!ramp) continue;
      if (ramp->getRampType() != AtomPropertyRampColorRule::ResidueNumber) continue;
      auto sel = ramp->getCompoundSelection();
      bool any = false;
      int lo = 0, hi = 0;
      for (gemmi::Chain &c : model->chains) {
         for (gemmi::Residue &r : c.residues) {
            if (r.atoms.empty() || !r.seqid.num) continue;
            // Only polymer residues are drawn by the ribbon, so only they should
            // set the ramp span - otherwise waters/ligands (often numbered in the
            // thousands) inflate endValue and squash the ramp into its start.
            gemmi::ResidueInfo ri = gemmi::find_tabulated_residue(r.name);
            if (!ri.is_amino_acid() && !ri.is_nucleic_acid()) continue;
            gemmi::CRA cra{&c, &r, &r.atoms[0]};
            if (sel && !sel->matches(cra)) continue;
            int n = *r.seqid.num;
            if (!any) { lo = hi = n; any = true; }
            else { if (n < lo) lo = n; if (n > hi) hi = n; }
         }
      }
      if (any && hi > lo) {
         ramp->setStartValue(float(lo));
         ramp->setEndValue(float(hi));
      }
   }
}

int coot::m2t::MolecularRepresentation::drawSpheres()
{
   gemmi::Model *model = myMolecule->getModel();
   if (!model) return 0;

   std::shared_ptr<BallsPrimitive> balls(new BallsPrimitive());
   int nAtoms = 0;
   float radiusMultiplier = floatParameters["ballsStyleRadiusMultiplier"];

   for (gemmi::Chain &chain : model->chains) {
      for (gemmi::Residue &residue : chain.residues) {
         for (gemmi::Atom &atom : residue.atoms) {
            gemmi::CRA cra{&chain, &residue, &atom};
            if (!selection->matches(cra)) continue;
            FCXXCoord atomColor = colorScheme->colorForAtom(cra);
            FCXXCoord atomCoord(atom.pos.x, atom.pos.y, atom.pos.z);
            float atomRadius = (float)coot::m2t::get_atom_radius(atom, residue) * radiusMultiplier;
            balls->addBall(atomCoord, atomColor, atomRadius);

            // index-array size limit - chunk balls (~100 per primitive)
            nAtoms++;
            if (nAtoms % 100 == 0) {
               displayPrimitives.push_back(balls);
               balls = std::make_shared<BallsPrimitive>();
            }
         }
      }
   }
   if (nAtoms % 100 != 0)
      displayPrimitives.push_back(balls);

   return 0;
}

// legacy (not in redraw() switch): the original used the line BondsPrimitive
int coot::m2t::MolecularRepresentation::drawBondsAsSticks()    { return 0; }

int coot::m2t::MolecularRepresentation::drawBondsAsCylinders()
{
   gemmi::Model *model = myMolecule->getModel();
   if (!model) return 0;
   float cylinderRadius = floatParameters["cylindersStyleCylinderRadius"];
   float ballRadius     = floatParameters["cylindersStyleBallRadius"];

   // a ball at each selected atom
   std::shared_ptr<BallsPrimitive> balls(new BallsPrimitive());
   int nAtoms = 0;
   for (gemmi::Chain &chain : model->chains)
      for (gemmi::Residue &residue : chain.residues)
         for (gemmi::Atom &atom : residue.atoms) {
            gemmi::CRA cra{&chain, &residue, &atom};
            if (!selection->matches(cra)) continue;
            FCXXCoord col = colorScheme->colorForAtom(cra);
            FCXXCoord coord(atom.pos.x, atom.pos.y, atom.pos.z);
            balls->addBall(coord, col, ballRadius);
            if (++nAtoms % 100 == 0) {
               displayPrimitives.push_back(balls);
               balls = std::make_shared<BallsPrimitive>();
            }
         }
   if (nAtoms % 100 != 0) displayPrimitives.push_back(balls);

   // a half-bond cylinder for each selected bond
   std::vector<gemmi::CRA> craIdx = build_cra_index(*model);
   int n = (int)craIdx.size();
   std::shared_ptr<CylindersPrimitive> cylinder(new CylindersPrimitive());
   cylinder->setAngularSampling(intParameters["cylindersStyleAngularSampling"]);
   int nBonds = 0;
   for (const coot::m2t::bond_t &b : myMolecule->bonds) {
      if (b.serial_1 < 0 || b.serial_1 >= n || b.serial_2 < 0 || b.serial_2 >= n) continue;
      gemmi::CRA &cra1 = craIdx[b.serial_1];
      gemmi::CRA &cra2 = craIdx[b.serial_2];
      if (!cra1.atom || !cra2.atom) continue;
      if (!selection->matches(cra1) || !selection->matches(cra2)) continue;
      // kludge (as original): cross-altconf atoms only bonded if close
      bool doContinue = (cra1.atom->altloc == '\0' && cra2.atom->altloc == '\0')
                        || atom_distance(cra1.atom, cra2.atom) < 1.9;
      if (!doContinue) continue;
      FCXXCoord c1(cra1.atom->pos.x, cra1.atom->pos.y, cra1.atom->pos.z);
      FCXXCoord c2(cra2.atom->pos.x, cra2.atom->pos.y, cra2.atom->pos.z);
      FCXXCoord col1 = colorScheme->colorForAtom(cra1);
      FCXXCoord col2 = colorScheme->colorForAtom(cra2);
      cylinder->addHalfAtomBondWithCoords(c1, col1, c2, col2, cylinderRadius);
      if (++nBonds % 1000 == 0) {
         displayPrimitives.push_back(cylinder);
         cylinder = std::make_shared<CylindersPrimitive>();
         cylinder->setAngularSampling(intParameters["cylindersStyleAngularSampling"]);
      }
   }
   if (nBonds % 1000 != 0) displayPrimitives.push_back(cylinder);
   return 0;
}

int coot::m2t::MolecularRepresentation::drawBondsAsNewSticks()
{
   gemmi::Model *model = myMolecule->getModel();
   if (!model) return 0;
   std::vector<gemmi::CRA> craIdx = build_cra_index(*model);
   int n = (int)craIdx.size();
   std::shared_ptr<SticksPrimitive> sticks(new SticksPrimitive());
   int nBonds = 0;
   for (const coot::m2t::bond_t &b : myMolecule->bonds) {
      if (b.serial_1 < 0 || b.serial_1 >= n || b.serial_2 < 0 || b.serial_2 >= n) continue;
      gemmi::CRA &cra1 = craIdx[b.serial_1];
      gemmi::CRA &cra2 = craIdx[b.serial_2];
      if (!cra1.atom || !cra2.atom) continue;
      if (!selection->matches(cra1) || !selection->matches(cra2)) continue;
      bool doAdd = (cra1.atom->altloc == '\0' && cra2.atom->altloc == '\0')
                   || atom_distance(cra1.atom, cra2.atom) < 1.9;
      if (!doAdd) continue;
      FCXXCoord c1(cra1.atom->pos.x, cra1.atom->pos.y, cra1.atom->pos.z);
      FCXXCoord c2(cra2.atom->pos.x, cra2.atom->pos.y, cra2.atom->pos.z);
      sticks->addAtom(cra1.atom->serial, c1, colorScheme->colorForAtom(cra1));
      sticks->addAtom(cra2.atom->serial, c2, colorScheme->colorForAtom(cra2));
      sticks->addBond(cra1.atom->serial, cra2.atom->serial);
      if (++nBonds % 1000 == 0) {
         displayPrimitives.push_back(sticks);
         sticks = std::make_shared<SticksPrimitive>();
      }
   }
   if (nBonds % 1000 != 0) displayPrimitives.push_back(sticks);
   return 0;
}

int coot::m2t::MolecularRepresentation::drawCalphas()
{
   for (DiscreteSegment *s : segments) delete s;
   segments.clear();
   myMolecule->identifySegments(segments, *selection);

   std::shared_ptr<SticksPrimitive> sticks(new SticksPrimitive());
   for (DiscreteSegment *segp : segments) {
      DiscreteSegment &segment = *segp;
      for (int i = 0; i < segment.nCalphas()-1; i++) {
         gemmi::CRA cra1 = segment.calpha(i);
         gemmi::CRA cra2 = segment.calpha(i+1);
         if (!cra1.atom || !cra2.atom) continue;
         FCXXCoord c1(cra1.atom->pos.x, cra1.atom->pos.y, cra1.atom->pos.z);
         FCXXCoord c2(cra2.atom->pos.x, cra2.atom->pos.y, cra2.atom->pos.z);
         sticks->addAtom(cra1.atom->serial, c1, colorScheme->colorForAtom(cra1));
         sticks->addAtom(cra2.atom->serial, c2, colorScheme->colorForAtom(cra2));
         sticks->addBond(cra1.atom->serial, cra2.atom->serial);
      }
   }
   displayPrimitives.push_back(sticks);
   return 0;
}

int coot::m2t::MolecularRepresentation::drawRibbon()
{
   constexpr int SSE_UNSET  = -32767;
   constexpr int SSE_NONE   = 0;
   constexpr int SSE_HELIX  = 1;
   constexpr int SSE_STRAND = 2;
   constexpr int SSE_DNARNA = 32767;

   // secondary structure from gemmi::Residue::flag ('H'/'E'/'L'), set by
   // coot::m2t::assign_secondary_structure; nucleic acids get the DNARNA path.
   auto sse_of = [](const gemmi::CRA &cra) -> int {
      if (!cra.residue) return SSE_NONE;
      if (gemmi::find_tabulated_residue(cra.residue->name).is_nucleic_acid()) return SSE_DNARNA;
      char f = cra.residue->flag;
      if (f == 'H') return SSE_HELIX;
      if (f == 'E') return SSE_STRAND;
      return SSE_NONE;
   };

   int subdivisionsPerCalpha = intParameters["ribbonStyleAxialSampling"];
   float stepPerSubdivision = 1./(float)subdivisionsPerCalpha;

   float radiusOneNone   = floatParameters["ribbonStyleCoilThickness"];
   float radiusTwoNone   = floatParameters["ribbonStyleCoilThickness"];
   float radiusOneHelix  = floatParameters["ribbonStyleHelixWidth"];
   float radiusTwoHelix  = floatParameters["ribbonStyleCoilThickness"];
   float radiusOneStrand = floatParameters["ribbonStyleStrandWidth"];
   float radiusTwoStrand = floatParameters["ribbonStyleCoilThickness"];
   float radiusOneDNARNA = floatParameters["ribbonStyleDNARNAWidth"];
   float radiusTwoDNARNA = floatParameters["ribbonStyleCoilThickness"];
   float radiusOneArrow  = floatParameters["ribbonStyleArrowWidth"];
   float radiusTwoArrow  = floatParameters["ribbonStyleCoilThickness"];

   for (DiscreteSegment *s : segments) delete s;
   segments.clear();
   myMolecule->identifySegments(segments, *selection);

   for (std::size_t iSegment = 0; iSegment < segments.size(); iSegment++) {
      DiscreteSegment &segment = *(segments[iSegment]);
      segment.evaluateNormals();
      if (boolParameters["smoothBetas"]) segment.smoothBetas();
      segment.evaluateSplines(subdivisionsPerCalpha);
      segment.evaluateColors(colorScheme);
      FCXXCoord color;

      std::shared_ptr<CylindersPrimitive> currentCylinder(new CylindersPrimitive());
      currentCylinder->setAngularSampling(intParameters["cylindersStyleAngularSampling"]);
      std::shared_ptr<BoxSectionPrimitive> currentBoxSection(new BoxSectionPrimitive());

      int lastSSE = SSE_UNSET;
      int currentSSE = SSE_UNSET;
      bool useArrowTipCoord = false;
      FCXXCoord arrowTipCoord;
      for (int iCalpha = 0; iCalpha < segment.nCalphas(); iCalpha++) {
         gemmi::CRA cra = segment.calpha(iCalpha);
         currentSSE = sse_of(cra);

         int nextSSE = SSE_NONE;
         if (iCalpha < segment.nCalphas()-1)
            nextSSE = sse_of(segment.calpha(iCalpha+1));

         if (lastSSE != currentSSE) {
            if (lastSSE == SSE_STRAND) {
               displayPrimitives.push_back(currentBoxSection);
               currentBoxSection = std::make_shared<BoxSectionPrimitive>();
            } else if (lastSSE != SSE_UNSET) {
               displayPrimitives.push_back(currentCylinder);
               currentCylinder = std::make_shared<CylindersPrimitive>();
               currentCylinder->setAngularSampling(intParameters["cylindersStyleAngularSampling"]);
            }
         }
         int endSubdivision = ((iCalpha == (segment.nCalphas()-1)) ? subdivisionsPerCalpha/2 : subdivisionsPerCalpha);
         int startSubdivision = ((iCalpha == 0) ? (subdivisionsPerCalpha/2) : 0);
         float xVal = 0.;

         bool inArrowTaper = false;
         FCXXCoord arrowBaseCoord, arrowBaseTangent, arrowBaseNormalOne, arrowBaseNormalTwo;
         float arrowBaseXVal = 0.f;

         for (int i = startSubdivision; i < endSubdivision; i++) {
            xVal = (iCalpha + (i*stepPerSubdivision)) - 0.5;
            FCXXCoord coord     = segment.coordFor(xVal);
            FCXXCoord normalOne = segment.normalOneFor(xVal);
            FCXXCoord normalTwo = segment.normalTwoFor(xVal);
            FCXXCoord tangentVec = segment.coordFor(xVal + stepPerSubdivision) - segment.coordFor(xVal - stepPerSubdivision);
            FCXXCoord tangent = tangentVec;
            tangent.normalise();
            normalOne = normalOne - tangent * (normalOne * tangent);
            normalOne.normalise();
            normalTwo = normalTwo - tangent * (normalTwo * tangent);
            normalTwo = normalTwo - normalOne * (normalTwo * normalOne);
            normalTwo.normalise();
            if (inArrowTaper) {
               coord = arrowBaseCoord + arrowBaseTangent * ((xVal - arrowBaseXVal) / (2.f * stepPerSubdivision));
               normalOne = arrowBaseNormalOne;
               normalTwo = arrowBaseNormalTwo;
            }
            if (useArrowTipCoord) { coord = arrowTipCoord; useArrowTipCoord = false; }
            color = segment.colorFor(xVal);
            float radiusMultiplier = (coord.r() > 0.0f) ? coord.r() : 1.0f;
            float radiusOne = radiusOneNone * radiusMultiplier;
            float radiusTwo = radiusTwoNone * radiusMultiplier;
            if (currentSSE == SSE_HELIX) {
               if (i < subdivisionsPerCalpha/2) {
                  if (currentSSE == lastSSE) { radiusOne = radiusOneHelix; radiusTwo = radiusTwoHelix; }
                  else {
                     float factor = (float)i / ((float)subdivisionsPerCalpha/2.f);
                     radiusOne = radiusOneNone + factor * (radiusOneHelix-radiusOneNone);
                     radiusTwo = radiusTwoNone + factor * (radiusTwoHelix-radiusTwoNone);
                  }
               } else {
                  if (currentSSE == nextSSE) { radiusOne = radiusOneHelix; radiusTwo = radiusTwoHelix; }
                  else {
                     float factor = (float)(i-(subdivisionsPerCalpha/2)) / ((float)subdivisionsPerCalpha/2.f);
                     radiusOne = radiusOneHelix - factor * (radiusOneHelix-radiusOneNone);
                     radiusTwo = radiusTwoHelix - factor * (radiusTwoHelix-radiusTwoNone);
                  }
               }
               radiusOne *= radiusMultiplier; radiusTwo *= radiusMultiplier;
               CylinderPoint cylinderPoint(coord, color, normalOne, normalTwo, radiusOne, radiusTwo);
               currentCylinder->addPoint(cylinderPoint);
            }
            else if (currentSSE == SSE_DNARNA) {
               if (i < subdivisionsPerCalpha/2) {
                  if (currentSSE == lastSSE) { radiusOne = radiusOneDNARNA; radiusTwo = radiusTwoDNARNA; }
                  else {
                     float factor = (float)i / ((float)subdivisionsPerCalpha/2.0f);
                     radiusOne = radiusOneNone + factor * (radiusOneDNARNA-radiusOneNone);
                     radiusTwo = radiusTwoNone + factor * (radiusTwoDNARNA-radiusTwoNone);
                  }
               } else {
                  if (currentSSE == nextSSE) { radiusOne = radiusOneDNARNA; radiusTwo = radiusTwoDNARNA; }
                  else {
                     float factor = (float)(i-(subdivisionsPerCalpha/2.0f)) / ((float)subdivisionsPerCalpha/2.0f);
                     radiusOne = radiusOneDNARNA - factor * (radiusOneDNARNA-radiusOneNone);
                     radiusTwo = radiusTwoDNARNA - factor * (radiusTwoDNARNA-radiusTwoNone);
                  }
               }
               radiusOne *= radiusMultiplier; radiusTwo *= radiusMultiplier;
               CylinderPoint cylinderPoint(coord, color, normalOne, normalTwo, radiusOne, radiusTwo);
               currentCylinder->addPoint(cylinderPoint);
            }
            else if (currentSSE == SSE_STRAND) {
               if (i == startSubdivision && lastSSE != currentSSE) {
                  CylinderPoint cylinderPoint(coord, color, normalOne, normalTwo,
                                              radiusOneNone * radiusMultiplier, radiusTwoNone * radiusMultiplier);
                  currentBoxSection->addPoint(cylinderPoint);
                  lastSSE = currentSSE;
               }
               if (i < subdivisionsPerCalpha/2) {
                  if (currentSSE == lastSSE) { radiusOne = radiusOneStrand; radiusTwo = radiusTwoStrand; }
               }
               if (i == subdivisionsPerCalpha/2) {
                  if (currentSSE != nextSSE) {
                     arrowBaseCoord = coord; arrowBaseTangent = tangentVec;
                     arrowBaseNormalOne = normalOne; arrowBaseNormalTwo = normalTwo;
                     arrowBaseXVal = xVal; inArrowTaper = true;
                     CylinderPoint cylinderPoint(coord, color, normalOne, normalTwo, radiusOneStrand * radiusMultiplier, radiusTwoStrand * radiusMultiplier);
                     currentBoxSection->addPoint(cylinderPoint);
                  }
               }
               if (i >= subdivisionsPerCalpha/2) {
                  if (currentSSE == nextSSE) { radiusOne = radiusOneStrand; radiusTwo = radiusTwoStrand; }
                  else {
                     float factor = (float)(i-(subdivisionsPerCalpha/2)) / ((float)subdivisionsPerCalpha/2.f);
                     radiusOne = radiusOneArrow - factor * (radiusOneArrow-radiusOneNone);
                     radiusTwo = radiusTwoArrow - factor * (radiusTwoArrow-radiusTwoNone);
                  }
               }
               radiusOne *= radiusMultiplier; radiusTwo *= radiusMultiplier;
               CylinderPoint cylinderPoint(coord, color, normalOne, normalTwo, radiusOne, radiusTwo);
               currentBoxSection->addPoint(cylinderPoint);
            }
            else {
               auto [ax, ay, az] = segment.anisoFor(xVal);
               float anisoMultiplierOne = std::sqrt(std::pow(normalOne.x()*ax,2) + std::pow(normalOne.y()*ay,2) + std::pow(normalOne.z()*az,2));
               float anisoMultiplierTwo = std::sqrt(std::pow(normalTwo.x()*ax,2) + std::pow(normalTwo.y()*ay,2) + std::pow(normalTwo.z()*az,2));
               radiusOne = radiusOneNone * anisoMultiplierOne;
               radiusTwo = radiusTwoNone * anisoMultiplierTwo;
               CylinderPoint cylinderPoint(coord, color, normalOne, normalTwo, radiusOne, radiusTwo);
               currentCylinder->addPoint(cylinderPoint);
            }
         }
         lastSSE = currentSSE;

         // extra point delimiting the end of the residue
         xVal += stepPerSubdivision;
         FCXXCoord coord = segment.coordFor(xVal);
         FCXXCoord normalOne = segment.normalOneFor(xVal);
         FCXXCoord normalTwo = segment.normalTwoFor(xVal);
         FCXXCoord tangent = segment.coordFor(xVal + stepPerSubdivision) - segment.coordFor(xVal - stepPerSubdivision);
         tangent.normalise();
         normalOne = normalOne - tangent * (normalOne * tangent);
         normalOne.normalise();
         normalTwo = normalTwo - tangent * (normalTwo * tangent);
         normalTwo = normalTwo - normalOne * (normalTwo * normalOne);
         normalTwo.normalise();
         if (inArrowTaper) {
            coord = arrowBaseCoord + arrowBaseTangent * ((xVal - arrowBaseXVal) / (2.f * stepPerSubdivision));
            normalOne = arrowBaseNormalOne; normalTwo = arrowBaseNormalTwo;
            arrowTipCoord = coord; useArrowTipCoord = true; inArrowTaper = false;
         }
         color = segment.colorFor(xVal);

         float radiusOne = radiusOneNone;
         float radiusTwo = radiusTwoNone;
         if (nextSSE == SSE_HELIX && lastSSE == SSE_HELIX) { radiusOne = radiusOneHelix; radiusTwo = radiusTwoHelix; }
         else if (nextSSE == SSE_DNARNA && lastSSE == SSE_DNARNA) { radiusOne = radiusOneDNARNA; radiusTwo = radiusTwoDNARNA; }
         else if (nextSSE == SSE_STRAND && lastSSE == SSE_STRAND) { radiusOne = radiusOneStrand; radiusTwo = radiusTwoStrand; }
         else {
            auto [ax, ay, az] = segment.anisoFor(xVal);
            float anisoMultiplierOne = std::sqrt(std::pow(normalOne.x()*ax,2) + std::pow(normalOne.y()*ay,2) + std::pow(normalOne.z()*az,2));
            float anisoMultiplierTwo = std::sqrt(std::pow(normalTwo.x()*ax,2) + std::pow(normalTwo.y()*ay,2) + std::pow(normalTwo.z()*az,2));
            radiusOne = radiusOneNone * anisoMultiplierOne;
            radiusTwo = radiusTwoNone * anisoMultiplierTwo;
         }
         CylinderPoint cylinderPoint(coord, color, normalOne, normalTwo, radiusOne, radiusTwo);
         if (currentSSE == SSE_STRAND) currentBoxSection->addPoint(cylinderPoint);
         else currentCylinder->addPoint(cylinderPoint);
      }
      if (currentSSE == SSE_STRAND) displayPrimitives.push_back(currentBoxSection);
      else if (currentSSE != SSE_UNSET) displayPrimitives.push_back(currentCylinder);
   }
   return 0;
}

int coot::m2t::MolecularRepresentation::drawDishyBases()
{
   std::map<gemmi::Chain *, DishyBaseContainer_t> dishy_bases_chain_map;
   myMolecule->identifyDishyBases(dishy_bases_chain_map, *selection);

   std::shared_ptr<BallsPrimitive> balls(new BallsPrimitive());
   std::shared_ptr<CylindersPrimitive> cylinder(new CylindersPrimitive());
   cylinder->setAngularSampling(intParameters["dishStyleAngularSampling"]);
   float cylinderRadius = floatParameters["cylindersStyleCylinderRadius"];
   float ballRadius     = floatParameters["cylindersStyleBallRadius"];

   for (auto &kv : dishy_bases_chain_map) {
      for (DishyBase_t &db : kv.second.bases) {
         gemmi::CRA atom1 = db.ribose_atoms[1];
         if (!atom1.atom) continue;
         FCXXCoord atom1Color = colorScheme->colorForAtom(atom1);

         // the dish itself
         float dishRadius = (float)db.radius;
         float dishThick  = 0.14f * dishRadius;
         balls->addBall(db.centre, atom1Color, dishRadius, db.normal, dishThick);
         if (balls->getBalls().size() % 100 == 0) {
            displayPrimitives.push_back(balls);
            balls = std::make_shared<BallsPrimitive>();
         }

         // flat fan over the ribose ring
         std::vector<FCXXCoord> fanPositions;
         for (gemmi::CRA &cra : db.ribose_atoms)
            if (cra.atom)
               fanPositions.push_back(FCXXCoord(cra.atom->pos.x, cra.atom->pos.y, cra.atom->pos.z));
         auto flatFan = std::make_shared<FlatFanPrimitive>(fanPositions, atom1Color);
         displayPrimitives.push_back(flatFan);

         // balls at the ribose atoms
         for (gemmi::CRA &cra : db.ribose_atoms) {
            if (!cra.atom) continue;
            FCXXCoord coord(cra.atom->pos.x, cra.atom->pos.y, cra.atom->pos.z);
            FCXXCoord atomColor = colorScheme->colorForAtom(cra);
            balls->addBall(coord, atomColor, ballRadius);
            if (balls->getBalls().size() % 100 == 0) {
               displayPrimitives.push_back(balls);
               balls = std::make_shared<BallsPrimitive>();
            }
         }

         // ribose ring bonds
         for (const std::pair<int,int> &bond : DishyBase_t::bondingPattern) {
            gemmi::CRA a1 = db.ribose_atoms[bond.first];
            gemmi::CRA a2 = db.ribose_atoms[bond.second];
            if (!a1.atom || !a2.atom) continue;
            FCXXCoord c1(a1.atom->pos.x, a1.atom->pos.y, a1.atom->pos.z);
            FCXXCoord c2(a2.atom->pos.x, a2.atom->pos.y, a2.atom->pos.z);
            FCXXCoord col1 = colorScheme->colorForAtom(a1);
            FCXXCoord col2 = colorScheme->colorForAtom(a2);
            cylinder->addHalfAtomBondWithCoords(c1, col1, c2, col2, cylinderRadius);
         }

         // stick from C1' a third of the way to the base centre
         FCXXCoord atom1Coord(atom1.atom->pos.x, atom1.atom->pos.y, atom1.atom->pos.z);
         FCXXCoord basePseudoAtomPosition = atom1Coord + (db.centre - atom1Coord) / 3.;
         cylinder->addHalfAtomBondWithCoords(atom1Coord, atom1Color, basePseudoAtomPosition, atom1Color, cylinderRadius);
      }
   }
   displayPrimitives.push_back(cylinder);
   if (balls->getBalls().size() % 100 != 0)
      displayPrimitives.push_back(balls);
   return 0;
}

int coot::m2t::MolecularRepresentation::drawStickBases()
{
   gemmi::Model *model = myMolecule->getModel();
   if (!model) return 0;

   auto is_nucleic_acid = [](const std::string &n) {
      return n=="G"||n=="A"||n=="T"||n=="C"||n=="U"||n=="DG"||n=="DA"||n=="DC"||n=="DT";
   };
   auto atom_name_pair = [](const std::string &n) -> std::pair<std::string,std::string> {
      if (n=="DG"||n=="DA"||n=="G"||n=="A") return {"C3'", "N1"};
      if (n=="DT"||n=="U")                  return {"C3'", "O4"};
      if (n=="DC"||n=="C")                  return {"C3'", "N4"};
      return {"", ""};
   };

   std::shared_ptr<CylindersPrimitive> cylinder(new CylindersPrimitive());
   cylinder->setAngularSampling(intParameters["cylindersStyleAngularSampling"]);
   float cylinderRadius = floatParameters["cylindersStyleCylinderRadius"];
   float ballRadius     = floatParameters["cylindersStyleBallRadius"];
   std::shared_ptr<BallsPrimitive> balls(new BallsPrimitive());

   for (gemmi::Chain &chain : model->chains) {
      for (gemmi::Residue &residue : chain.residues) {
         if (!is_nucleic_acid(residue.name)) continue;
         auto names = atom_name_pair(residue.name);
         if (names.first.empty()) continue;
         gemmi::Atom *atom_1 = nullptr, *atom_2 = nullptr;
         for (gemmi::Atom &atom : residue.atoms) {
            if (!atom_1 && atom.name == names.first)  atom_1 = &atom;
            if (!atom_2 && atom.name == names.second) atom_2 = &atom;
         }
         if (!atom_1 || !atom_2) continue;
         gemmi::CRA cra1{&chain, &residue, atom_1};
         gemmi::CRA cra2{&chain, &residue, atom_2};
         if (!selection->matches(cra1) || !selection->matches(cra2)) continue;
         FCXXCoord atom1Color = colorScheme->colorForAtom(cra1);
         FCXXCoord atom2Color = colorScheme->colorForAtom(cra2);
         FCXXCoord c1(atom_1->pos.x, atom_1->pos.y, atom_1->pos.z);
         FCXXCoord c2(atom_2->pos.x, atom_2->pos.y, atom_2->pos.z);
         cylinder->addHalfAtomBondWithCoords(c1, atom1Color, c2, atom2Color, cylinderRadius);
         FCXXCoord ballCoord(atom_2->pos.x, atom_2->pos.y, atom_2->pos.z);
         balls->addBall(ballCoord, atom1Color, ballRadius);
      }
   }
   displayPrimitives.push_back(cylinder);
   displayPrimitives.push_back(balls);
   return 0;
}
int coot::m2t::MolecularRepresentation::drawHydrogenBonds()
{
   gemmi::Model *model = myMolecule->getModel();
   if (!model) return 0;

   // Matches the original: rMax comes from a parameter with no installed default,
   // so an unset (0) value means "draw nothing".
   float hbondDistance = floatParameters[std::string("hBondsStyleRMax")];
   if (hbondDistance <= 0.0f) return 0;

   double search_radius = hbondDistance > 4.0 ? hbondDistance : 4.0;
   gemmi::NeighborSearch ns(*model, myMolecule->getStructure().cell, search_radius);
   ns.populate(true);

   std::shared_ptr<BallsPrimitive> balls(new BallsPrimitive());
   int nDots = 0;

   for (gemmi::Chain &chain : model->chains) {
      for (gemmi::Residue &residue : chain.residues) {
         for (gemmi::Atom &atom : residue.atoms) {
            if (atom.element == gemmi::El::C) continue;           // non-carbon only
            gemmi::CRA cra1{&chain, &residue, &atom};
            if (!selection->matches(cra1)) continue;

            std::vector<gemmi::NeighborSearch::Mark*> marks = ns.find_neighbors(atom, 1.0, hbondDistance);
            for (gemmi::NeighborSearch::Mark *m : marks) {
               gemmi::CRA cra2 = m->to_cra(*model);
               if (!cra2.atom) continue;
               if (cra2.atom->serial <= atom.serial) continue;   // each pair once
               if (cra2.atom->element == gemmi::El::C) continue;
               if (!selection->matches(cra2)) continue;

               int rn1 = residue.seqid.num ? *residue.seqid.num : 0;
               int rn2 = (cra2.residue && cra2.residue->seqid.num) ? *cra2.residue->seqid.num : 0;
               std::string n1 = atom.name, n2 = cra2.atom->name;
               bool aa1 = gemmi::find_tabulated_residue(residue.name).is_amino_acid();
               bool aa2 = cra2.residue && gemmi::find_tabulated_residue(cra2.residue->name).is_amino_acid();
               bool name_pair = (n1=="N"&&n2=="N") || (n1=="O"&&n2=="O") ||
                                (n1=="N"&&n2=="O") || (n1=="O"&&n2=="N");
               // faithful to the original (draws unless it is a close backbone N/O pair)
               bool draw = (std::abs(rn1 - rn2) > 2) || !aa1 || !aa2 || !name_pair;
               if (!draw) continue;

               FCXXCoord c1(atom.pos.x, atom.pos.y, atom.pos.z);
               FCXXCoord c2(cra2.atom->pos.x, cra2.atom->pos.y, cra2.atom->pos.z);
               FCXXCoord diff = c2 - c1;
               for (int iStep = 1; iStep < 7; iStep++) {
                  float step = (float)iStep / 8.;
                  FCXXCoord ballCoord = c1 + diff*step;
                  FCXXCoord atomColor = (iStep > 5) ? colorScheme->colorForAtom(cra2)
                                                    : colorScheme->colorForAtom(cra1);
                  balls->addBall(ballCoord, atomColor, 0.1);
                  nDots++;
                  if (nDots % 100 == 0) {
                     displayPrimitives.push_back(balls);
                     balls = std::make_shared<BallsPrimitive>();
                  }
               }
            }
         }
      }
   }
   if (nDots % 100 != 0)
      displayPrimitives.push_back(balls);

   return 0;
}

// ---- Surfaces (CXXSurface gemmi twin) ----
int coot::m2t::MolecularRepresentation::drawSurfaceOfKind(int surfaceKind)
{
   gemmi::Model *model = myMolecule->getModel();
   if (!model) return 0;
   gemmi::Structure *structure = &myMolecule->getStructure();

   float probeRadius = floatParameters["surfaceStyleProbeRadius"];
   float radiusMultiplier = floatParameters["ballsStyleRadiusMultiplier"];
   if (probeRadius <= 0.f) probeRadius = 1.4f;
   if (radiusMultiplier <= 0.f) radiusMultiplier = 1.0f;

   // The full selection is the context; the per-chunk subset is "central".
   std::set<const gemmi::Atom*> fullSet;
   std::vector<const gemmi::Atom*> selAtoms;
   for (gemmi::Chain &chain : model->chains)
      for (gemmi::Residue &res : chain.residues)
         for (gemmi::Atom &atom : res.atoms) {
            gemmi::CRA cra{&chain, &res, &atom};
            if (selection->matches(cra)) { fullSet.insert(&atom); selAtoms.push_back(&atom); }
         }

   // Chop into chunks if GLIndexType is only a short (vertex index limit).
   int chunkSize = (sizeof(GLIndexType) == sizeof(short)) ? 100 : 1000000;
   for (size_t start = 0; start < selAtoms.size(); start += chunkSize) {
      std::set<const gemmi::Atom*> chunkSet;
      size_t end = std::min(selAtoms.size(), start + (size_t)chunkSize);
      for (size_t i = start; i < end; i++) chunkSet.insert(selAtoms[i]);
      auto surfacePrimitive = std::make_shared<coot::m2t::SurfacePrimitive>(
         structure, *model, chunkSet, fullSet, colorScheme,
         (coot::m2t::SurfacePrimitive::SurfaceType) surfaceKind, probeRadius, radiusMultiplier);
      if (surfacePrimitive->getCXXSurfaceMaker()) displayPrimitives.push_back(surfacePrimitive);
   }
   return 0;
}
int coot::m2t::MolecularRepresentation::drawMolecularSurface() { return drawSurfaceOfKind(coot::m2t::SurfacePrimitive::MolecularSurface); }
int coot::m2t::MolecularRepresentation::drawVdWSurface()       { return drawSurfaceOfKind(coot::m2t::SurfacePrimitive::VdWSurface); }
int coot::m2t::MolecularRepresentation::drawAccessibleSurface(){ return drawSurfaceOfKind(coot::m2t::SurfacePrimitive::AccessibleSurface); }
