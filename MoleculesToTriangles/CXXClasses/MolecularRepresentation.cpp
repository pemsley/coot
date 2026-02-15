/*
 * MoleculesToTriangles/CXXClasses/MolecularRepresentation.cpp
 *
 * Copyright 2009 by Martin Noble, University of Oxford
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

#include <cstring>

#include "MolecularRepresentation.h"
#include "DishyBase.h"
#include "SpherePrimitive.h"
#include "CylindersPrimitive.h"
#include "BoxSectionPrimitive.h"
#include "BallsPrimitive.h"
//#include "LinesPrimitive.h"
#include "SticksPrimitive.h"
#include "BondsPrimitive.h"
#include "SurfacePrimitive.h"
#include "FlatFanPrimitive.h"
#include "Renderer.h"
#include "MoleculesToTriangles/CXXSurface/CXXUtils.h"
#include "MoleculesToTriangles/CXXSurface/CXXSurfaceVertex.h"

std::string MolecularRepresentation::renderStyles[] = {"Ribbon", "Sticks", "Calpha", "Spheres", "Cylinders", "MolecularSurface", "VdWSurface", "AccessibleSurface", "HydrogenBonds"};

int MolecularRepresentation::nStyles() {
    return sizeof(renderStyles)/ sizeof(std::string);
};

int MolecularRepresentation::drawSpheres()
{
    mmdb::Manager *mmdb = myMolecule->getMmdb();
    //selection->describe();
    selHnd = selection->handleInMMDB(mmdb);
    auto handles = colorScheme->prepareForMMDB(mmdb);

    mmdb::Atom** SelAtom;
    int nSelAtoms;
    mmdb->GetSelIndex(selHnd, SelAtom, nSelAtoms);
    std::shared_ptr<BallsPrimitive>balls(new BallsPrimitive());
    int nAtoms = 0;
    float radiusMultiplier = floatParameters["ballsStyleRadiusMultiplier"];

    for (int i=0; i<nSelAtoms; i++){
        mmdb::Atom* atom1 = SelAtom[i];
        FCXXCoord atom1Color = colorScheme->colorForAtom(atom1, handles);
        FCXXCoord atom1Coord(atom1->x,atom1->y, atom1->z);
        float atomRadius = CXXUtils::getAtomRadius(mmdb, atom1) * radiusMultiplier;

        balls->addBall(atom1Coord, atom1Color, atomRadius);

        //Because we might be restricted to GL_SHORT_INTs in the index array (don't ask !)
        //we cannot reliably put all of the balls  into the same primitive: sensible limit
        //would be about 1000

        nAtoms ++;
        if (nAtoms%100 == 0) {
            displayPrimitives.push_back(balls);
            std::shared_ptr<BallsPrimitive> balls(new BallsPrimitive());
        }
    }
    if (nAtoms%100 != 0) {
        displayPrimitives.push_back(balls);
    }
    colorScheme->freeSelectionHandles(mmdb, handles);
    mmdb->DeleteSelection(selHnd);

    return 0;
}

int MolecularRepresentation::drawBondsAsSticks()
{
    mmdb::Manager *mmdb = myMolecule->getMmdb();
    //selection->describe();
    selHnd = selection->handleInMMDB(mmdb);
    std::map<std::shared_ptr<ColorRule>,int>handles = colorScheme->prepareForMMDB(mmdb);

    mmdb::Atom** selAtoms;
    int nSelAtoms;

    std::shared_ptr<BondsPrimitive>bonds(new BondsPrimitive());
    displayPrimitives.push_back(bonds);
    bonds->setColorScheme(colorScheme);

    mmdb->GetSelIndex(selHnd, selAtoms, nSelAtoms);
    for (int iAtom = 0; iAtom < nSelAtoms; iAtom++){
        mmdb::Atom* atom1 = selAtoms[iAtom];

        int nBondedAtoms;
        mmdb::AtomBond* bondedAtoms;
        atom1->GetBonds(bondedAtoms, nBondedAtoms);
        for (int iOtherAtom = 0; iOtherAtom < nBondedAtoms; iOtherAtom++){
            mmdb::Atom* atom2 = bondedAtoms[iOtherAtom].atom;
            if (atom2->isInSelection(selHnd)){
                bonds->addPair(atom1, atom2);
            }
        }
    }
    bonds->evaluateGLPrimitives(handles);
    colorScheme->freeSelectionHandles(mmdb, handles);
    mmdb->DeleteSelection(selHnd);
    return 0;
}

int MolecularRepresentation::drawBondsAsCylinders()
{
    mmdb::Manager *mmdb = myMolecule->getMmdb();
    //selection->describe();
    selHnd = selection->handleInMMDB(mmdb);
    std::map<std::shared_ptr<ColorRule>,int>handles = colorScheme->prepareForMMDB(mmdb);

    FCXXCoord xAxis (1., 0., 0., 0.);
    FCXXCoord yAxis (0., 1., 0., 0.);
    FCXXCoord zAxis (0., 0., 1., 0.);
    
    float cylinderRadius = floatParameters[std::string("cylindersStyleCylinderRadius")];
    float ballRadius = floatParameters[std::string("cylindersStyleBallRadius")];
    
    int nBonds = 0;
    shared_ptr<CylindersPrimitive>cylinder(new CylindersPrimitive());
    cylinder->setAngularSampling(intParameters["cylindersStyleAngularSampling"]);
    int nAtoms = 0;
    std::shared_ptr<BallsPrimitive>balls(new BallsPrimitive());
    
    for (int iAtom = 1; iAtom <= myMolecule->getMmdb()->GetNumberOfAtoms(); iAtom++){
        mmdb::Atom* atom1 = myMolecule->getMmdb()->GetAtomI(iAtom);
        
        if (atom1->isInSelection(selHnd)){
            FCXXCoord atom1Color = colorScheme->colorForAtom(atom1, handles);
            
            FCXXCoord atom1Coord(atom1->x,atom1->y, atom1->z);
            balls->addBall(atom1Coord, atom1Color, ballRadius);
            
            //Because we might be restricted to GL_SHORT_INTs in the index array (don't ask !)
            //we cannot reliably put all of the balls  into the same primitive: sensible limit
            //would be about 1000
            
            nAtoms ++;
            if (nAtoms%100 == 0) {
                displayPrimitives.push_back(balls);
                balls = shared_ptr<BallsPrimitive>(new BallsPrimitive());
            }
            
            //displayPrimitives.push_back(new SpherePrimitive(FCXXCoord (atom1->x,atom1->y, atom1->z), 0.3, atom1Color));
            
            int nBondedAtoms;
            mmdb::AtomBond* bondedAtoms;
            atom1->GetBonds(bondedAtoms, nBondedAtoms);
            FCXXCoord coord1(atom1->x, atom1->y, atom1->z);
            for (int iOtherAtom = 0; iOtherAtom < nBondedAtoms; iOtherAtom++){
                mmdb::Atom* atom2 = bondedAtoms[iOtherAtom].atom;
                if (atom2->GetIndex()>iAtom && atom2->isInSelection(selHnd)){
                    
                    //Nasty kludge here...MMDB's MakeBonds screws up where there are multiple conformations
                    bool doContinue = false;
                    if (!strcmp(atom1->altLoc,"") &&
                        !strcmp(atom2->altLoc,"")) {
                        doContinue = true;
                    }
                    else {
                        float dx = atom1->x - atom2->x;
                        float dy = atom1->y - atom2->y;
                        float dz = atom1->z - atom2->z;
                        float dist = sqrtf (dx*dx + dy*dy + dz*dz);
                        if (dist < 1.9f) doContinue = true;
                    }
                    
                    if (doContinue){
                        FCXXCoord atom2Color =  colorScheme->colorForAtom(atom2, handles);
                        cylinder->addHalfAtomBond(atom1, atom1Color, atom2, atom2Color, cylinderRadius);
                        
                        //Because we might be restricted to GL_SHORT_INTs in the index array (don't ask !)
                        //we cannot reliably put all of the bonds  into the same primitive: sensible limit
                        //would be about 1000
                        
                        nBonds ++;
                        if (nBonds%1000 == 0) {
                            displayPrimitives.push_back(cylinder);
                            cylinder = shared_ptr<CylindersPrimitive>(new CylindersPrimitive());
                            cylinder->setAngularSampling(intParameters["cylindersStyleAngularSampling"]);
                        }
                    }
                }
            }
        }
    }
    if (nAtoms%100 != 0) {
        displayPrimitives.push_back(balls);
    }
    if (nBonds%1000 != 0) {
        displayPrimitives.push_back(cylinder);
    }
    
    colorScheme->freeSelectionHandles(mmdb, handles);
    mmdb->DeleteSelection(selHnd);
	
    return 0;		
}

int MolecularRepresentation::drawHydrogenBonds()
{
    mmdb::Manager *mmdb = myMolecule->getMmdb();
	//selection->describe();
	selHnd = selection->handleInMMDB(mmdb);
    std::map<std::shared_ptr<ColorRule>,int>handles = colorScheme->prepareForMMDB(mmdb);
    const char *allText= "/*/*/*/[!C]";
    mmdb->Select(selHnd, mmdb::STYPE_ATOM, *allText, mmdb::SKEY_AND);
    int nSelAtoms;
    mmdb::Atom** selAtoms;
    mmdb->GetSelIndex(selHnd, selAtoms, nSelAtoms);
    
    mmdb::Contact* contactsArray = NULL;
    int nContacts;
    float hbondDistance = floatParameters[std::string("hBondsStyleRMax")];
    mmdb->SeekContacts(selAtoms, nSelAtoms, selAtoms, nSelAtoms, 1.0, hbondDistance, 1, contactsArray, nContacts);
    
    int nDots = 0;
    std::shared_ptr<BallsPrimitive>balls(new BallsPrimitive());
    
    for (int iContact = 0; iContact < nContacts; iContact++){
        mmdb::Contact &contact = contactsArray[iContact];
        mmdb::Atom* atom1 = selAtoms[contact.id1];
        mmdb::Atom* atom2 = selAtoms[contact.id2];
        mmdb::Residue* residue1 = atom1->GetResidue();
        mmdb::Residue* residue2 = atom2->GetResidue();
        std::string atom1Name = std::string(atom1->name);
        std::string atom2Name = std::string(atom2->name);
#ifdef DEBUG_MINE
        std::cout << residue1->GetSeqNum() << " [" << atom1Name << "]" << residue2->GetSeqNum() << "[" << atom2Name << "[" << std::string(" N  ") << "]\n";
#endif
        if ((abs(residue1->GetSeqNum() - residue2->GetSeqNum()) > 2) ||
            !residue1->isAminoacid() ||
            !residue2->isAminoacid() || !(
                                          (atom1Name.compare(std::string(" N  "))==0 && atom2Name.compare(std::string(" N  "))==0) ||
                                          (atom1Name.compare(std::string(" O  "))==0 && atom2Name.compare(std::string(" O  "))==0) ||
                                          (atom1Name.compare(std::string(" N  "))==0 && atom2Name.compare(std::string(" O  "))==0) ||
                                          (atom1Name.compare(std::string(" O  "))==0 && atom2Name.compare(std::string(" N  "))==0) 
                                          )
            ){
            FCXXCoord atom1Coord(atom1->x, atom1->y, atom1->z);
            FCXXCoord atom2Coord(atom2->x, atom2->y, atom2->z);
            FCXXCoord diff = atom2Coord - atom1Coord;
            for (int iStep = 1; iStep < 7; iStep++){
                float step = (float)iStep / 8.;
                FCXXCoord ballCoord = atom1Coord + diff*step;
                FCXXCoord atomColor(colorScheme->colorForAtom(atom1, handles));
                if (iStep>5) atomColor = colorScheme->colorForAtom(atom2, handles);
                balls->addBall(ballCoord, atomColor, 0.1);
                nDots++;
                if (nDots%100 == 0) {
                    displayPrimitives.push_back(balls);
                    balls = std::shared_ptr<BallsPrimitive>(new BallsPrimitive());
                }
                
            }
        }
    }
    if (nDots%100 != 0) {
        displayPrimitives.push_back(balls);
    }

    mmdb->DeleteSelection(selHnd);
    colorScheme->freeSelectionHandles(mmdb, handles);

    return 0;
}

int MolecularRepresentation::drawBondsAsNewSticks()
{
    mmdb::Manager *mmdb = myMolecule->getMmdb();
	//selection->describe();
	selHnd = selection->handleInMMDB(mmdb);

    int nBonds = 0;
    std::shared_ptr<SticksPrimitive>sticks(new SticksPrimitive());
    sticks->setColorScheme(colorScheme);
    sticks->setMmdb(myMolecule->getMmdb());

    for (int iAtom = 1; iAtom <= myMolecule->getMmdb()->GetNumberOfAtoms(); iAtom++){
        mmdb::Atom* atom1 = myMolecule->getMmdb()->GetAtomI(iAtom);

        if (atom1->isInSelection(selHnd)){
            int nBondedAtoms;
            mmdb::AtomBond * bondedAtoms;
            atom1->GetBonds(bondedAtoms, nBondedAtoms);
            for (int iOtherAtom = 0; iOtherAtom < nBondedAtoms; iOtherAtom++){
                mmdb::Atom* atom2 = bondedAtoms[iOtherAtom].atom;
                if (atom2->GetIndex()>iAtom &&
                    atom2->isInSelection(selHnd)){

                    //Nasty kludge here...MMDB's MakeBonds screws up where there are multiple conformations
                    if (!strcmp(atom1->altLoc,"") &&
                        !strcmp(atom2->altLoc,"")) {
                        sticks->addPair(atom1, atom2);
                        nBonds++;
                    }
                    else {
                        float dx = atom1->x - atom2->x;
                        float dy = atom1->y - atom2->y;
                        float dz = atom1->z - atom2->z;
                        float dist = sqrtf (dx*dx + dy*dy + dz*dz);
                        if (dist < 1.9f) {
                            sticks->addPair(atom1, atom2);
                            nBonds++;
                        }
                    }

                    if (nBonds%1000 == 0) {
                        displayPrimitives.push_back(sticks);
                        sticks = std::shared_ptr<SticksPrimitive>(new SticksPrimitive());
                        sticks->setColorScheme(colorScheme);
                        sticks->setMmdb(myMolecule->getMmdb());
                    }
                }
            }
        }
    }
    if (nBonds%1000 != 0) {
        displayPrimitives.push_back(sticks);
    }
    mmdb->DeleteSelection(selHnd);
    return 0;
}

int MolecularRepresentation::drawDishyBases()
{
    std::map<mmdb::Chain *, DishyBaseContainer_t> dishy_bases_chain_map;
    mmdb::Manager *mmdb = myMolecule->getMmdb();
    //selection->describe();
    selHnd = selection->handleInMMDB(mmdb);
    int resultCode = myMolecule->identifyDishyBases(dishy_bases_chain_map, selHnd);

    FCXXCoord xAxis (1., 0., 0., 0.);
    FCXXCoord yAxis (0., 1., 0., 0.);
    FCXXCoord zAxis (0., 0., 1., 0.);

    std::shared_ptr<BallsPrimitive>balls(new BallsPrimitive());

    shared_ptr<CylindersPrimitive>cylinder(new CylindersPrimitive());
    cylinder->setAngularSampling(intParameters["dishStyleAngularSampling"]);

    float cylinderRadius = floatParameters[std::string("cylindersStyleCylinderRadius")];
    float ballRadius = floatParameters[std::string("cylindersStyleBallRadius")];

    auto dbContainerIter = dishy_bases_chain_map.begin();
    std::map<std::shared_ptr<ColorRule>,int>handles = colorScheme->prepareForMMDB(mmdb);

    for (; dbContainerIter!=dishy_bases_chain_map.end();dbContainerIter++){
        auto dbCont = dbContainerIter->second;
        auto dishyBaseIter = (dbContainerIter->second).bases.begin();
        for (; dishyBaseIter != (dbContainerIter->second).bases.end(); ++dishyBaseIter){

            auto atom1 = dishyBaseIter->ribose_atoms[1];
            FCXXCoord atom1Color =  colorScheme->colorForAtom(atom1, handles);

            balls->addBall(dishyBaseIter->centre, atom1Color, dishyBaseIter->radius, dishyBaseIter->normal, 0.14*dishyBaseIter->radius);

            if (balls->getBalls().size()%100 == 0){
                displayPrimitives.push_back(balls);
                balls = std::shared_ptr<BallsPrimitive>(new BallsPrimitive);
            }

            auto flatFanPrimitive = std::shared_ptr<FlatFanPrimitive>(new FlatFanPrimitive(dishyBaseIter->ribose_atoms, atom1Color));
            displayPrimitives.push_back(flatFanPrimitive);

            auto riboseAtomIter = dishyBaseIter->ribose_atoms.begin();
            for (; riboseAtomIter!= dishyBaseIter->ribose_atoms.end(); ++ riboseAtomIter){
                FCXXCoord coord((*riboseAtomIter)->x, (*riboseAtomIter)->y, (*riboseAtomIter)->z);
                FCXXCoord atomColor =  colorScheme->colorForAtom(*riboseAtomIter, handles);
                balls->addBall(coord, atomColor, ballRadius);
                if (balls->getBalls().size()%100 == 0){
                    displayPrimitives.push_back(balls);
                    balls = std::shared_ptr<BallsPrimitive>(new BallsPrimitive);
                }
            }

            auto bond = DishyBase_t::bondingPattern.begin();
            for (; bond!= DishyBase_t::bondingPattern.end(); ++bond){
                auto atom1 = dishyBaseIter->ribose_atoms[bond->first];
                FCXXCoord atom1Color =  colorScheme->colorForAtom(atom1, handles);
                auto atom2 = dishyBaseIter->ribose_atoms[bond->second];
                FCXXCoord atom2Color =  colorScheme->colorForAtom(atom2, handles);
                cylinder->addHalfAtomBond(atom1, atom1Color, atom2, atom2Color, cylinderRadius);
            }
            // Draw a stick from ribose_atoms[1] to 1/3 of the way to
            // centre.
            FCXXCoord atom1Coord(dishyBaseIter->ribose_atoms[1]->x,
                                 dishyBaseIter->ribose_atoms[1]->y,
                                 dishyBaseIter->ribose_atoms[1]->z);
            FCXXCoord basePseudoAtomPosition = atom1Coord + (dishyBaseIter->centre - atom1Coord) / 3.;
            cylinder->addHalfAtomBondWithCoords(atom1Coord, dishyBaseIter->ribose_atoms[1], atom1Color,
                                                basePseudoAtomPosition, dishyBaseIter->ribose_atoms[1], atom1Color,
                                                cylinderRadius);

        }
    }
    displayPrimitives.push_back(cylinder);
    if (balls->getBalls().size()%100 != 0){
        displayPrimitives.push_back(balls);
    }
    colorScheme->freeSelectionHandles(mmdb, handles);
    mmdb->DeleteSelection(selHnd);

    return 0;
}

int MolecularRepresentation::drawStickBases() {

   auto is_nucleic_acid = [] (const std::string &base_name) {
      if (base_name == "G") return true;
      if (base_name == "A") return true;
      if (base_name == "T") return true;
      if (base_name == "C") return true;
      if (base_name == "U") return true;
      if (base_name == "DG") return true;
      if (base_name == "DA") return true;
      if (base_name == "DC") return true;
      if (base_name == "DT") return true;
      // others
      return false;
   };

   auto get_atom_name_pair = [] (const std::string &base_name) {
      if (base_name == "DG" || base_name == "DA" || base_name == "G" || base_name == "A")
         return std::pair<std::string, std::string> (" C3'", " N1 ");
      if (base_name == "DT" || base_name == "U")
         return std::pair<std::string, std::string> (" C3'", " O4 ");
      if (base_name == "DC" || base_name == "C")
         return std::pair<std::string, std::string> (" C3'", " N4 ");
      return std::pair<std::string, std::string> ("", "");
   };

   mmdb::Manager *mmdb = myMolecule->getMmdb();
   selHnd = selection->handleInMMDB(mmdb);
   shared_ptr<CylindersPrimitive>cylinder(new CylindersPrimitive());
   cylinder->setAngularSampling(intParameters["cylindersStyleAngularSampling"]);
   float cylinderRadius = floatParameters[std::string("cylindersStyleCylinderRadius")];
   std::map<std::shared_ptr<ColorRule>,int>handles = colorScheme->prepareForMMDB(mmdb);
   float ballRadius = floatParameters[std::string("cylindersStyleBallRadius")];
   std::shared_ptr<BallsPrimitive>balls(new BallsPrimitive());

   {
      for(int imod = 1; imod<=mmdb->GetNumberOfModels(); imod++) {
         mmdb::Model *model_p = mmdb->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               int n_res = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     if (residue_p->isDNARNA()) {
                        std::string base_name(residue_p->GetResName());
                        if (is_nucleic_acid(base_name)) {

                           std::pair<std::string, std::string> atom_name_pair = get_atom_name_pair(base_name);
                           const std::string atom_name_1 = atom_name_pair.first;
                           const std::string atom_name_2 = atom_name_pair.second;
                           mmdb::Atom *atom_1 = nullptr;
                           mmdb::Atom *atom_2 = nullptr;

                           // iterate over all alt-confs in the residue
                           // (which I am not doing at the moment)

                           mmdb::PPAtom residue_atoms = nullptr;
                           int nResidueAtoms = 0;
                           residue_p->GetAtomTable(residue_atoms, nResidueAtoms);
                           for(int i = 0; i < nResidueAtoms; i++) {
                              mmdb::Atom *at = residue_atoms[i];
                              std::string atom_name(at->GetAtomName());
                              if (! atom_1)
                                 if (atom_name == atom_name_1)
                                    atom_1 = at;
                              if (! atom_2)
                                 if (atom_name == atom_name_2)
                                    atom_2 = at;
                           }
                           if (atom_1 && atom_2) {
                               if(atom_1->isInSelection(selHnd)&&atom_2->isInSelection(selHnd)) {
                                   FCXXCoord atom1Color = colorScheme->colorForAtom(atom_1, handles);
                                   FCXXCoord atom2Color = colorScheme->colorForAtom(atom_2, handles);
                                   cylinder->addHalfAtomBond(atom_1, atom1Color, atom_2, atom2Color, cylinderRadius);
                                   FCXXCoord atom1Coord(atom_2->x,atom_2->y, atom_2->z);
                                   balls->addBall(atom1Coord, atom1Color, ballRadius);
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

   mmdb->DeleteSelection(selHnd);
   displayPrimitives.push_back(cylinder);
   displayPrimitives.push_back(balls);
   colorScheme->freeSelectionHandles(mmdb, handles);
   return 0;
}

int MolecularRepresentation::drawRibbon()
{
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

    mmdb::Manager *mmdb = myMolecule->getMmdb();
    //selection->describe();
    selHnd = selection->handleInMMDB(mmdb);
    std::map<std::shared_ptr<ColorRule>,int>handles = colorScheme->prepareForMMDB(mmdb);

    std::vector<DiscreteSegment *> segments;
    myMolecule->identifySegments(segments, selHnd);

    for (std::size_t iSegment = 0; iSegment < segments.size(); iSegment++){
        DiscreteSegment &segment = *(segments[iSegment]);
        segment.evaluateNormals();
        if (boolParameters["smoothBetas"]) segment.smoothBetas();
        segment.evaluateSplines(subdivisionsPerCalpha);
        FCXXCoord color;

        std::shared_ptr<CylindersPrimitive>currentCylinder(new CylindersPrimitive());
        currentCylinder->setAngularSampling(intParameters["cylindersStyleAngularSampling"]);
        std::shared_ptr<BoxSectionPrimitive>currentBoxSection(new BoxSectionPrimitive());

        int lastSSE = -32767;
        int currentSSE = -32767;
        bool useArrowTipCoord = false;
        FCXXCoord arrowTipCoord;
        for (int iCalpha = 0; iCalpha<segment.nCalphas(); iCalpha++){
            const mmdb::Atom *calpha = (segment.calpha(iCalpha));
            auto residue = const_cast<mmdb::Atom*>(calpha)->GetResidue();
            if (residue->isDNARNA()) {
                currentSSE = 32767;
            }
            else currentSSE = const_cast<mmdb::Atom*>(calpha)->GetResidue()->SSE;
            //std::cout << const_cast<mmdb::Atom*>(calpha)->GetResidue()->GetResidueNo() <<":"<<currentSSE<<std::endl;
            if (currentSSE == mmdb::SSE_Bulge) currentSSE = mmdb::SSE_Strand;

            int nextSSE = mmdb::SSE_None;
            if (iCalpha < segment.nCalphas()-1){
                auto nextResidue = segment.calpha(iCalpha+1)->GetResidue();
                if (nextResidue->isDNARNA()) {
                    nextSSE = 32767;
                }
                else nextSSE = segment.calpha(iCalpha+1)->GetResidue()->SSE;
                if (nextSSE == mmdb::SSE_Bulge) nextSSE = mmdb::SSE_Strand;
            }
            if (lastSSE != currentSSE){
                if (lastSSE == mmdb::SSE_Strand) {
                    displayPrimitives.push_back(currentBoxSection);
                    currentBoxSection = std::shared_ptr<BoxSectionPrimitive>(new BoxSectionPrimitive());
                }
                else if (lastSSE != -32767){
                    displayPrimitives.push_back(currentCylinder);
                    currentCylinder = std::shared_ptr<CylindersPrimitive>(new CylindersPrimitive());
                    currentCylinder->setAngularSampling(intParameters["cylindersStyleAngularSampling"]);
                }
            }
            color = colorScheme->colorForAtom(calpha, handles);

            int endSubdivision = ((iCalpha == (segment.nCalphas()-1))?subdivisionsPerCalpha/2:subdivisionsPerCalpha);
            int startSubdivision = ((iCalpha == 0) ? (subdivisionsPerCalpha/2):0);
            float xVal = 0.;

            //Arrow head extrapolation state
            bool inArrowTaper = false;
            FCXXCoord arrowBaseCoord, arrowBaseTangent, arrowBaseNormalOne, arrowBaseNormalTwo;
            float arrowBaseXVal = 0.f;

            for (int i=startSubdivision; i<endSubdivision; i++){
                //for (int i=0; i<subdivisionsPerCalpha; i++){
                xVal = (iCalpha + (i*stepPerSubdivision)) - 0.5;
                FCXXCoord coord = segment.coordFor(xVal);
                FCXXCoord normalOne = segment.normalOneFor(xVal);
                FCXXCoord normalTwo = segment.normalTwoFor(xVal);
                //Re-orthonormalise against tangent after spline interpolation
                FCXXCoord tangentVec = segment.coordFor(xVal + stepPerSubdivision) - segment.coordFor(xVal - stepPerSubdivision);
                FCXXCoord tangent = tangentVec;
                tangent.normalise();
                normalOne = normalOne - tangent * (normalOne * tangent);
                normalOne.normalise();
                normalTwo = normalTwo - tangent * (normalTwo * tangent);
                normalTwo = normalTwo - normalOne * (normalTwo * normalOne);
                normalTwo.normalise();
                //Override with extrapolated frame during arrow taper
                if (inArrowTaper) {
                    coord = arrowBaseCoord + arrowBaseTangent * ((xVal - arrowBaseXVal) / (2.f * stepPerSubdivision));
                    normalOne = arrowBaseNormalOne;
                    normalTwo = arrowBaseNormalTwo;
                }
                //Start next SSE element at arrow tip for continuity
                if (useArrowTipCoord) {
                    coord = arrowTipCoord;
                    useArrowTipCoord = false;
                }
                //Widths at which things should be drawn is a matter of painful heuristics
                float radiusOne = radiusOneNone;
                float radiusTwo = radiusTwoNone;
                if (currentSSE == mmdb::SSE_Helix) {
                    if (i<subdivisionsPerCalpha/2){
                        if (currentSSE == lastSSE) {
                            radiusOne = radiusOneHelix;
                            radiusTwo = radiusTwoHelix;
                        }
                        else {
                            float factor = (float)i / ((float)subdivisionsPerCalpha/2.f);
                            radiusOne = radiusOneNone + factor * (radiusOneHelix-radiusOneNone);
                            radiusTwo = radiusTwoNone + factor * (radiusTwoHelix-radiusTwoNone);
                        }
                    }
                    else {
                        if (currentSSE == nextSSE) {
                            radiusOne = radiusOneHelix;
                            radiusTwo = radiusTwoHelix;
                        }
                        else {
                            float factor = (float)(i-(subdivisionsPerCalpha/2)) / ((float)subdivisionsPerCalpha/2.f);
                            radiusOne = radiusOneHelix - factor * (radiusOneHelix-radiusOneNone);
                            radiusTwo = radiusTwoHelix - factor * (radiusTwoHelix-radiusTwoNone);
                        }
                    }
                    CylinderPoint cylinderPoint(coord, color, normalOne, normalTwo, radiusOne, radiusTwo, calpha);
                    currentCylinder->addPoint(cylinderPoint);
                }
                else if (currentSSE == 32767) {
                    if (i<subdivisionsPerCalpha/2){
                        if (currentSSE == lastSSE) {
                            radiusOne = radiusOneDNARNA;
                            radiusTwo = radiusTwoDNARNA;
                        }
                        else {
                            float factor = (float)i / ((float)subdivisionsPerCalpha/2.f);
                            radiusOne = radiusOneNone + factor * (radiusOneDNARNA-radiusOneNone);
                            radiusTwo = radiusTwoNone + factor * (radiusTwoDNARNA-radiusTwoNone);
                        }
                    }
                    else {
                        if (currentSSE == nextSSE) {
                            radiusOne = radiusOneDNARNA;
                            radiusTwo = radiusTwoDNARNA;
                        }
                        else {
                            float factor = (float)(i-(subdivisionsPerCalpha/2)) / ((float)subdivisionsPerCalpha/2.f);
                            radiusOne = radiusOneDNARNA - factor * (radiusOneDNARNA-radiusOneNone);
                            radiusTwo = radiusTwoDNARNA - factor * (radiusTwoDNARNA-radiusTwoNone);
                        }
                    }
                    CylinderPoint cylinderPoint(coord, color, normalOne, normalTwo, radiusOne, radiusTwo, calpha);
                    currentCylinder->addPoint(cylinderPoint);
                }
                else if (currentSSE == mmdb::SSE_Strand) {
                    if (i==startSubdivision && lastSSE != currentSSE) {
                        CylinderPoint cylinderPoint(coord, color, normalOne, normalTwo, radiusOneNone, radiusTwoNone, calpha);
                        currentBoxSection->addPoint(cylinderPoint);
                        lastSSE = currentSSE;
                    }
                    if (i<subdivisionsPerCalpha/2){
                        if (currentSSE == lastSSE) {
                            radiusOne = radiusOneStrand;
                            radiusTwo = radiusTwoStrand;
                        }
                    }
                    if (i == subdivisionsPerCalpha/2){
                        if (currentSSE != nextSSE){
                            //Capture frame at arrow base for extrapolation
                            arrowBaseCoord = coord;
                            arrowBaseTangent = tangentVec;
                            arrowBaseNormalOne = normalOne;
                            arrowBaseNormalTwo = normalTwo;
                            arrowBaseXVal = xVal;
                            inArrowTaper = true;
                            CylinderPoint cylinderPoint(coord, color, normalOne, normalTwo, radiusOneStrand, radiusTwoStrand, calpha);
                            currentBoxSection->addPoint(cylinderPoint);
                        }
                    }
                    if (i>=subdivisionsPerCalpha/2){
                        if (currentSSE == nextSSE) {
                            radiusOne = radiusOneStrand;
                            radiusTwo = radiusTwoStrand;
                        }
                        else {
                            float factor = (float)(i-(subdivisionsPerCalpha/2)) / ((float)subdivisionsPerCalpha/2.f);
                            radiusOne = radiusOneArrow - factor * (radiusOneArrow-radiusOneNone);
                            radiusTwo = radiusTwoArrow - factor * (radiusTwoArrow-radiusTwoNone);
                        }
                    }
                    CylinderPoint cylinderPoint(coord, color, normalOne, normalTwo, radiusOne, radiusTwo, calpha);
                    currentBoxSection->addPoint(cylinderPoint);
                }
                else {
                    CylinderPoint cylinderPoint(coord, color, normalOne, normalTwo, radiusOne, radiusTwo, calpha);
                    currentCylinder->addPoint(cylinderPoint);
                }
            }
            lastSSE = currentSSE;

            //Add an extra segment to delimit the end of a residue
            xVal += stepPerSubdivision;
            FCXXCoord coord = segment.coordFor(xVal);
            FCXXCoord normalOne = segment.normalOneFor(xVal);
            FCXXCoord normalTwo = segment.normalTwoFor(xVal);
            //Re-orthonormalise against tangent after spline interpolation
            FCXXCoord tangent = segment.coordFor(xVal + stepPerSubdivision) - segment.coordFor(xVal - stepPerSubdivision);
            tangent.normalise();
            normalOne = normalOne - tangent * (normalOne * tangent);
            normalOne.normalise();
            normalTwo = normalTwo - tangent * (normalTwo * tangent);
            normalTwo = normalTwo - normalOne * (normalTwo * normalOne);
            normalTwo.normalise();
            //Extrapolate arrow tip and store for next SSE continuity
            if (inArrowTaper) {
                coord = arrowBaseCoord + arrowBaseTangent * ((xVal - arrowBaseXVal) / (2.f * stepPerSubdivision));
                normalOne = arrowBaseNormalOne;
                normalTwo = arrowBaseNormalTwo;
                arrowTipCoord = coord;
                useArrowTipCoord = true;
                inArrowTaper = false;
            }

            float radiusOne;
            float radiusTwo;
            radiusOne = radiusOneNone;
            radiusTwo = radiusTwoNone;
            if (nextSSE == mmdb::SSE_Helix && lastSSE == mmdb::SSE_Helix) {
                radiusOne = radiusOneHelix;
                radiusTwo = radiusTwoHelix;
            }
            else if (nextSSE == 32767 && lastSSE == 32767){
                radiusOne = radiusOneDNARNA;
                radiusTwo = radiusTwoDNARNA;
            }
            else if (nextSSE == mmdb::SSE_Strand && lastSSE == mmdb::SSE_Strand){
                radiusOne = radiusOneStrand;
                radiusTwo = radiusTwoStrand;
            }
            CylinderPoint cylinderPoint(coord, color, normalOne, normalTwo, radiusOne, radiusTwo, calpha);
            if (currentSSE == mmdb::SSE_Strand) currentBoxSection->addPoint(cylinderPoint);
            else currentCylinder->addPoint(cylinderPoint);
        }
        if (currentSSE == mmdb::SSE_Strand) displayPrimitives.push_back(currentBoxSection);
        else if (currentSSE != -32767) displayPrimitives.push_back(currentCylinder);
    }
    colorScheme->freeSelectionHandles(mmdb, handles);
    mmdb->DeleteSelection(selHnd);

    return 0;
}

int MolecularRepresentation::drawCalphas()
{
    mmdb::Manager *mmdb = myMolecule->getMmdb();
	//selection->describe();
	selHnd = selection->handleInMMDB(mmdb);
    std::map<std::shared_ptr<ColorRule>,int>handles = colorScheme->prepareForMMDB(mmdb);

    std::vector<DiscreteSegment *> segments;
    myMolecule->identifySegments(segments, selHnd);

    std::shared_ptr<BondsPrimitive>bonds(new BondsPrimitive());
    bonds->setColorScheme(colorScheme);
    displayPrimitives.push_back(bonds);

    for (std::size_t iSegment = 0; iSegment < segments.size(); iSegment++){
        DiscreteSegment &segment = *(segments[iSegment]);

        for (int i=0; i<(segment.nCalphas()-1); i++){
            mmdb::Atom* atom1 = segment.calpha(i);
            mmdb::Atom* atom2 = segment.calpha(i+1);
            bonds->addPair(atom1, atom2);
            bonds->addPair(atom2, atom1);
        }
    }
    bonds->evaluateGLPrimitives(handles);

    colorScheme->freeSelectionHandles(mmdb, handles);
    mmdb->DeleteSelection(selHnd);

    return 0;
}

int MolecularRepresentation::drawMolecularSurface()
{
    return drawSurfaceOfKind(SurfacePrimitive::MolecularSurface);
}

int MolecularRepresentation::drawSurfaceOfKind(int surfaceKind)
{
    //return 0;

    mmdb::Manager *mmdb = myMolecule->getMmdb();
	//selection->describe();
	selHnd = selection->handleInMMDB(mmdb);
	std::map<std::shared_ptr<ColorRule>,int>handles = colorScheme->prepareForMMDB(mmdb);

    //If we are limited to short int for indexes, we wil have to chop this surface
    //into bite-sized chunks
    mmdb::Atom** SelAtoms;
    int nAtoms;

    int chunkSize = 1000000;
    if (sizeof(GLIndexType) == sizeof(short)) chunkSize = 100;
    mmdb->GetSelIndex(selHnd, SelAtoms, nAtoms);

    int iAtom = 0;
    int chunkHndl = mmdb->NewSelection();

    float probeRadius = floatParameters["surfaceStyleProbeRadius"];
    float radiusMultiplier = floatParameters["ballsStyleRadiusMultiplier"];

    for (; iAtom<nAtoms; iAtom++){
        mmdb->SelectAtom(chunkHndl, SelAtoms[iAtom], mmdb::SKEY_OR);
        if ((iAtom+1)%chunkSize == 0){ 
            std::shared_ptr<SurfacePrimitive>surfacePrimitive(new SurfacePrimitive(mmdb, chunkHndl, selHnd, colorScheme, (SurfacePrimitive::SurfaceType) surfaceKind, probeRadius, radiusMultiplier));
            //surfacePrimitive->generateArrays();
            if (surfacePrimitive->getCXXSurfaceMaker()) {
                displayPrimitives.push_back(surfacePrimitive);
            }
            mmdb->DeleteSelection(chunkHndl);
            chunkHndl = mmdb->NewSelection();
        }
        if (redrawProgressCallback){
            redrawProgressCallback(redrawProgressCallbackUserInfo, (float)iAtom / (float)nAtoms);
        }

    }
    if (iAtom%chunkSize != 0){
        std::shared_ptr<SurfacePrimitive>surfacePrimitive(new SurfacePrimitive(mmdb, chunkHndl, selHnd, colorScheme, (SurfacePrimitive::SurfaceType) surfaceKind, probeRadius, radiusMultiplier));
        if (surfacePrimitive->getCXXSurfaceMaker()) {
            displayPrimitives.push_back(surfacePrimitive);
        }
        mmdb->DeleteSelection(chunkHndl);
    }

    colorScheme->freeSelectionHandles(mmdb, handles);
    mmdb->DeleteSelection(selHnd);

    return 0;
}


int MolecularRepresentation::drawVdWSurface()
{
    return drawSurfaceOfKind(SurfacePrimitive::VdWSurface);
}

int MolecularRepresentation::drawAccessibleSurface()
{
    return drawSurfaceOfKind(SurfacePrimitive::AccessibleSurface);
}

#include "MoleculesToTriangles/CXXSurface/CXXChargeTable.h"
#include "MoleculesToTriangles/CXXSurface/CXXCreator.h"
#include "AtomPropertyRampColorRule.h"

void MolecularRepresentation::colorByOwnPotential()
{
    colorByPotential(getCompoundSelection()->getSelectionString(), getMolecule());
}

void MolecularRepresentation::colorByPotential(std::string chargingAtomString, std::shared_ptr<MyMolecule> theMolecule)
{
    CompoundSelection chargingAtoms(chargingAtomString,"ChargingAtoms");
    int selHnd = chargingAtoms.handleInMMDB(theMolecule->getMmdb());
    //Instantiate an electrostatics map and cause it to calculate itself
    CXXChargeTable theChargeTable;
    CXXUtils::assignCharge(theMolecule->getMmdb(), selHnd, &theChargeTable);
    CXXCreator *theCreator = new CXXCreator(theMolecule->getMmdb(), selHnd);
    theCreator->calculate();
    clipper::Cell cell;
    auto theClipperNXMap = theCreator->coerceToClipperMap(cell);
    // Now bring the surface and the map together
    double coords[4];

    //Create a color ramp rule
    AtomPropertyRampColorRule rampRule;
    rampRule.setStartRGB(FCXXCoord (1.,0.,0.,1.));
    rampRule.setStartValue(-0.5);
    rampRule.setMiddleRGB(FCXXCoord (1.,1.,1.,1.));
    rampRule.setEndRGB(FCXXCoord (0.,0.,1.,1.));
    rampRule.setEndValue( 0.5);

    auto primitivePntr = getDisplayPrimitives().begin();
    for (; primitivePntr != getDisplayPrimitives().end(); ++primitivePntr){
        if (VertexColorNormalPrimitive *p = dynamic_cast<VertexColorNormalPrimitive *>(primitivePntr->get())){
            double pMin = 1e30;
            double pMax = -1e30;
            FCXXCoord  colMin = FCXXCoord (0.,0.,0.,0.);
            FCXXCoord  colMax = FCXXCoord (0.,0.,0.,0.);
            for (size_t iVertex =0; iVertex < p->nVertices(); iVertex++){
                auto vertex = p->getVertexColorNormalArray()[iVertex];
                clipper::Coord_orth orthogonals(vertex.vertex[0], vertex.vertex[1], vertex.vertex[2]);
                double potential;
                const clipper::Coord_map mapUnits(theClipperNXMap.coord_map(orthogonals));
                potential = theClipperNXMap.interp<clipper::Interp_cubic>( mapUnits );
                FCXXCoord  correspondingColor = rampRule.colorForValue(potential);
                if (potential<pMin) {
                    pMin = potential;
                    colMin = correspondingColor;
                }
                if (potential>pMax) {
                    pMax = potential;
                    colMax = correspondingColor;
                }
                for (int i=0; i<3; i++) p->getVertexColorNormalArray()[iVertex].color[i] = correspondingColor[i];
                p->getVertexColorNormalArray()[iVertex].color[3] = 1.;
            }
            std::cout << "pMin was " << pMin << " pMax was " << pMax << colMin << colMax;
            //Any renderers will have to update the buffers in which color information is stored
            //p->liberateAllHandles();
        }
    }
}


