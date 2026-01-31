/*
 * MoleculesToTriangles/CXXClasses/MyMolecule.cpp
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

#include <string>
#include <cstring>
#include "MyMolecule.h"
#include "DiscreteSegment.h"
//#include "DisplayPrimitive.h"
#include "MoleculesToTriangles/CXXSurface/CXXUtils.h"
#include "MoleculesToTriangles/CXXSurface/TokenIterator.h"
//#include "ColorScheme.h"
//#include "ObjectSelection.h"
#include "mmdb2/mmdb_manager.h"
#include "CompoundSelection.h"
#include "DishyBase.h"

//#include "mmdb_sbase.h"
#include <set>

MyMolecule::MyMolecule() : doDraw(true), ownsMMDB(true) {
    mmdb = new mmdb::Manager();
};
MyMolecule::~MyMolecule(){
    if (ownsMMDB) delete mmdb;
};

int MyMolecule::loadCoords( char *data, int length, int secondaryStructureUsageFlag)
{
    Delimiters delimiters("\n\r");
    TokenIterator<char*, Delimiters> charIter(data, data + length, delimiters), end2;
    std::string token;
    mmdb->SetFlag( mmdb::MMDBF_AutoSerials | mmdb::MMDBF_IgnoreDuplSeqNum );
    char tokenBuffer[256];
    bool isDamminModel = false;
    int iDamminAtom = 1;
    while (charIter != end2){
        token = *charIter++;
        strcpy (tokenBuffer, token.c_str());
        //Here a funky kludge to recognise and handle dummy atom molecules from Dammin
        if (!strncmp(tokenBuffer, " Project description:", 21)){
            isDamminModel = true;
        }
        if (!isDamminModel) {
            int RC = mmdb->PutPDBString ( tokenBuffer );
        }
        else {
            char formattedCard[256 + 7];
            if (!strncmp ("ATOM", tokenBuffer, 4)){
                AtomCard atomCard;
                
                atomCard.selected = 0;
                
                strncpy (atomCard.atname, " DAM", 4);
                atomCard.atname[4] = '\0';
                
                strncpy (atomCard.alt, tokenBuffer+16, 1);
                atomCard.alt[1] = '\0';
                
                strncpy (atomCard.restype, "DAM ", 4);
                atomCard.restype[4] = '\0';
                
                strncpy (atomCard.chainid, tokenBuffer+21, 1);
                atomCard.chainid[1] = '\0';
                
                strncpy (atomCard.resname, tokenBuffer+22, 4);
                atomCard.resname[4] = '\0';
                char field[12];
                
                field[8] = '\0';
                
                strncpy (field, tokenBuffer+30, 8);
                atomCard.crd.xyz[0] = (float) atof(field);
                strncpy (field, tokenBuffer+38, 8);
                atomCard.crd.xyz[1] = atof(field);
                strncpy (field, tokenBuffer+46, 8);
                atomCard.crd.xyz[2] = atof(field);
                field[6] = '\0';
                strncpy (field, tokenBuffer+54, 6);
                atomCard.properties[0] = atof(field);
                strncpy (field, tokenBuffer+60, 6);
                atomCard.properties[1] = atof(field);
                
                FormatPDBCard(atomCard, formattedCard, iDamminAtom++);
            }
            else {
               snprintf(formattedCard, 256 + 7, "REMARK %s", tokenBuffer);
            }
            int RC = mmdb->PutPDBString ( formattedCard );
        }
    }
    return processCoords(secondaryStructureUsageFlag);
}

void
secondary_structure_header_to_residue_sse(mmdb::Manager *mol) {

   // add the SSE atribute for residue of the mol, extracted
   // from the HELIX and SHEET records

   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_sheets  = model_p->GetNumberOfSheets();
         int n_helices = model_p->GetNumberOfHelices();

         for (int ih=1; ih<=n_helices; ih++) {
            mmdb:: Helix *helix_p = model_p->GetHelix(ih);
            if (helix_p) {
               std::string chain_id = helix_p->initChainID;

               int n_chains = model_p->GetNumberOfChains();
               for (int ichain=0; ichain<n_chains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  if (chain_p) {
                     std::string this_chain_id = chain_p->GetChainID();
                     if (this_chain_id == chain_id) {
                        int start_res_no = helix_p->initSeqNum;
                        int   end_res_no = helix_p->endSeqNum;
                        int n_res = chain_p->GetNumberOfResidues();
                        for (int ires=0; ires<n_res; ires++) {
                           mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                           if (residue_p) {
                              if (residue_p->GetSeqNum() >= start_res_no) {
                                 if (residue_p->GetSeqNum() <= end_res_no) {
                                    residue_p->SSE = mmdb::SSE_Helix;
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }

         for (int isheet=1; isheet<=n_sheets; isheet++) {
            mmdb::Sheet *sheet_p = model_p->GetSheet(isheet);
            if (sheet_p) {
               int n_strands = model_p->GetNumberOfStrands(isheet);
               for (int istrand=1; istrand<=n_strands; istrand++) {
                  mmdb::Strand *strand_p = model_p->GetStrand(isheet, istrand);
                  if (strand_p) {
                     std::string chain_id = strand_p->initChainID;
                     int start_res_no = strand_p->initSeqNum;
                     int   end_res_no = strand_p->endSeqNum;
                     int n_chains = model_p->GetNumberOfChains();
                     for (int ichain=0; ichain<n_chains; ichain++) {
                        mmdb::Chain *chain_p = model_p->GetChain(ichain);
                        if (chain_p) {
                           std::string this_chain_id = chain_p->GetChainID();
                           if (this_chain_id == chain_id) {
                              int n_res = chain_p->GetNumberOfResidues();
                              for (int ires=0; ires<n_res; ires++) {
                                 mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                                 if (residue_p) {
                                    if (residue_p->GetSeqNum() >= start_res_no) {
                                       if (residue_p->GetSeqNum() <= end_res_no) {
                                          residue_p->SSE = mmdb::SSE_Strand;
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
   }
}


int MyMolecule::processCoords(int secondaryStructureUsageFlag){

   bool calc_secondary_structure = false;
   bool transfer_secondary_structure = false;
   if (secondaryStructureUsageFlag == CALC_SECONDARY_STRUCTURE)
      calc_secondary_structure = true;
   if (secondaryStructureUsageFlag == USE_HEADER_INFO)
      transfer_secondary_structure = true;

   //mmdb->Cryst.CalcCoordTransforms();
   identifyBonds();
   CXXUtils::assignUnitedAtomRadius(mmdb);

   if (calc_secondary_structure) {
      int nModels = mmdb->GetNumberOfModels();
      for (int iModel = 1; iModel <= nModels; iModel++){
         mmdb::Model *model = mmdb->GetModel(iModel);
        model->CalcSecStructure(true);
      }
   }

   if (transfer_secondary_structure)
      secondary_structure_header_to_residue_sse(mmdb);

   return 0;
}

int MyMolecule::loadFromPDB(const char *filePath, int secondaryStructureUsageFlag){
    int RC;
    mmdb::InitMatType();
    mmdb = new mmdb::Manager();

    //Now read the MMDB for purpose of calculating surface
    mmdb->SetFlag( mmdb::MMDBF_PrintCIFWarnings );
    RC = mmdb->ReadCoorFile (filePath);
    // if RC > 0 reading in file has failed - say so and quit
    if (RC) {
        std::cout << "error could not read file " <<
        filePath << endl;
    }
    else {
        std::cout << processCoords(secondaryStructureUsageFlag);
    }
    return RC;
}

MyMolecule::MyMolecule(const char *filePath, int secondaryStructureUsageFlag) : MyMolecule()
{
    int RC = loadFromPDB(filePath, secondaryStructureUsageFlag);
}

MyMolecule::MyMolecule(std::string filePathString, int secondaryStructureUsageFlag) : MyMolecule()
{
    int RC = loadFromPDB(filePathString.c_str(), secondaryStructureUsageFlag);
}

int MyMolecule::FormatPDBCard (AtomCard theAtom, char *card,int count){
    char token[200];
    
    strcpy (card,"ATOM  ");
    snprintf(token, 198, "%5d ", count + 1);
    strcat (card, token);
    strcat (card, theAtom.atname);
    strcat (card, theAtom.alt);
    strcat (card, theAtom.restype);
    strcat (card, theAtom.chainid);
    strcat (card, theAtom.resname);
    snprintf(token, 198, "    %8.3f%8.3f%8.3f%6.2f%6.2f               ",
            theAtom.crd.xyz[0], theAtom.crd.xyz[1], theAtom.crd.xyz[2],
            theAtom.properties[0], theAtom.properties[1]);
    strcat (card, token);
    return (1);
}


int MyMolecule::identifySegments(std::vector<DiscreteSegment *> &segments, int selHnd)
{
    
    bool do_dishy_bases = true; // I don't know where this goes - it's a user setting
    // this is where they go if they were calculated
    std::map<mmdb::Chain *, DishyBaseContainer_t> dishy_bases_chain_map;
    
    int nModels = mmdb->GetNumberOfModels();
    // std::cout << "debug:: identifySegments(): Protein consists of " << nModels << " Models\n";
    mmdb::Model* model = mmdb->GetModel(1);
    
    //Here is an approach for (bloomin) segIDs....repeat the segment determination for each segID
    // simples :-)
    
    std::set<std::string> segIDs;
    
    {
        // Establish a list of segIDs
        int nAtoms = model->GetNumberOfAtoms(false);
        mmdb::Atom** atoms = model->GetAllAtoms();
        
        for (int iAtom = 0; iAtom < nAtoms; iAtom++){
            segIDs.insert(std::string(atoms[iAtom]->segID));
        }
    }
    // std::cout << "debug:: identifySegments(): there are " << segIDs.size() << " segIDs\n";

    std::set<std::string>::iterator segIDIter = segIDs.begin();
    for (;segIDIter != segIDs.end(); ++segIDIter){
        int          nChains;
        mmdb::Chain**     chains;
        model->GetChainTable(chains, nChains);
        // std::cout << "debug:: identifySegments(): Model 1 consists of " << nChains << " Chains\n";
        //Pass through and identify any chain breaks
        for (int iChain = 0; iChain < nChains; iChain ++){
           mmdb::Chain* chain = model->GetChain(iChain);
           FCXXCoord lastCoord (-1e30, -1e30, -1e30);
           int nResidues;
           mmdb::Residue**   residues;
           chain->GetResidueTable(residues, nResidues);
           // std::cout << "debug:: identifySegments(): Chain " << chain->GetChainID() << " consists of " << nResidues << " Residues\n";
            if (nResidues > 0){
                for (int iResidue = 0; iResidue < nResidues; iResidue++){
                    if (residues[iResidue]->isAminoacid() ||
                        !strcmp(residues[iResidue]->GetResName(),"MSE") ||
                        !strcmp(residues[iResidue]->GetResName(),"TPO") ||
                        !strcmp(residues[iResidue]->GetResName(),"THP") ||
                        !strcmp(residues[iResidue]->GetResName(),"SEP") ||
                        !strcmp(residues[iResidue]->GetResName(),"PTR")
                        ){
                        mmdb::Residue* residue = chain->GetResidue(iResidue);
                        mmdb::Atom** atomsOfResidue;
                        int nAtoms;
                        residue->GetAtomTable(atomsOfResidue, nAtoms);
                        for (int iAtom=0; iAtom < nAtoms; iAtom++){
                            mmdb::Atom* calpha = atomsOfResidue[iAtom];
                            if (std::string(calpha->segID) == *segIDIter){
                                //std::cout << calpha->segID << "oops\n";
                                if (!strcmp(calpha->name," CA ") &&
                                    calpha->isInSelection(selHnd)){
                                    //Consider only the main alternative location
                                    if (!strcmp(calpha->altLoc,"") ||
                                        !strcmp(calpha->altLoc,"A") ||
                                        calpha->occupancy > 0.5){
                                        FCXXCoord calphaPosition(calpha->x, calpha->y, calpha->z);
                                        FCXXCoord difference = calphaPosition - lastCoord;
                                        float distance = difference.get3DLength();
                                        if (distance > 4.1){
                                            // std::cout << "New chain start "<< calpha->GetResName() << calpha->GetSeqNum () << calpha->segID << "\n";
                                            DiscreteSegment *segment = new DiscreteSegment();
                                            auto [ax, ay, az] = getRadiiForResidue(calpha->GetResidue());
                                            segment->addCalpha(calpha, ax, ay, az);
                                            segments.push_back(segment);
                                        }
                                        else {
                                            auto [ax, ay, az] = getRadiiForResidue(calpha->GetResidue());
                                            segments.back()->addCalpha(calpha, ax, ay, az);
                                        }
                                        lastCoord = calphaPosition;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            // shall we do nucleic acid?
            if (chain->isNucleotideChain()) {
                int nResidues;
                mmdb::Residue**   residues;
                chain->GetResidueTable(residues, nResidues);
                if (nResidues > 0){
                    for (int iResidue = 0; iResidue < nResidues; iResidue++){
                        mmdb::Residue *residue_p = residues[iResidue];
                        if (residue_p) {
                            if (residue_p->isNucleotide()) {
                                mmdb::Atom **residue_atoms = 0;
                                int nAtoms;
                                residue_p->GetAtomTable(residue_atoms, nAtoms);
                                for (int iAtom=0; iAtom < nAtoms; iAtom++){
                                    mmdb::Atom* atom_p = residue_atoms[iAtom];
                                    std::string atom_name(atom_p->name);
                                    // if (atom_name == " P  ") {
                                    if (atom_name == " C3'") {
                                        if (atom_p->isInSelection(selHnd)) {
                                            FCXXCoord atom_pos(atom_p->x, atom_p->y, atom_p->z);
                                            FCXXCoord difference = atom_pos - lastCoord;
                                            float distance = difference.get3DLength();
                                            // std::cout << "distance " << distance << std::endl; almost all less than 7.5A
                                            if (distance > 8.0){ // heuristic
                                                // std::cout << "New chain nucleotide start "<< atom_p->GetResName() << " "
                                                // << atom_p->GetSeqNum () << atom_p->GetChainID() << "\n";
                                                DiscreteSegment *segment = new DiscreteSegment();
                                                auto [ax, ay, az] = getRadiiForResidue(residue_p);
                                                segment->addCalpha(atom_p, ax, ay, az);
                                                segments.push_back(segment);
                                            } else {
                                                auto [ax, ay, az] = getRadiiForResidue(residue_p);
                                                segments.back()->addCalpha(atom_p, ax, ay, az);
                                            }
                                            lastCoord = atom_pos;
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
    return 0;
}

int MyMolecule::identifyDishyBases(std::map<mmdb::Chain *, DishyBaseContainer_t> &dishy_bases_chain_map, int selHnd)
{
    int nModels = mmdb->GetNumberOfModels();
    // std::cout << "debug:: identifyDishyBases(): Protein consists of " << nModels << " Models\n";
    mmdb::Model* model = mmdb->GetModel(1);
    
    //Here is an approach for (bloomin) segIDs....repeat the segment determination for each segID
    // simples :-)
    
    std::set<std::string> segIDs;
    
    {
        // Establish a list of segIDs
        int nAtoms = model->GetNumberOfAtoms(false);
        mmdb::Atom** atoms = model->GetAllAtoms();
        
        for (int iAtom = 0; iAtom < nAtoms; iAtom++){
            segIDs.insert(std::string(atoms[iAtom]->segID));
        }
    }
    
    std::set<std::string>::iterator segIDIter = segIDs.begin();
    for (;segIDIter != segIDs.end(); ++segIDIter){
        int          nChains;
        mmdb::Chain**     chains;
        model->GetChainTable(chains, nChains);
        // std::cout << "debug:: identifyDishyBases(): Model 1 consists of " << nChains << " Chains\n";
        //Pass through and identify any chain breaks
        for (int iChain = 0; iChain < nChains; iChain ++){
            mmdb::Chain* chain = model->GetChain(iChain);
            FCXXCoord lastCoord (-1e30, -1e30, -1e30);
            int nResidues;
            mmdb::Residue**   residues;
            chain->GetResidueTable(residues, nResidues);
            // std::cout << "debug:: identifyDishyBases(): Chain " << chain->GetChainID() << " consists of " << nResidues << " Residues\n";
            
            // shall we do nucleic acid?
            if (chain->isNucleotideChain()) {
                int nResidues = 0;
                mmdb::Residue**   residues = nullptr;
                chain->GetResidueTable(residues, nResidues);
                if (nResidues > 0){
                    DishyBaseContainer_t dishybases;
                    int res_no_min =  9999;
                    int res_no_max = -9999;
                    int n_residues = 0;
                    
                    for (int iResidue = 0; iResidue < nResidues; iResidue++){
                        mmdb::Residue *residue_p = residues[iResidue];
                        if (residue_p) {
                            n_residues++;
                            mmdb::Atom **residue_atoms = 0;
                            int nAtoms;
                            residue_p->GetAtomTable(residue_atoms, nAtoms);
                            std::string res_name = residue_p->GetResName();
                            
                            if(nAtoms>0){
                               if(!residue_atoms[0]->isInSelection(selHnd)) {
                                   continue;
                               }
                            }
                            // setup the alt confs for this residue
                            //
                            std::set<std::string> residue_alt_confs_set;
                            for (int iAtom=0; iAtom < nAtoms; iAtom++){
                                mmdb::Atom* atom_p = residue_atoms[iAtom];
                                if (! atom_p->isTer()) {
                                    std::string atom_alt_conf(atom_p->altLoc);
                                    residue_alt_confs_set.insert(atom_alt_conf);
                                }
                            }
                            
                            // for each alt conf
                            //
                            for (std::size_t ialt=0; ialt<residue_alt_confs_set.size(); ialt++) {
                                
                                std::vector<std::string> ref_base_names;
                                std::vector<mmdb::Atom *> base_atoms;
                                if (res_name == "C" ) ref_base_names = dishybases.cytidine_base_names;
                                if (res_name == "DC") ref_base_names = dishybases.cytidine_base_names;
                                if (res_name == "T" ) ref_base_names = dishybases.thymine_base_names;
                                if (res_name == "DT") ref_base_names = dishybases.thymine_base_names;
                                if (res_name == "U" ) ref_base_names = dishybases.uracil_base_names;
                                if (res_name == "DU") ref_base_names = dishybases.uracil_base_names;
                                if (res_name == "A" ) ref_base_names = dishybases.adenine_base_names;
                                if (res_name == "DA") ref_base_names = dishybases.adenine_base_names;
                                if (res_name == "G" ) ref_base_names = dishybases.guanine_base_names;
                                if (res_name == "DG") ref_base_names = dishybases.guanine_base_names;
                                float radius = 2.7; // guess - as is 4.3 below
                                if (res_name == "A" || res_name == "DA" || res_name == "G" || res_name == "DG")
                                    radius = 3.4;
                                
                                std::vector<mmdb::Atom *> ribose_atoms(5,0);
                                for (int iAtom=0; iAtom < nAtoms; iAtom++){
                                    mmdb::Atom* atom_p = residue_atoms[iAtom];
                                    std::string atom_name(atom_p->name);
                                    std::string atom_alt_conf(atom_p->altLoc);
                                    if (atom_alt_conf.empty() ||
                                        (residue_alt_confs_set.find(atom_alt_conf) != residue_alt_confs_set.end())) {
                                        if (std::find(ref_base_names.begin(), ref_base_names.end(), atom_name) != ref_base_names.end()) {
                                            base_atoms.push_back(atom_p);
                                        }
                                        // was it a ribose atom?
                                        if (atom_name == " O4'") ribose_atoms[0] = atom_p;
                                        if (atom_name == " C1'") ribose_atoms[1] = atom_p;
                                        if (atom_name == " C2'") ribose_atoms[2] = atom_p;
                                        if (atom_name == " C3'") ribose_atoms[3] = atom_p;
                                        if (atom_name == " C4'") ribose_atoms[4] = atom_p;
                                    }
                                }
                                // we have all the ribose atoms?
                                for (std::size_t i=0; i<5; i++)
                                    if (!ribose_atoms[i])
                                        continue;
                                // we have some base atoms?
                                if (base_atoms.size() < 4)
                                    continue;
                                
                                // Moving on - standard, non-pathalogical case
                                //
                                FCXXCoord ribose_centre;
                                for (std::size_t i=0; i<5; i++) {
                                    FCXXCoord pos(ribose_atoms[i]->x, ribose_atoms[i]->y, ribose_atoms[i]->z);
                                    ribose_centre += pos;
                                }
                                ribose_centre *= 0.2;
                                FCXXCoord base_centre;
                                for (std::size_t i=0; i<base_atoms.size(); i++) {
                                    FCXXCoord pos(base_atoms[i]->x, base_atoms[i]->y, base_atoms[i]->z);
                                    base_centre += pos;
                                }
                                base_centre /= float(base_atoms.size());
                                std::vector<FCXXCoord> base_atom_positions(base_atoms.size());
                                for (unsigned int i=0; i<base_atoms.size(); i++)
                                    base_atom_positions[i] = FCXXCoord(base_atoms[i]->x, base_atoms[i]->y, base_atoms[i]->z);
                                DishyPlaneLSQ_t lsq(base_atom_positions);
                                FCXXCoord base_normal = lsq.normal();
                                DishyBase_t db(base_centre, base_normal, radius, ribose_atoms, ribose_centre);
                                dishybases.add(db);
                            }
                        }
                    }
                    dishy_bases_chain_map[chain] = dishybases;
                }
            }
        }
    }
    return 0;
}

int MyMolecule::identifyBonds()
{
    //PCSBase      sbase = new CSBase();
    /*
     int RC = sbase->LoadIndex (path);
     if (RC!=SBASE_Ok)  {
     printf ( " ****** compbase files not found.\n" );
     return result;
     }
     PCFile structFile = sbase->GetStructFile();
     */
    mmdb->MakeBonds(true);
    int nModels = mmdb->GetNumberOfModels();
    // std::cout << "debug:: identifyBonds(): Protein consists of " << nModels << " Models\n";
    mmdb::Model* model = mmdb->GetModel(1);
    
    int          nChains;
    mmdb::Chain**     chains;
    model->GetChainTable(chains, nChains);
    // std::cout << "debug:: identifyBonds(): Model 1 consists of " << nChains << " Chains\n";
    //Pass through and add bonds corresponding to peptide bonds
    for (int iChain = 0; iChain < nChains; iChain ++){
        mmdb::Chain* chain = model->GetChain(iChain);
        int nResidues;
        mmdb::Residue**   residues;
        chain->GetResidueTable(residues, nResidues);
        // std::cout << "debug:: identifyBonds(): Chain " << chain->GetChainID() << " consists of " << nResidues << " Residues\n";
        mmdb::Residue* lastResidue = 0;
        for (int iResidue = 0; iResidue < nResidues; iResidue++){
            mmdb::Residue* residue = chain->GetResidue(iResidue);
            if (residue->isAminoacid()){
                if (lastResidue != 0){
                    mmdb::Atom* CA_i = residue->GetAtom("CA", " C", "*");
                    mmdb::Atom* CA_i_minus_1 = lastResidue->GetAtom("CA", " C", "*");
                    if (CA_i != 0 && CA_i_minus_1 != 0) {
                        FCXXCoord Coord_CA_i( CA_i->x, CA_i->y, CA_i->z);
                        FCXXCoord Coord_CA_i_minus_1( CA_i_minus_1->x, CA_i_minus_1->y, CA_i_minus_1->z);
                        FCXXCoord delta = (Coord_CA_i-Coord_CA_i_minus_1);
                        if(delta.get3DLength()<4.1){
                            mmdb::Atom* N_i = residue->GetAtom("N", " N", "*");
                            mmdb::Atom* C_i_minus_1 = lastResidue->GetAtom("C", " C", "*");
                            if(N_i!=0 && C_i_minus_1!=0){
                                N_i->AddBond(C_i_minus_1, 1, 1);
                                C_i_minus_1->AddBond(N_i, 1, 1);
                            }
                        }
                    }
                }
                lastResidue = residue;
            }
        }
    }
    return 0;
}

FCXXCoord MyMolecule::getCentre()
{
    return centreOfSelectionString(std::string("*/*/*/*"));
}

FCXXCoord MyMolecule::centreOfSelectionString(std::string selectionString)
{
    CompoundSelection selection(selectionString);
    int selHnd = selection.handleInMMDB(mmdb);
    FCXXCoord result = centreOfSelectionHandle(selHnd);
    mmdb->DeleteSelection(selHnd);
    return result;
}

FCXXCoord MyMolecule::centreOfSelectionHandle(int selHnd)
{
    mmdb::AtomStat sAtomStat;
    mmdb->GetAtomStatistics(selHnd, sAtomStat);
    return FCXXCoord (sAtomStat.xm, sAtomStat.ym, sAtomStat.zm);
}

void MyMolecule::writePDB(const std::string &filePath)
{
    mmdb->WritePDBASCII(filePath.c_str());
}

void MyMolecule::setResidueRadii(const std::map<std::tuple<std::string, int, std::string>, std::tuple<float, float, float>> &radii) {
    residueRadii = radii;
}

std::tuple<float, float, float> MyMolecule::getRadiiForResidue(mmdb::Residue *res) const {
    if (!res) return std::make_tuple(1.0f, 1.0f, 1.0f);
    std::string chain_id = res->GetChainID();
    int res_no = res->GetSeqNum();
    std::string ins_code = res->GetInsCode();
    auto key = std::make_tuple(chain_id, res_no, ins_code);
    auto it = residueRadii.find(key);
    if (it != residueRadii.end()) {
        return it->second;
    }
    return std::make_tuple(1.0f, 1.0f, 1.0f);
}

std::ostream& operator<<(std::ostream& o, const MyMolecule &myMolecule)
{
    o<< "Original name:"<<myMolecule.getPDBCode()<<"\n"<< "nAtoms:" << myMolecule.getMmdb()->GetNumberOfAtoms("/*/*/*.*/*:*");
    return o;
}
