/*
 *  CompoundSelection.mm
 *  Aesop
 *
 *  Created by Martin Noble on 16/03/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "CompoundSelection.h"
#include "MoleculesToTriangles/CXXSurface/TokenIterator.h"
#include <iostream>
#include <algorithm>
#include <string>
#include <map>

CompoundSelection::~CompoundSelection(){
    std::vector<std::pair<mmdb::SELECTION_KEY, SelectionPrimitive *> >::iterator pair = pairs.begin();
    for (;pair != pairs.end(); ++pair){
        if(CompoundSelection *p = dynamic_cast<CompoundSelection *>(pair->second)){
            delete p;
        }
        else if(MMDBStringPrimitive *p = dynamic_cast<MMDBStringPrimitive *>(pair->second)){
            delete p;
        }
        else if(MMDBSecondaryTypePrimitive *p = dynamic_cast<MMDBSecondaryTypePrimitive *>(pair->second)){
            delete p;
        }
        else if(MMDBSubsetTypePrimitive *p = dynamic_cast<MMDBSubsetTypePrimitive *>(pair->second)){
            delete p;
        }
    }
}

int CompoundSelection::handleInMMDB(mmdb::Manager *mmdb){
    int mmdbHandle = mmdb->NewSelection();

    std::cout << "In CompoundSelection::handleInMMDB " << " " << mmdbHandle << " " << selectionString << std::endl;
    
    std::vector<std::pair<mmdb::SELECTION_KEY, SelectionPrimitive *> >::iterator pair = pairs.begin();
    while (pair != pairs.end()){
        int primitiveHandle = pair->second->handleInMMDB(mmdb);
        mmdb->Select(mmdbHandle, mmdb::STYPE_ATOM, primitiveHandle, pair->first);
        mmdb->DeleteSelection(primitiveHandle);
        pair++;
    }
    invertIfNeccessary(mmdbHandle, mmdb);
    return mmdbHandle;
}

inline void SelectionPrimitive::invertIfNeccessary(int &mmdbHandle, mmdb::Manager *mmdb){
    const char *all = "/*/*/*/*";
    if (invert) {
        int allAtoms = mmdb->NewSelection();
        //std::cout << "In invert " << allAtoms << std::endl;
        mmdb->Select(allAtoms, mmdb::STYPE_ATOM, all, mmdb::SKEY_NEW);
        mmdb->Select(allAtoms, mmdb::STYPE_ATOM, mmdbHandle, mmdb::SKEY_CLR);
        mmdb->DeleteSelection(mmdbHandle);
        mmdbHandle = allAtoms;
    }
}

CompoundSelection::CompoundSelection(std::string _selectionText){
    name = std::string("blankname");
    setSelectionString(_selectionText);
}

CompoundSelection::CompoundSelection(std::string _selectionText, std::string _name){
    setSelectionString(_selectionText);
    setName(_name);
}

std::string CompoundSelection::trimString(std::string _aString){
    size_t firstNonSpace = _aString.find_first_not_of(' ');
    size_t lastNonSpace = _aString.find_last_not_of(' ');
    std::string result = "";
    if (firstNonSpace != std::string::npos && lastNonSpace != std::string::npos){
        result = _aString.substr(firstNonSpace, (lastNonSpace-firstNonSpace)+1);
    }
    return result;
}

void CompoundSelection::setSelectionString(const std::string &_selectionString) {
    selectionString = trimString(_selectionString);
    //std::cout << "Parsing compoundSelectionString with text" <<selectionString << std::endl;
    invert = false;
    int subClauseLevel = 0;
    std::string subClause = "";
    //Inch forward through string
    bool doInvert = false;
    bool inCID = false;
    mmdb::SELECTION_KEY combineRule = mmdb::SKEY_NEW;
    for (auto cIter = selectionString.begin(); cIter != selectionString.end(); ++cIter){
        if (subClauseLevel > 0){
            //Check for termination
            if (*cIter == '}'){
                //std::cout << "Finishing subclause level" << subClauseLevel <<" " << subClause << std::endl;
                subClauseLevel--;
                if (subClauseLevel == 0){
                    CompoundSelection *primitive = new CompoundSelection(subClause);
                    primitive->setInvert(doInvert);
                    pairs.push_back(std::pair<mmdb::SELECTION_KEY, SelectionPrimitive *>(combineRule, primitive));
                    doInvert = false;
                    subClause = "";
                }
                else subClause += *cIter;
            }
            else if (*cIter == '{'){
                subClauseLevel++;
                subClause += *cIter;
            }
            else subClause += *cIter;
        }
        else if (inCID){
            //Check for termination
            if (*cIter == '&' || *cIter == '|' || *cIter == '\0' || *cIter == ' '){
                //std::cout << "Reached subclause terminator :" << subClause << std::endl;
                auto subString = trimString(subClause);
                std::string upperCaseSubString = subString;
                if (MMDBSecondaryTypePrimitive::secondaryTypes.find(subString) !=
                    MMDBSecondaryTypePrimitive::secondaryTypes.end() ){
                    MMDBSecondaryTypePrimitive *primitive = new MMDBSecondaryTypePrimitive(subString);
                    primitive->setInvert(doInvert);
                    pairs.push_back(std::pair<mmdb::SELECTION_KEY, SelectionPrimitive *>(combineRule, primitive));
                    doInvert = false;
                }
                else if (MMDBSubsetTypePrimitive::subsetTypes.find(upperCaseSubString) !=
                         MMDBSubsetTypePrimitive::subsetTypes.end() ){
                    MMDBSubsetTypePrimitive *primitive = new MMDBSubsetTypePrimitive(upperCaseSubString);
                    primitive->setInvert(doInvert);
                    pairs.push_back(std::pair<mmdb::SELECTION_KEY, SelectionPrimitive *>(combineRule, primitive));
                    doInvert = false;
                }
                else {
                    MMDBStringPrimitive *primitive = new MMDBStringPrimitive(subString);
                    primitive->setInvert(doInvert);
                    pairs.push_back(std::pair<mmdb::SELECTION_KEY, SelectionPrimitive *>(combineRule, primitive));
                    doInvert = false;
                }
                inCID = false;
            }
            else {
                subClause += *cIter;
            }
        }
        else if (*cIter == '{'){
            subClause = "";
            subClauseLevel++;
        }
        else if (*cIter == '!'){
            doInvert = true;
        }
        else if (*cIter == '&'){
            combineRule = mmdb::SKEY_AND;
        }
        else if (*cIter == '|'){
            combineRule = mmdb::SKEY_OR;
        }
        else if (*cIter == ' ' || *cIter == '\t' || *cIter == '\0' || *cIter == '\n' || *cIter == '\r'){
        }
        else {
            subClause = *cIter;
            inCID=true;
        }
    }
    if (trimString(subClause).length() > 0){
        auto subString = trimString(subClause);
        std::string upperCaseSubString = subString;
        if (MMDBSecondaryTypePrimitive::secondaryTypes.find(subString) !=
            MMDBSecondaryTypePrimitive::secondaryTypes.end() ){
            MMDBSecondaryTypePrimitive *primitive = new MMDBSecondaryTypePrimitive(subString);
            primitive->setInvert(doInvert);
            pairs.push_back(std::pair<mmdb::SELECTION_KEY, SelectionPrimitive *>(combineRule, primitive));
            doInvert = false;
        }
        else if (MMDBSubsetTypePrimitive::subsetTypes.find(upperCaseSubString) !=
                 MMDBSubsetTypePrimitive::subsetTypes.end() ){
            MMDBSubsetTypePrimitive *primitive = new MMDBSubsetTypePrimitive(upperCaseSubString);
            primitive->setInvert(doInvert);
            pairs.push_back(std::pair<mmdb::SELECTION_KEY, SelectionPrimitive *>(combineRule, primitive));
            doInvert = false;
        }
        else {
            MMDBStringPrimitive *primitive = new MMDBStringPrimitive(subString);
            primitive->setInvert(doInvert);
            pairs.push_back(std::pair<mmdb::SELECTION_KEY, SelectionPrimitive *>(combineRule, primitive));
            doInvert = false;
        }
    }
}

void CompoundSelection::describe()
{
    std::cout << "Compound selection with following subclauses : " << std::endl;
    std::vector<std::pair<mmdb::SELECTION_KEY, SelectionPrimitive *> >::iterator pair = pairs.begin();
    while (pair != pairs.end()){
        std::string ruleString;
        switch (pair->first) {
            case mmdb::SKEY_NEW:
                ruleString = "SKEY_NEW : ";
                break;
            case mmdb::SKEY_OR:
                ruleString = "SKEY_OR  : ";
                break;
            case mmdb::SKEY_AND:
                ruleString = "SKEY_AND : ";
                break;
            default:
                break;
        };
        std::cout << ruleString;
        pair->second->describe();
        pair++;
    }
}

void CompoundSelection::deleteInMMDB(mmdb::Manager *mmdb)
{
    //  get new selection handle
    auto selHnd = handleInMMDB(mmdb);
    
    //  get selected set of atoms
    int nAtoms;
    mmdb::Atom** atom;
    mmdb->GetSelIndex ( selHnd,atom,nAtoms );
    std::cout << "Selection contained "<<nAtoms<<std::endl;;
    
    //  deletion loop over all selected atoms.
    //  NOTE that the atom array itself must
    //  be neither deleted nor changed!
    for (int i=0;i<nAtoms;i++)
        delete atom[i];
    
    //  update internal references
    mmdb->FinishStructEdit();
}


MMDBStringPrimitive::MMDBStringPrimitive(std::string _selectionString){
    selectionString = _selectionString;
}

int MMDBStringPrimitive::handleInMMDB(mmdb::Manager *mmdb){
    int result = mmdb->NewSelection();
    //std::cout << "In MMDBSelection::handleInMMDB " << " " << result << " " << selectionString << std::endl;
    mmdb->Select(result, mmdb::STYPE_ATOM, selectionString.c_str(), mmdb::SKEY_NEW);
    invertIfNeccessary(result, mmdb);
    return result;
}

void MMDBStringPrimitive::describe()
{
    std::cout << (invert?"not ":"") << "MMDBStringPrimitive with selection Text -" << selectionString << "-\n";
}

void CompoundSelection::updateUsingMMDBSelection(mmdb::Manager *mmdb, int selHnd){
    std::string newString;// = "{";
    int nModels;
    mmdb::PPModel modelTable;
    mmdb->GetModelTable(modelTable, nModels);
    
    bool firstModel = true;
    for (int iModel = 0; iModel < nModels; iModel++){
        mmdb::Model *model = modelTable[iModel];
        if (modelHandleIntersection(selHnd, model) != IsDisjoint){
            if (!firstModel) newString += " | ";
            firstModel = false;
            
            char modelStringBuffer[256];
            sprintf(modelStringBuffer, " { /%d/*/*.*/*:* } & {", iModel+1);
            std::string modelString(modelStringBuffer);
            newString += modelString;
            
            mmdb::Chain** chainTable;
            int nChains;
            
            bool firstChain = true;
            model->GetChainTable(chainTable, nChains);
            for (int iChain = 0; iChain < nChains; iChain++){
                mmdb::Chain *chain = chainTable[iChain];
                if (chainHandleIntersection(selHnd, chain) != IsDisjoint){
                    if (!firstChain) newString += " | ";
                    firstChain = false;
                    
                    char chainStringBuffer[256];
                    sprintf(chainStringBuffer, "{ { /*/%s/*.*/*:* } & { ", chain->GetChainID());
                    std::string chainString(chainStringBuffer);
                    newString += chainString;
                    
                    mmdb::Residue** residueTable;
                    int nResidues;
                    chain->GetResidueTable(residueTable, nResidues);
                    
                    bool inSection = false;
                    int startResidue=0, stopResidue=0;
                    
                    bool firstResidueRange = true;
                    for (int iResidue = 0; iResidue < nResidues; iResidue++){
                        mmdb::Residue*residue = residueTable[iResidue];
                        if (residueHandleIntersection(selHnd, residue) != IsDisjoint) {
                            stopResidue = residue->GetSeqNum();
                            if (!inSection) {
                                startResidue = residue->GetSeqNum();
                            }
                            inSection = true;
                        }
                        else {
                            if (inSection) { //Handle section termination
                                if (!firstResidueRange) newString += " | ";
                                firstResidueRange = false;
                                
                                char residueStringBuffer[256];
                                sprintf(residueStringBuffer, " { /*/*/%d-%d/*:* } ", startResidue, stopResidue);
                                std::string residueString(residueStringBuffer);
                                newString += residueString;
                            }
                            inSection = false;
                        }
                    }
                    if (inSection){ //Handle section termination
                        if (!firstResidueRange) newString += " | ";
                        char residueStringBuffer[256];
                        sprintf(residueStringBuffer, " { /*/*/%d-%d/*:* } ", startResidue, stopResidue);
                        std::string residueString(residueStringBuffer);
                        newString += residueString;
                    }
                    newString += " }  } ";
                }
            }
            newString += " } ";
        }
    }
    //newString += " } ";
    //std::cout << newString;
    setSelectionString(newString);
}

std::pair<std::string, int> MMDBSecondaryTypePrimitive::secondaryTypesArray[] =
{
    std::pair<std::string,int>(std::string("SSE_None"), mmdb::SSE_None),
    std::pair<std::string,int>(std::string("SSE_Helix"), mmdb::SSE_Helix),
    std::pair<std::string,int>(std::string("SSE_Strand"), mmdb::SSE_Strand)
};
std::map<std::string, int> MMDBSecondaryTypePrimitive::secondaryTypes =
std::map<std::string, int> (MMDBSecondaryTypePrimitive::secondaryTypesArray,
                            MMDBSecondaryTypePrimitive::secondaryTypesArray +
                            (sizeof(MMDBSecondaryTypePrimitive::secondaryTypesArray)/sizeof(std::pair<std::string, int>)));

int MMDBSecondaryTypePrimitive::handleInMMDB(mmdb::Manager *mmdb){
    int handle = mmdb->NewSelection();
    //std::cout << "In SecondarySelection::handleInMMDB " << " " << handle << " " << selectionString << std::endl;
    mmdb::Atom **atoms;
    int nAtoms;
    mmdb->GetAtomTable(atoms, nAtoms);
    for (int i=0; i<nAtoms; i++){
        mmdb::Atom* atom = atoms[i];
        if (atom){
            if (atom->GetResidue()->SSE == type){
                mmdb->SelectAtom(handle, atom, mmdb::SKEY_OR, false);
            }
        }
    }
    mmdb->MakeSelIndex(handle);
    return handle;
}

void MMDBSecondaryTypePrimitive::describe()
{
    std::cout << "MMDBSecondaryTypePrimitive with selection Text -" << selectionString << "-\n";
}


std::pair<std::string, std::string> MMDBSubsetTypePrimitive::subsetTypesArray[] =
{
    std::pair<std::string,std::string>(std::string("MAIN"), std::string("/*/*/*.*/N,CA,C,O,H")),
    std::pair<std::string,std::string>(std::string("SIDE"), std::string("/*/*/*.*/!N,C,O,H")),
    std::pair<std::string,std::string>(std::string("WATER"), std::string("/*/*/(WAT,HOH,OH2,H2O)")),
    std::pair<std::string,std::string>(std::string("MONOMERS"), std::string("/*/*/(!ALA,CYS,ASP,GLU,PHE,GLY,HIS,ILE,LYS,LEU,MET,ASN,PRO,GLN,ARG,SER,THR,VAL,TRP,TYR,WAT,HOH,THP,SEP,TPO,TYP,PTR,OH2,H2O)")),
    std::pair<std::string,std::string>(std::string("AMINOACIDS"), std::string("/*/*/(ALA,CYS,ASP,GLU,PHE,GLY,HIS,ILE,LYS,LEU,MET,ASN,PRO,GLN,ARG,SER,THR,VAL,TRP,TYR,THP,SEP,TPO,TYP,PTR,MSE)")),
    std::pair<std::string,std::string>(std::string("NUCLEICACIDS"), std::string("/*/*/(DG,DA,DC,DT,DU,A,G,T,C,U)")),
    std::pair<std::string,std::string>(std::string("ALL"), std::string("/*/*/*.*/*:*")),
    
};
std::map<std::string, std::string> MMDBSubsetTypePrimitive::subsetTypes =
std::map<std::string, std::string> (MMDBSubsetTypePrimitive::subsetTypesArray,
                                    MMDBSubsetTypePrimitive::subsetTypesArray +
                                    (sizeof(MMDBSubsetTypePrimitive::subsetTypesArray)/sizeof(std::pair<std::string, std::string>)));



int MMDBSubsetTypePrimitive::handleInMMDB(mmdb::Manager *mmdb){
    int result = mmdb->NewSelection();
    //std::cout << "In MMDBSelection::handleInMMDB " << " " << result << " " << selectionString << std::endl;
    mmdb->Select(result, mmdb::STYPE_ATOM, selectionString.c_str(), mmdb::SKEY_NEW);
    invertIfNeccessary(result, mmdb);
    return result;
}

void MMDBSubsetTypePrimitive::describe()
{
    std::cout << "MMDBSubsetTypePrimitive with selection Text -" << selectionString << "-\n";
}

