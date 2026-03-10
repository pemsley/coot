/*
 * MoleculesToTriangles/CXXClasses/CompoundSelection.h
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

#ifndef CompoundSelection_h
#define CompoundSelection_h

#include <vector>
#include <memory>
#include <map>
#include <utility>
#include <string>
#include <iostream>
#include "mmdb2/mmdb_manager.h"
#include "MyMolecule.h"

class SelectionPrimitive {
//Modelled
protected:
	std::string selectionString;
    std::string name;
public:
	int invert;
    virtual ~SelectionPrimitive() = default;
    const std::string getName() const {
        return
        name;
    };
    void setName(const std::string _name){
        name = _name;
    };

    virtual int handleInMMDB(mmdb::Manager *mmdb) = 0;
	virtual void describe() = 0;
	void setInvert(int _invert) {
		invert = _invert;
	};
    const std::string &getSelectionString() const {
        return selectionString;
    };
    virtual void setSelectionString(const std::string &_selectionString) = 0;
    virtual std::string format() = 0;
    void invertIfNeccessary(int &mmdbHandle, mmdb::Manager *mmdb);
};

class CompoundSelection : public SelectionPrimitive {
private:
public:
    static std::string trimString(std::string _aString);
    std::vector<std::pair<mmdb::SELECTION_KEY, SelectionPrimitive *> > pairs;
    CompoundSelection() {
        name = std::string("blankname");
    };
    CompoundSelection(std::string selectionText);
    CompoundSelection(std::string selectionText, std::string _name);
    virtual void setSelectionString(const std::string &_selectionString);
	virtual ~CompoundSelection();
    virtual int handleInMMDB(mmdb::Manager *mmdb);
    void deleteInMMDB(mmdb::Manager *mmdb);
	virtual void describe();
    
    /*static std::shared_ptr<CompoundSelection> create(std::string _selectionText){
        return std::shared_ptr<CompoundSelection>(new CompoundSelection(_selectionText));
    };
     */
    static std::shared_ptr<CompoundSelection> create(std::string _selectionText, std::string _name){
        return std::shared_ptr<CompoundSelection>(new CompoundSelection(_selectionText, _name));
    };
    
    enum Intersection {IsDisjoint, IsSubset, IsIntersector};
    static enum Intersection handleAtomsIntersection(int selHnd, mmdb::Atom **atoms, int nAtoms){
        bool noneIn = true;
        bool allIn = true;
        for (int iAtom = 0; iAtom < nAtoms && (allIn || noneIn); iAtom++){
            mmdb::Atom *atom = atoms[iAtom];
            if (!atom->isTer()){
                if (atoms[iAtom]->isInSelection(selHnd)) noneIn = false;
                else allIn = false;
            }
        }
        if (allIn) return IsSubset;
        else if (noneIn) return IsDisjoint;
        else return IsIntersector;
    };
    static enum Intersection residueHandleIntersection(int selHnd, mmdb::Residue *residue){
        int nAtoms;
        mmdb::Atom**  atomTable;
        residue->GetAtomTable(atomTable, nAtoms);
        return handleAtomsIntersection(selHnd, atomTable, nAtoms);
    };
    static enum Intersection chainHandleIntersection(int selHnd, mmdb::Chain *chain){
        bool noneIn = true;
        bool allIn = true;
        int nResidues;
        mmdb::Residue** resTable;
        chain->GetResidueTable(resTable, nResidues);
        for (int i=0; i<nResidues && (allIn || noneIn); i++){
            mmdb::Residue*residue = resTable[i];
            enum Intersection result= residueHandleIntersection(selHnd, residue);
            if (result != IsDisjoint) noneIn = false;
            if (result != IsSubset) allIn = false;
        }
        if (allIn) return IsSubset;
        else if (noneIn) return IsDisjoint;
        else return IsIntersector;
    };
    static enum Intersection modelHandleIntersection(int selHnd, mmdb::Model *model ){
        bool noneIn = true;
        bool allIn = true;
        int nChains;
        mmdb::Chain** chainTable;
        model->GetChainTable(chainTable, nChains);
        for (int i=0; i<nChains && (allIn || noneIn); i++){
            mmdb::Chain *chain = chainTable[i];
            enum Intersection result= chainHandleIntersection(selHnd, chain);
            if (result != IsDisjoint) noneIn = false;
            if (result != IsSubset) allIn = false;
        }
        if (allIn) return IsSubset;
        else if (noneIn) return IsDisjoint;
        else return IsIntersector;
    };
    void updateUsingMMDBSelection(mmdb::Manager *mmdb, int selHnd);
    virtual std::string format(){
        std::string combineRuleString;        
        std::vector<std::pair<mmdb::SELECTION_KEY, SelectionPrimitive *> >::iterator pair = pairs.begin();
        for (; pair != pairs.end(); ++pair){
            if (CompoundSelection *compoundSelection = dynamic_cast<CompoundSelection *>(pair->second)){
                if (compoundSelection->pairs.size() > 1) combineRuleString += std::string("{");
            }
            int combineRule = pair->first;
            if (combineRule == mmdb::SKEY_NEW) combineRuleString += std::string(" ");
            else if (combineRule == mmdb::SKEY_AND) combineRuleString += std::string(" & ");
            else if (combineRule == mmdb::SKEY_OR) combineRuleString += std::string(" | ");
            if (pair->second->invert) combineRuleString += std::string(" ! ");
            combineRuleString += pair->second->format();
            if (CompoundSelection *compoundSelection = dynamic_cast<CompoundSelection *>(pair->second)){
                if (compoundSelection->pairs.size() > 1) combineRuleString += std::string("}");
            }
        }
        return combineRuleString;
    };

};

class MMDBStringPrimitive : public SelectionPrimitive{
private:
public:
	MMDBStringPrimitive(std::string _selectionString);
	virtual int handleInMMDB(mmdb::Manager *mmdb);
	virtual void describe();
    virtual void setSelectionString(const std::string &_selectionString) {
        selectionString = _selectionString;
    };
    virtual std::string format(){
        return selectionString;  
    };
};

class MMDBSecondaryTypePrimitive : public SelectionPrimitive {
private:
    int type;
public:
    static std::pair<std::string, int> secondaryTypesArray[];
    static std::map<std::string, int> secondaryTypes;
	MMDBSecondaryTypePrimitive(std::string _selectionText){
        setSelectionString(_selectionText);
    };
    virtual int handleInMMDB(mmdb::Manager *mmdb);
    virtual void describe();
    virtual void setSelectionString(const std::string &_selectionString) {
        selectionString = _selectionString;
        type = secondaryTypes[selectionString];
    };
    virtual std::string format(){
        return selectionString;  
    };
};

class MMDBSubsetTypePrimitive : public SelectionPrimitive {
private:
    int type;
public:
    static std::pair<std::string, std::string> subsetTypesArray[];
    static std::map<std::string, std::string> subsetTypes;
	MMDBSubsetTypePrimitive(std::string _selectionText){
        setSelectionString(_selectionText);
    };
    virtual int handleInMMDB(mmdb::Manager *mmdb);
    virtual void describe();
    virtual void setSelectionString(const std::string &_selectionString) {
        name = _selectionString;
        selectionString = subsetTypes[_selectionString];
    };
    enum SubsetSelectionType {
        SubsetSelectionTypeMain, SubsetSelectionTypeSide
    };
    virtual std::string format(){
        return name;  
    };
};

#endif
