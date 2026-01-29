/*
 * MoleculesToTriangles/CXXClasses/MyMolecule.h
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

#ifndef MyMolecule_h
#define MyMolecule_h
#include <memory>
#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"

#include <vector>
#include <string>
#include <map>
#include <tuple>

class DiscreteSegment;
class DisplayPrimitive;
class ObjectSelection;
class DishyBaseContainer_t;

#include "mmdb2/mmdb_manager.h"

enum SecondaryStructureUsageType { USE_HEADER_INFO, DONT_USE, CALC_SECONDARY_STRUCTURE };

typedef struct _PdbCoords {
    float xyz[4];
} PdbCoords;
typedef struct _AtomCard {
    int selected;
    char atname[8];
    char alt[4];
    char restype[8];
    char chainid[4];
    char resname[8];
    PdbCoords crd; 
    float properties[4];
} AtomCard;

class MyMolecule {
private:
    bool doDraw;
    bool ownsMMDB;
    std::map<std::tuple<std::string, int, std::string>, float> residueRadii;
    int loadFromPDB(const char *filepath, int secondaryStructureUsageFlag);
    int processCoords(int secondaryStructureUsageFlag);
public:
    mmdb::Manager *mmdb;
	std::string PDBCode;
    MyMolecule();
    MyMolecule(const char *filePath, int secondaryStructureUsageFlag);
    MyMolecule(std::string filePathString, int secondaryStructureUsageFlag);
    MyMolecule(mmdb::Manager *fromMMDBManager, int secondaryStructureUsageFlag) : doDraw(true), mmdb(fromMMDBManager) {
        ownsMMDB = false;
        processCoords(secondaryStructureUsageFlag);
    };
    ~MyMolecule();
    int loadCoords( char *fileData, int length, int secondaryStructureUsageFlag);
    int identifySegments(std::vector<DiscreteSegment *> &segments, int selHnd);
    int identifyDishyBases(std::map<mmdb::Chain *, DishyBaseContainer_t> &dishy_bases_chain_map, int selHnd);
    int identifyBonds();
    FCXXCoord getCentre();
    FCXXCoord centreOfSelectionString(std::string selectionString);
    FCXXCoord centreOfSelectionHandle(int selHnd);
    mmdb::Manager *getMmdb() const{
        return mmdb;
    };
	void setPDBCode(std::string _code){
		PDBCode = _code;
	};
	const std::string getPDBCode() const{
		return PDBCode;
	};
    bool getDoDraw() const {
        return doDraw;
    };
    void setDoDraw(const bool &yesOrNo) {
        doDraw = yesOrNo;
    };
    int FormatPDBCard (AtomCard theAtom, char *card,int count);
    void writePDB(const std::string &filePath);
    void setResidueRadii(const std::map<std::tuple<std::string, int, std::string>, float> &radii);
    float getRadiusForResidue(mmdb::Residue *res, float default_radius = 1.0f) const;
    static std::shared_ptr<MyMolecule> create(std::string(filePathString), int secondaryStructureUsageFlag){
        auto newMolecule = new MyMolecule(filePathString, secondaryStructureUsageFlag);
        std::cout << newMolecule;
        return std::shared_ptr<MyMolecule>(newMolecule);
    };
    static std::shared_ptr<MyMolecule> createFromString(std::string(contentsAsString)){
        auto newMolecule = new MyMolecule();
        //20240912-PE maybe the secondary structure usage status should be passed to this function
        // (I don't know what this function is or if it is used)
        newMolecule->loadCoords(&(contentsAsString[0]), contentsAsString.size(), CALC_SECONDARY_STRUCTURE);
        return std::shared_ptr<MyMolecule>(newMolecule);
    };

};
std::ostream& operator<<(std::ostream&,const MyMolecule &);

#endif
