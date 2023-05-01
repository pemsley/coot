/*
 *  CXXChargeTable.cpp
 *  lpbSolver
 *
 *  Created by gruber on Mon Jul 12 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include "CXXChargeTable.h"

CXX_mot::CXXChargeTable::CXXChargeTable() {
	
	map<string, stringDoubleMap> atomNameToChargeMap;
	
	// this is ugly - should be able to modify these in situ ....
	// put in each empty map and then modify the inserted map !!

	stringDoubleMap HetMap;
	addChargeToResidueMap(HetMap," OT1", -0.5);
	addChargeToResidueMap(HetMap," OT2", -0.5);
	addChargeToResidueMap(HetMap," CU ", 2.0);
	addChargeToResidueMap(HetMap," ZN ", 2.0);
	addResidueMap("HET", HetMap);
	
	stringDoubleMap ALAMap;
	addChargeToResidueMap(ALAMap," OT1", -0.5);
	addChargeToResidueMap(ALAMap," OT2", -0.5);
	addResidueMap("ALA", ALAMap);

	stringDoubleMap CYSMap;
	addChargeToResidueMap(CYSMap," OT1", -0.5);
	addChargeToResidueMap(CYSMap," OT2", -0.5);
	addResidueMap("CYS", CYSMap);

	stringDoubleMap ASPMap;
	addChargeToResidueMap(ASPMap," OD1", -0.5);
	addChargeToResidueMap(ASPMap," OD2", -0.5);
	addChargeToResidueMap(ASPMap," OT1", -0.5);
	addChargeToResidueMap(ASPMap," OT2", -0.5);
	addResidueMap("ASP", ASPMap);

	stringDoubleMap GLUMap;
	addChargeToResidueMap(GLUMap," OE1", -0.5);
	addChargeToResidueMap(GLUMap," OE2", -0.5);
	addChargeToResidueMap(GLUMap," OT1", -0.5);
	addChargeToResidueMap(GLUMap," OT2", -0.5);
	addResidueMap("GLU", GLUMap);

	stringDoubleMap PHEMap;
	addChargeToResidueMap(PHEMap," OT1", -0.5);
	addChargeToResidueMap(PHEMap," OT2", -0.5);
	addResidueMap("PHE", PHEMap);

	stringDoubleMap GLYMap;
	addChargeToResidueMap(GLYMap," OT1", -0.5);
	addChargeToResidueMap(GLYMap," OT2", -0.5);
	addResidueMap("GLY", GLYMap);

	stringDoubleMap HISMap;
	addChargeToResidueMap(HISMap," OT1", -0.5);
	addChargeToResidueMap(HISMap," OT2", -0.5);
	addResidueMap("HIS", HISMap);

	stringDoubleMap ILEMap;
	addChargeToResidueMap(ILEMap," OT1", -0.5);
	addChargeToResidueMap(ILEMap," OT2", -0.5);
	addResidueMap("ILE", ILEMap);

	stringDoubleMap LYSMap;
	addChargeToResidueMap(LYSMap," NZ ", 1.0);
	addChargeToResidueMap(LYSMap," OT1", -0.5);
	addChargeToResidueMap(LYSMap," OT2", -0.5);
	addResidueMap("LYS", LYSMap);

	stringDoubleMap LEUMap;
	addChargeToResidueMap(LEUMap," OT1", -0.5);
	addChargeToResidueMap(LEUMap," OT2", -0.5);
	addResidueMap("LEU", LEUMap);

	stringDoubleMap METMap;
	addChargeToResidueMap(METMap," OT1", -0.5);
	addChargeToResidueMap(METMap," OT2", -0.5);
	addResidueMap("MET", METMap);

	stringDoubleMap ASNMap;
	addChargeToResidueMap(ASNMap," OT1", -0.5);
	addChargeToResidueMap(ASNMap," OT2", -0.5);
	addResidueMap("ASN", ASNMap);

	stringDoubleMap PROMap;
	addChargeToResidueMap(PROMap," OT1", -0.5);
	addChargeToResidueMap(PROMap," OT2", -0.5);
	addResidueMap("PRO", PROMap);

	stringDoubleMap GLNMap;
	addChargeToResidueMap(GLNMap," OT1", -0.5);
	addChargeToResidueMap(GLNMap," OT2", -0.5);
	addResidueMap("GLN", GLNMap);

	stringDoubleMap ARGMap;
	addChargeToResidueMap(ARGMap," NH1", 0.5);
	addChargeToResidueMap(ARGMap," NH2", 0.5);
	addChargeToResidueMap(ARGMap," OT1", -0.5);
	addChargeToResidueMap(ARGMap," OT2", -0.5);
	addResidueMap("ARG", ARGMap);

	stringDoubleMap SERMap;
	addChargeToResidueMap(SERMap," OT1", -0.5);
	addChargeToResidueMap(SERMap," OT2", -0.5);
	addResidueMap("SER", SERMap);

	stringDoubleMap THRMap;
	addChargeToResidueMap(THRMap," OT1", -0.5);
	addChargeToResidueMap(THRMap," OT2", -0.5);
	addResidueMap("THR", THRMap);

	stringDoubleMap VALMap;
	addChargeToResidueMap(VALMap," OT1", -0.5);
	addChargeToResidueMap(VALMap," OT2", -0.5);
	addResidueMap("VAL", VALMap);

	stringDoubleMap TRPMap;
	addChargeToResidueMap(TRPMap," OT1", -0.5);
	addChargeToResidueMap(TRPMap," OT2", -0.5);
	addResidueMap("TRP", TRPMap);

	stringDoubleMap TYRMap;
	addChargeToResidueMap(TYRMap," OT1", -0.5);
	addChargeToResidueMap(TYRMap," OT2", -0.5);
	addResidueMap("TYR", TYRMap);
	
}



int CXX_mot::CXXChargeTable::addResidueMap (string theResidueName, stringDoubleMap theAtomCharges) {
	
	atomNameToChargeMap.insert(make_pair(theResidueName, theAtomCharges));
	return 0;
	
}

int CXX_mot::CXXChargeTable::addChargeToResidueMap (stringDoubleMap& theResidueMap, string theAtomName, double theCharge) {
	
	
	theResidueMap.insert(make_pair(theAtomName, theCharge));
	return 0;	
}



double CXX_mot::CXXChargeTable::getCharge(string residueName, string atomName) {
	
	double charge = 0;
	stringStringDoubleMap::iterator positionInChargeMap;
	stringDoubleMap::iterator positionInResidueMap;
	
	positionInChargeMap = atomNameToChargeMap.find(residueName);
	if (positionInChargeMap != atomNameToChargeMap.end()) {

		positionInResidueMap = positionInChargeMap->second.find(atomName);
		if (positionInResidueMap != positionInChargeMap->second.end()) {
			
			charge = positionInResidueMap->second;
		}
	}

	return charge;
}

