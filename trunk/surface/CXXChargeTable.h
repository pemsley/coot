/*
 *  CXXChargeTable.h
 *  lpbSolver
 *
 *  Created by gruber on Mon Jul 12 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CXXChargeTable_included
#define CXXChargeTable_included

#include <iostream>
#ifndef  __STRING_H
#include <string>
#endif

#include <map>
using namespace std;

typedef map<string, double> stringDoubleMap;
typedef map<string, stringDoubleMap> stringStringDoubleMap;

class CXXChargeTable {
	
private:
	
	stringStringDoubleMap atomNameToChargeMap;
		
	int addResidueMap (string residueName, stringDoubleMap theAtomCharges);
	int addChargeToResidueMap (stringDoubleMap& theMap, string theAtomName, double theCharge);
	
public:
		
	CXXChargeTable();
	double getCharge(string residueName, string atomName);
	
	
};

#endif
