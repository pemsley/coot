/* 
 * 
 * Copyright 2004 by The University of Oxford
 * Author: Martin Noble, Jan Gruber
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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
#include "CXXUtils.h"

int CXXUtils::assignUnitedAtomRadius  (PCMMDBManager theManager, int selHnd) {
	 // Add a radius property to the atoms
	 int iRadiusHandle = theManager->RegisterUDReal(UDR_ATOM, "PerAtomRadius");
	 if (!iRadiusHandle) {
		 printf ( " registration failed.\n" );
		 exit ( 1 );
	 }
	 if (iRadiusHandle==UDDATA_WrongUDRType) {
		 printf ( " wrong registration type used.\n" );
		 exit ( 2 );
	 }
	 
	 CAtom **SelAtom;
	 int nSelAtoms;
	 theManager->GetSelIndex(selHnd, SelAtom, nSelAtoms);
	 
	 for (int iAtom = 0; iAtom < nSelAtoms; iAtom++) {
		 CAtom *anAtom = SelAtom[iAtom];
		 string atomName(anAtom->name);
		 string residueName(anAtom->residue->name);
		 int unassigned = 1;
		 double radius = 1.8;
		 for (int iAtomRadius = 0; iAtomRadius < CXXUtils::nAtomRadii && unassigned; iAtomRadius++){
			if (atomName == CXXUtils::unitedAtomRadii[iAtomRadius].atomName){
				if (residueName == string(CXXUtils::unitedAtomRadii[iAtomRadius].residueName) ||
					"*  "== string(CXXUtils::unitedAtomRadii[iAtomRadius].residueName)){
					unassigned = 0;
					radius = CXXUtils::unitedAtomRadii[iAtomRadius].radius;
				}
			}
		 }
		 anAtom->PutUDData(iRadiusHandle, radius);
		 if (unassigned) cout << "Not able to match atom :" << atomName << ":" << residueName << ":\n";
	 }
	 return 0;
}

void CXXUtils::reformatAtomRadii(){
	for (int iAtomRadius = 0; iAtomRadius < CXXUtils::nAtomRadii; iAtomRadius++){
		string paddedAtomName, paddedResidueName;
		string testAtomName(CXXUtils::unitedAtomRadii[iAtomRadius].atomName);
		string testResidueName(CXXUtils::unitedAtomRadii[iAtomRadius].residueName);
		float &testAtomRadius(CXXUtils::unitedAtomRadii[iAtomRadius].radius);
		
		
		cout << "{ \"";
		int atomNameLength=0;
		if (testAtomName.substr(0,1)=="H" ||
			testAtomName.substr(0,1)=="D" ||
			testAtomName.substr(0,1)=="T" 
			) {
		}
		else {
			cout << " ";
			atomNameLength++;
		}
		cout << testAtomName;
		atomNameLength+=testAtomName.length();
		while (atomNameLength < 4){
			cout << " ";
			atomNameLength++;
		}
		cout << "\", \"";
		
		int residueNameLength=0;
		cout << testResidueName;
		residueNameLength += testResidueName.length();
		while (residueNameLength < 3){
			cout << " ";
			residueNameLength++;
		}
		cout << "\", ";
		cout << testAtomRadius;
		cout << "},\n";
	}
}
