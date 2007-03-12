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
