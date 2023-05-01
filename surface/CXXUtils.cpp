#include "CXXUtils.h"
#include <stdlib.h>
#include <string.h>

int CXXUtils_old::assignCharge(mmdb::PManager theManager, int selHnd,
			       CXX_mot::CXXChargeTable *theChargeTable) {
	mmdb::Atom **SelAtom;
	int nSelAtoms;
	theManager->GetSelIndex(selHnd, SelAtom, nSelAtoms);
	
	//Assign atom charges
	for (int iAtom = 0; iAtom < nSelAtoms; iAtom++) {
	   mmdb::Atom *theAtom = SelAtom[iAtom];
	   string atomName(theAtom->name);
	   string residueName(theAtom->residue->name);
	   double theAtomCharge;
	   theAtomCharge = theChargeTable->getCharge(residueName, atomName);		
	   theAtom->charge = theAtomCharge;
	}
	return 0;
}

int CXXUtils_old::assignUnitedAtomRadius  (mmdb::PManager theManager, int selHnd) {
	// Add a radius property to the atoms
	int iRadiusHandle = theManager->RegisterUDReal(mmdb::UDR_ATOM, "PerAtomRadius");
	if (!iRadiusHandle) {
	   printf ( "ERROR:: registration failed.\n" );
	   // exit ( 1 );  // no exit from libraries
	   return -100;
	}
	if (iRadiusHandle==mmdb::UDDATA_WrongUDRType) {
	   printf ( "ERROR:: wrong registration type used.\n" );
	   // exit ( 2 );
	   return -100;
	}
	
	mmdb::Atom **SelAtom;
	int nSelAtoms;
	theManager->GetSelIndex(selHnd, SelAtom, nSelAtoms);
	
	for (int iAtom = 0; iAtom < nSelAtoms; iAtom++) {
		mmdb::Atom *anAtom = SelAtom[iAtom];
		string atomName(anAtom->name);
		string residueName(anAtom->residue->name);
		int unassigned = 1;
		double radius = 1.8;
		for (int iAtomRadius = 0; iAtomRadius < CXXUtils_old::nAtomRadii && unassigned; iAtomRadius++){
			if (atomName == CXXUtils_old::unitedAtomRadii[iAtomRadius].atomName){
				if (residueName == string(CXXUtils_old::unitedAtomRadii[iAtomRadius].residueName) ||
					"*  "== string(CXXUtils_old::unitedAtomRadii[iAtomRadius].residueName)){
					unassigned = 0;
					radius = CXXUtils_old::unitedAtomRadii[iAtomRadius].radius;
				}
			}
		}
		anAtom->PutUDData(iRadiusHandle, radius);
		//		 if (unassigned) cout << "Not able to match atom :" << atomName << ":" << residueName << ":\n";
	}
	return 0;
}

void CXXUtils_old::reformatAtomRadii(){
	for (int iAtomRadius = 0; iAtomRadius < CXXUtils_old::nAtomRadii; iAtomRadius++){
		string paddedAtomName, paddedResidueName;
		string testAtomName(CXXUtils_old::unitedAtomRadii[iAtomRadius].atomName);
		string testResidueName(CXXUtils_old::unitedAtomRadii[iAtomRadius].residueName);
		float &testAtomRadius(CXXUtils_old::unitedAtomRadii[iAtomRadius].radius);
		
		
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

int CXXUtils_old::selectionStringToSelHnd(mmdb::PManager allAtomsManager_in, std::string selectionString, int existingSelection, mmdb::SELECTION_KEY selKeyRequest){
   int selHnd;
   mmdb::SELECTION_KEY selKey;
	if (existingSelection == -1) {
		selHnd = allAtomsManager_in->NewSelection();
	}
	else {
		selHnd = existingSelection;
	}
	selKey = selKeyRequest;
	char *pstring = (char *) malloc (sizeof(selectionString.c_str())+1);
	strcpy (pstring, selectionString.c_str());
	allAtomsManager_in->Select ( selHnd, mmdb::STYPE_ATOM, pstring, selKey);
	mmdb::PPAtom SelAtoms; 
	int nSelAtoms;
	allAtomsManager_in->GetSelIndex(selHnd, SelAtoms, nSelAtoms);
	cout << "Selection now contains " << nSelAtoms << " atoms\n";
	free (pstring);
	return selHnd;
}

int CXXUtils_old::unCharge(mmdb::PManager theManager, int selHnd){
	mmdb::Atom **SelAtom;
	int nSelAtoms;
	theManager->GetSelIndex(selHnd, SelAtom, nSelAtoms);
	
	//Assign atom charges
	for (int iAtom = 0; iAtom < nSelAtoms; iAtom++) {
		mmdb::Atom *theAtom = SelAtom[iAtom];
		theAtom->charge = 0.;
	}
	return 0;
}

