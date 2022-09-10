#include "CXXUtils.h"
#include <stdlib.h>
#include <string.h>

int CXXUtils::assignCharge(mmdb::Manager* theManager, int selHnd, CXXChargeTable *theChargeTable){
	mmdb::Atom* *SelAtom;
	int nSelAtoms;
	theManager->GetSelIndex(selHnd, SelAtom, nSelAtoms);
	
	//Assign atom charges
	for (int iAtom = 0; iAtom < nSelAtoms; iAtom++) {
		mmdb::Atom* theAtom = SelAtom[iAtom];
		string atomName(theAtom->name);
		string residueName(theAtom->residue->name);
		double theAtomCharge;
		theAtomCharge = theChargeTable->getCharge(residueName, atomName);
		theAtom->charge = theAtomCharge;
	}
	return 0;
}

int CXXUtils::assignUnitedAtomRadius  (mmdb::Manager* theManager) {

   int allAtoms = theManager->NewSelection();
    const char* all = "/*/*/*/*.*";
    theManager->Select(allAtoms, mmdb::STYPE_ATOM, all, mmdb::SKEY_NEW);
    assignUnitedAtomRadius(theManager, allAtoms);
    theManager->DeleteSelection(allAtoms);
    return 0;
}

int CXXUtils::assignUnitedAtomRadius  (mmdb::Manager* theManager, int selHnd) {
	// Add a radius property to the atoms
	int iRadiusHandle = theManager->RegisterUDReal(mmdb::UDR_ATOM, "PerAtomRadius");
	if (!iRadiusHandle) {
		printf ( " registration failed.\n" );
		exit ( 1 );
	}
	if (iRadiusHandle==mmdb::UDDATA_WrongUDRType) {
		printf ( " wrong registration type used.\n" );
		exit ( 2 );
	}
	
	mmdb::Atom* *SelAtom;
	int nSelAtoms;
	theManager->GetSelIndex(selHnd, SelAtom, nSelAtoms);
	
    
    std::map<std::string,std::map<std::string, float> >mappedRadii;
    for (int iAtomRadius = 0; iAtomRadius < CXXUtils::nAtomRadii; iAtomRadius++){
        std::string residueName(CXXUtils::unitedAtomRadii[iAtomRadius].residueName);
        
        if (mappedRadii.find(residueName) == mappedRadii.end()) {
            std::map<std::string, float> atomMap;
            mappedRadii[residueName] = atomMap;
        }
        std::string atomName(CXXUtils::unitedAtomRadii[iAtomRadius].atomName);
        mappedRadii[residueName][atomName] = CXXUtils::unitedAtomRadii[iAtomRadius].radius;
    }
    
    for (int iAtom = 0; iAtom < nSelAtoms; iAtom++) {
        mmdb::Atom* anAtom = SelAtom[iAtom];
        std::string atomName(anAtom->name);
        std::string residueName("*  ");
        if (anAtom->residue != NULL){
            residueName = std::string(anAtom->residue->name);
        }
        
        std::map<std::string,std::map<std::string, float> >::iterator residueMapIter = mappedRadii.find(residueName);
        std::map<std::string,std::map<std::string, float> >::iterator generalMapIter = mappedRadii.find("*  ");
        
        double radius = 1.8;
        if (residueMapIter != mappedRadii.end() && residueMapIter->second.find(atomName) != residueMapIter->second.end()){
            radius = residueMapIter->second[atomName];
        }
        else if (generalMapIter->second.find(atomName) != generalMapIter->second.end()){
            radius = generalMapIter->second[atomName];
        }
        
		anAtom->PutUDData(iRadiusHandle, radius);
		//		 if (unassigned) cout << "Not able to match atom :" << atomName << ":" << residueName << ":\n";
	}
	return 0;
}

double CXXUtils::getAtomRadius(mmdb::Manager* theManager, mmdb::Atom* theAtom){
    //Here get handle of a radius data type from MMDB if such has been stored
    int iRadiusHandle = theManager->GetUDDHandle(mmdb::UDR_ATOM, "PerAtomRadius");
    double theRadius;
    if (iRadiusHandle>0){
        int success = theAtom->GetUDData (iRadiusHandle, theRadius);
        if (success != mmdb::UDDATA_Ok) theRadius = 1.8;
    }
    else theRadius = mmdb::getVdWaalsRadius(theAtom->element);
    return theRadius;
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

int CXXUtils::selectionStringToSelHnd(mmdb::Manager* allAtomsManager_in, std::string selectionString, int existingSelection, mmdb::SELECTION_KEY selKeyRequest){
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
	mmdb::Atom** SelAtoms;
	int nSelAtoms;
	allAtomsManager_in->GetSelIndex(selHnd, SelAtoms, nSelAtoms);
	cout << "Selection now contains " << nSelAtoms << " atoms\n";
	free (pstring);
	return selHnd;
}

int CXXUtils::unCharge(mmdb::Manager* theManager, int selHnd){
	mmdb::Atom* *SelAtom;
	int nSelAtoms;
	theManager->GetSelIndex(selHnd, SelAtom, nSelAtoms);
	
	//Assign atom charges
	for (int iAtom = 0; iAtom < nSelAtoms; iAtom++) {
		mmdb::Atom* theAtom = SelAtom[iAtom];
		theAtom->charge = 0.;
	}
	return 0;
}

