/*
 *  AtomPropertyRampColorRule.mm
 *  Aesop
 *
 *  Created by Martin Noble on 27/02/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include <string>
#include "AtomPropertyRampColorRule.h"

AtomPropertyRampColorRule *AtomPropertyRampColorRule::defaultRampRule(){
	AtomPropertyRampColorRule
	*result = new AtomPropertyRampColorRule();
	return result;
}
/*
void AtomPropertyRampColorRule::prepareForSelectionInMMDB(int handle, mmdb::Manager *mmdb){
    if (rampType == ResidueNumber){
        //Count residues in selection, and set start and end of ramp accordingly
        mmdb::Atom** selectedAtoms;
        int nSelAtoms;
        mmdb->GetSelIndex(handle, selectedAtoms, nSelAtoms);
        if (nSelAtoms > 0){
            mmdb::Atom* firstAtom = selectedAtoms[0];
            mmdb::Atom* lastAtom = selectedAtoms[nSelAtoms-1];
            float stv = firstAtom->GetResidueNo();
            setStartValue(stv);
            float ltv = lastAtom->GetResidueNo();
            setEndValue(ltv);
        }
    }
}
 */


