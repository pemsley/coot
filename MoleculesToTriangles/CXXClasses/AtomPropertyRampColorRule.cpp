/*
 * MoleculesToTriangles/CXXClasses/AtomPropertyRampColorRule.cpp
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


