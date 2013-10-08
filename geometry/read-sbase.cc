/* geometry/protein-geometry.cc
 * 
 * Copyright 2010 The University of Oxford
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


#include <algorithm>
#include <stdlib.h>

#include "utils/coot-utils.hh"
#include "geometry/protein-geometry.hh"

#ifdef HAVE_CCP4SRS
#include "ccp4srs/ccp4srs_defs.h"
#include "srs-interface.hh"
#endif 

#include <string.h> // strlen, strcpy for interface to sbase.

void
coot::protein_geometry::read_ccp4srs_residues() {

#ifdef HAVE_CCP4SRS

   if (SBase) { 
      CCP4SRSMonomer *monomer_p = NULL;
      std::vector<std::string> local_residue_codes;
      local_residue_codes.push_back("ASN");
      local_residue_codes.push_back("TRP");

      for (int i=0; i<local_residue_codes.size(); i++) {
	 monomer_p = SBase->getMonomer(local_residue_codes[i].c_str());
	 if (monomer_p) { 
	    std::cout << monomer_p->ID() << " " << monomer_p->chem_name() << std::endl;
	    for (int iat=0; iat<monomer_p->n_atoms(); iat++) {
	       CCP4SRSAtom *at = monomer_p->atom(iat);
	       std::cout << "    " << at->name() << " ("
			 << at->x() << "," << at->y() << ","
			 << at->z() <<")\n";
	    } 
	 } else {
	    std::cout << "WARNING:: structure " << local_residue_codes[i]
		      << " not found in SBase " << std::endl;
	 } 
      }
   } else {
      std::cout << "WARNING:: CCP4SRS not initialised"  << std::endl;
   }
#endif    
}

// Here res_name is the tlc/comp_id.
// 
CResidue *
coot::protein_geometry::get_ccp4srs_residue(const std::string &res_name) const {

   CResidue *residue_p = NULL;

#ifdef HAVE_CCP4SRS   
   if (SBase) {
      CCP4SRSMonomer *monomer_p = SBase->getMonomer(res_name.c_str());
      if (monomer_p) {
	 residue_p = new CResidue;
	 for (int iat=0; iat<monomer_p->n_atoms(); iat++) {
	    CCP4SRSAtom *at = monomer_p->atom(iat);

	    std::string new_atom_name = coot::atom_id_mmdb_expand(at->name(), at->element());

	    // we need to add a sanity check on the position of the
	    // atoms (e.g. S in SO4 is at (-1.7976e+308 -1.7976e+308
	    // -1.7976e+308)).

	    bool add_atom_ok = 0;
	    if (fabs(at->x()) < 2000) { 
	       if (fabs(at->y()) < 2000) { 
		  if (fabs(at->z()) < 2000) {
		     add_atom_ok = 1;
		  }
	       }
	    }

	    if (add_atom_ok) { 
	       CAtom *new_atom = new CAtom;
	       new_atom->SetCoordinates(at->x(), at->y(), at->z(), 1.0, 30.0);
	       new_atom->SetAtomName(new_atom_name.c_str());
	       new_atom->SetElementName(at->element());
	       new_atom->Het = true; 
	       residue_p->AddAtom(new_atom);
	    } else {
	       std::cout << "WARNING:: rejecting " << res_name << " SBase atom :" << new_atom_name
			 << ": at (" << at->x() << " " << at->y() << " " << at->z() << ")" << std::endl;
	    }
	 }
      }
   }

#endif // HAVE_CCP4SRS   
   return residue_p;
}


std::vector<std::pair<std::string, std::string> >
coot::protein_geometry::matching_ccp4srs_residues_names(const std::string &compound_name_frag) const {

   std::vector<std::pair<std::string,std::string> > v;
#ifdef HAVE_CCP4SRS   
   std::string compound_name = coot::util::upcase(compound_name_frag);
   if (SBase) {
      int n_entries = SBase->n_entries();
      // std::pair<std::string,std::string> p(monomer_p->ID(), monomer_p->name());
   }
#endif // HAVE_CCP4SRS   
   return v;
}


// return mmdb sbase return codes
//
// Try to use the MONOMER_DIR_STR, ie. COOT_SBASE_DIR first, if that
// fails then use the fallback directory sbase_monomer_dir_in
// 
int
coot::protein_geometry::init_ccp4srs(const std::string &ccp4srs_monomer_dir_in) {

#ifdef HAVE_CCP4SRS

   int RC = CCP4SRS_FileNotFound; // initial status.
   // std::cout << "init_ccp4srs() with " << ccp4srs_monomer_dir_in << std::endl;
   std::string dir;
   const char *d1 = getenv(MONOMER_DIR_STR);
   const char *d2 = getenv("CCP4_LIB");

   if (d1) {
      if (file_exists(d1))
	 dir = d1;
   } else {
      if (d2) {
	 std::string dir_a = util::append_dir_dir(d2, "ccp4srs");
	 std::string dir_b = util::append_dir_dir(dir_a, "srsdata");
	 if (file_exists(dir_b))
	    dir = dir_b;
      } else {
	 if (file_exists(ccp4srs_monomer_dir_in))
	    dir = ccp4srs_monomer_dir_in;
      }
   }
   
   if (dir.length()) {
      std::cout << "about to loadIndex with dir " << dir << std::endl;
      SBase = new CCP4SRSBase;
      RC = SBase->loadIndex(dir.c_str());
      std::cout << "... loadIndex() returned " << RC << std::endl;
      if (RC != CCP4SRS_Ok) {
         std::cout << "CCP4SRS init problem." << std::endl;
	 delete SBase;
	 SBase = NULL;
      }
   } 
   return RC;
#else
   return -1;
#endif   
}


// used to created data from sbase to put into protein_geometry
// object.
//
// return a flag to signify success.
//
// This relies on SBase being setup before we get to make this
// call.
// 
bool
coot::protein_geometry::fill_using_ccp4srs(const std::string &monomer_type) {

#ifdef HAVE_CCP4SRS

   bool success = false;
   coot::dictionary_residue_restraints_t rest(true); // constructor for SBase
   rest.residue_info.comp_id = monomer_type;
   
   if (SBase) {
      CCP4SRSMonomer *monomer_p = SBase->getMonomer(monomer_type.c_str());
      if (monomer_p) {

	 // molecule info
	 std::string comp_id = monomer_type;
	 std::string three_letter_code = comp_id;
	 std::string name  = util::downcase(monomer_p->chem_name());
	 std::string group = util::downcase(monomer_p->chem_type());
	 int number_atoms_all = monomer_p->n_atoms();
	 int number_atoms_nh = 0;
	 std::string description_level = ".";

	 // count non-hydrogens
	 for (int iat=0; iat<monomer_p->n_atoms(); iat++)
	    if (strcmp(monomer_p->atom(iat)->element(), "H"))
	       number_atoms_nh++;


	 dict_chem_comp_t d(comp_id,
			    three_letter_code,
			    name,
			    group,
			    number_atoms_all,
			    number_atoms_nh,
			    description_level);
	 rest.residue_info = d;
      
	 

	 // atoms
	 
	 for (int iat=0; iat<monomer_p->n_atoms(); iat++) {
	    CCP4SRSAtom *at = monomer_p->atom(iat);
	    std::pair<bool, float> pc(0,0); // partial charge.
	    std::string type_symbol = at->element();
	    std::string type_energy = at->energy_type();
	    std::string atom_id = at->name();
	    std::string atom_id_4c = atom_id_mmdb_expand(atom_id, type_symbol);
	    coot::dict_atom dict_at(atom_id,
				    atom_id_4c,
				    type_symbol,
				    type_energy,
				    pc);
	    // std::cout << "pushing back an atom with name :" << atom_id_4c << ":" << std::endl;
	    rest.atom_info.push_back(dict_at);
	 }
      }

      if (rest.atom_info.size() > 1) {

	 // bonds
	 
	 for (int ib=0; ib<monomer_p->n_bonds(); ib++) {
	    CCP4SRSBond *bond = monomer_p->bond(ib);
	    int ind_1 = bond->atom1();
	    int ind_2 = bond->atom2();
	    int order = bond->order();
	    double dist = bond->length();
	    double esd = bond->length_esd();
	    if (1)
	       std::cout << "... "<< monomer_type << " atom index " << ind_1 << " " << ind_2
			 << " order " << order << std::endl;
	    std::string type = "single";
	    switch (order) {
	    case CCP4SRSBond::noOrder:
	       type = "single";
	       break;
	    case CCP4SRSBond::Single:
	       type = "single";
	       break;
	    case CCP4SRSBond::Aromatic:
	       type = "aromatic";
	       break;
	    case CCP4SRSBond::Double:
	       type = "double";
	       break;
	    case CCP4SRSBond::Triple:
	       type = "triple";
	       break;
	    case CCP4SRSBond::Deloc:
	       type = "deloc";
	       break;
	    case CCP4SRSBond::Covalent:
	       type = "covalent";
	       break;
	    case CCP4SRSBond::Metal:
	       type = "metal";
	       break;
	    default:
	       type = "single";
	    } 
	    std::string atom_name_1 = monomer_p->atom(ind_1)->name();
	    std::string atom_name_2 = monomer_p->atom(ind_2)->name();
	    std::cout << "adding bond between " << atom_name_1 << " and " << atom_name_2 << " type " << type << std::endl;
	    coot::dict_bond_restraint_t dict_bond(atom_name_1, atom_name_2, type, dist, esd);
	    success = true;
	    rest.bond_restraint.push_back(dict_bond);
	 }

	 // angles
	 
	 for (int ia=0; ia<monomer_p->n_angles(); ia++) {
	    CCP4SRSAngle *angle = monomer_p->angle(ia);
	    int ind_1 = angle->atom1();
	    int ind_2 = angle->atom2();
	    int ind_3 = angle->atom3();
	    double value = angle->value();
	    double esd = angle->esd();
	    std::string atom_name_1 = monomer_p->atom(ind_1)->name();
	    std::string atom_name_2 = monomer_p->atom(ind_2)->name();
	    std::string atom_name_3 = monomer_p->atom(ind_3)->name();
	    coot::dict_angle_restraint_t dict_angle(atom_name_1, atom_name_2, atom_name_3,
						    value, esd);
	    rest.angle_restraint.push_back(dict_angle);
	 }

	 // torsions
	 for (int ia=0; ia<monomer_p->n_torsions(); ia++) {
	    CCP4SRSTorsion *torsion = monomer_p->torsion(ia);
	    int ind_1 = torsion->atom1();
	    int ind_2 = torsion->atom2();
	    int ind_3 = torsion->atom3();
	    int ind_4 = torsion->atom4();
	    double value =   torsion->value();
	    double esd =     torsion->esd();
	    double period =  torsion->esd();
	    std::string id = torsion->id();
	    std::string atom_name_1 = monomer_p->atom(ind_1)->name();
	    std::string atom_name_2 = monomer_p->atom(ind_2)->name();
	    std::string atom_name_3 = monomer_p->atom(ind_3)->name();
	    std::string atom_name_4 = monomer_p->atom(ind_4)->name();
	    coot::dict_torsion_restraint_t dict_torsion(id,
							atom_name_1, atom_name_2, atom_name_3, atom_name_4,
							value, esd, period);
	    rest.torsion_restraint.push_back(dict_torsion);
	 }

	 // chirals
	 
	 for (int ic=0; ic<monomer_p->n_chicenters(); ic++) {
	    CCP4SRSChiCenter *chiral = monomer_p->chicenter(ic);
	    int ind_c = chiral->center();
	    int ind_1 = chiral->atom1();
	    int ind_2 = chiral->atom2();
	    int ind_3 = chiral->atom3();
	    int sign  = chiral->sign();
	    std::string id = chiral->id();
	    std::string atom_name_c = monomer_p->atom(ind_c)->name();
	    std::string atom_name_1 = monomer_p->atom(ind_1)->name();
	    std::string atom_name_2 = monomer_p->atom(ind_2)->name();
	    std::string atom_name_3 = monomer_p->atom(ind_3)->name();
	    coot::dict_chiral_restraint_t dict_chiral(id,
						      atom_name_c, atom_name_1, atom_name_2, atom_name_3, 
						      sign);
	    rest.chiral_restraint.push_back(dict_chiral);
	 }

	 // planes

	 for (int ip=0; ip<monomer_p->n_planes(); ip++) {
	    CCP4SRSPlane *plane = monomer_p->plane(ip);
	    std::string id = plane->id();
	    std::vector<std::string> atom_names;
	    std::vector<double> esds;
	    for (unsigned int iat=0; iat<plane->size(); iat++) {
	       std::string atom_name = monomer_p->atom(plane->atom(iat))->name();
	       atom_names.push_back(atom_name);
	    }
	    double esd = plane->esds()[0];
	    coot::dict_plane_restraint_t dict_plane(id, atom_names, esd);
	    rest.plane_restraint.push_back(dict_plane);
	 }
      }
   }
   if (success) {
      std::cout << "adding restraint " << rest.atom_info.size() << " atoms and "
		<< rest.bond_restraint.size() << " " << rest.angle_restraint.size()
		<< " " << rest.torsion_restraint.size() << " " << rest.chiral_restraint.size()
		<< " " << rest.plane_restraint.size() << std::endl;
      add(rest);
   }

   return success;
#else
   return false;
#endif      
}


// A new pdb file has been read in (say).  The residue types
// have been compared to the dictionary.  These (comp_ids) are
// the types that are not in the dictionary.  Try to load an
// sbase description at least so that we can draw their bonds
// correctly.  Use fill_using_sbase().
// 
bool
coot::protein_geometry::try_load_sbase_description(const std::vector<std::string> &comp_ids_with_duplicates) {

#ifdef HAVE_CCP4SRS   
   bool status = false; // none added initially.

   std::vector<std::string> uniques;
   for (unsigned int ic=0; ic<comp_ids_with_duplicates.size(); ic++) { 
      if (std::find(uniques.begin(), uniques.end(), comp_ids_with_duplicates[ic]) ==
	  uniques.end())
	 uniques.push_back(comp_ids_with_duplicates[ic]);
   }

   
   if (SBase) {
      for (unsigned int i=0; i<uniques.size(); i++) {
	 std::cout << i << " " << uniques[i] << std::endl;
	 const std::string &comp_id = uniques[i];
	 if (is_non_auto_load_ligand(comp_id)) {
	    std::cout << "INFO:: sbase-descriptions: comp-id: " << comp_id
		      << " is marked for non-autoloading - ignore " << std::endl;
	 } else { 
	    bool s = fill_using_ccp4srs(comp_id);
	    std::cout << "DEBUG:: sbase dictionary for " << comp_id << " successfully loaded "
		      << std::endl;
	    if (s)
	       status = true;
	 }
      }
   }
   return status; // return true of something was added.
#else
   return false;
#endif // HAVE_CCP4SRS   
}

#ifdef HAVE_CCP4SRS
coot::match_results_t
coot::protein_geometry::residue_from_best_match(CGraph &graph1, CGraph &graph2,
						CGraphMatch &match, int n_match,
						CCP4SRSMonomer *monomer_p) const {

   match_results_t r("", "", NULL);
   int best_match = -1;
   for (int imatch=0; imatch<n_match; imatch++) {
      r.success = 1;
      r.name = monomer_p->chem_name();
      r.comp_id = graph2.GetName();
      int n;
      realtype p1, p2;
      ivector FV1, FV2;
      match.GetMatch(imatch, FV1, FV2, n, p1, p2); // n p1 p2 set
      if (0)
	 std::cout << "   match " << imatch << " " << " set n pairs " << n << std::endl;
      int n_type_match = 0;
      for (int ipair=1; ipair<=n; ipair++) {
	 PCVertex V1 = graph1.GetVertex ( FV1[ipair] );
	 PCVertex V2 = graph2.GetVertex ( FV2[ipair] );
	 if ((!V1) || (!V2))  {
	    std::cout << "Can't get vertices for match " << ipair << std::endl;
	 } else {
	    int type_1 = V1->GetType();
	    int type_2 = V2->GetType();
	    if (type_1 == type_2) {
	       // std::cout << "   type match on " << type_1 << std::endl;
	       n_type_match++;
	    }
	 }
      }
      if (0)
	 std::cout << "This match matches " << n_type_match << " types out of " << n << std::endl;
   }
   return r;
}

#endif // HAVE_CCP4SRS


#ifdef HAVE_CCP4SRS
std::vector<coot::match_results_t>
coot::protein_geometry::compare_vs_ccp4srs(CGraph *graph_1, float similarity, int n_vertices) const {

   std::vector<coot::match_results_t> v;

   if (! SBase) {
      std::cout << "CCP4SRS is not initialized" << std::endl;
   } else {
      //  4.  Run the query through all databsae

      //  There are several methods for retrieving graphs
      //  from the sbase, here we use one most convenient
      //  for serial extractions.
      PCFile graphFile = SBase->getGraphFile();
      if (!graphFile)  {
	 printf ( "\n CCP4SRS graph file not found.\n" );
      }
  
      int exclude_H_flag = 1;  // neglect hydrogens
      CGraph *graph_2 = NULL;
      int min_match = coot::get_min_match(n_vertices, similarity);

      if (0) // debug
	 std::cout << "min_match " << min_match
		   << " n_vertices: " << n_vertices << " "
		   <<     (similarity * float(n_vertices)) << " " 
		   << int (similarity * float(n_vertices)) << " " 
		   << std::endl;

      std::cout << "INFO:: Close fragments must match " << min_match << " atoms of "
		<< n_vertices << std::endl;
      

      int nStructures = SBase->n_entries();
      int n_match = 0;
      std::cout << "INFO:: Searching " << nStructures << " CCP4SRS structures\n";
      for (int is=0; is<nStructures; is++)  {
	 int rc = SBase->getGraph(graphFile, graph_2, exclude_H_flag);
	 if (graph_2 == NULL) {
	    std::cout << "bad status on get graph " << is << std::endl;
	 } else {

	    int n2 = graph_2->GetNofVertices();
	    if ((n2 >= int (double(similarity) * double(n_vertices))) &&
		(n2 < (2.0 - double(similarity)) * double(n_vertices))) { 

	       graph_2->MakeVertexIDs();

	       // 20131008: We can't make the arg True, because that
	       // removes ~90% of the "hits" and returns molecules
	       // with the wrong bond orders.
	       // 
	       graph_2->Build(False); // 20100608 was True

	       CGraphMatch *match  = new CGraphMatch();
	       if (min_match > 0) { 

		  match->MatchGraphs(graph_1, graph_2, min_match);
		  int nMatches = match->GetNofMatches();
		  if (nMatches > 0) {
		     if (1)
			std::cout << "found " << nMatches << " match(es) for query in structure "
				  << is << " " << graph_1->GetName() << " vs " << graph_2->GetName()
				  << std::endl;
		     CFile *sf = SBase->getStructFile();
		     CCP4SRSMonomer *monomer_p = SBase->getMonomer(is, sf);
		     if (monomer_p) {
			if (monomer_p->chem_name()) {
			   std::cout << "    " << n_match << " " << graph_2->GetName() << " : "
				     << monomer_p->chem_name() << "\n";
			   coot::match_results_t res =
			      residue_from_best_match(*graph_1, *graph_2, *match, nMatches, monomer_p);
			   v.push_back(res);
			   n_match++;
			}
		     }
		  }
	       }
	       delete match;
	    }
	 }
      }
      std::cout << "Search complete" << std::endl;
  
      graphFile->shut();
      delete graphFile;

      delete graph_2;
   }
   return v;
}
#endif // HAVE_CCP4SRS   


#ifdef HAVE_CCP4SRS
// return 0 on strangeness, to pass in search.
// 
int coot::get_min_match(const int &n1, const float similarity) {
   int most_1 = int (similarity * float(n1));
   return most_1;
   // return n1;
}
#endif
