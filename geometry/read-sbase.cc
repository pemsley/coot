/* geometry/protein-geometry.cc
 * 
 * Copyright 2010 The University of Oxford
 * Copyright 2013, 2014, 2015 by Medical Research Council
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */


#include <algorithm>
#include <iomanip>
#include <stdlib.h>
#include <string.h> // strlen, strcpy for interface to CCP4SRS.

#include "utils/coot-utils.hh"
#include "geometry/protein-geometry.hh"

#ifdef HAVE_CCP4SRS
#include "srs-interface.hh"
#endif 


// This is a toy/test function.
// 
void
coot::protein_geometry::read_ccp4srs_residues() {

#ifdef HAVE_CCP4SRS

   if (ccp4srs) { 
      ccp4srs::Monomer *monomer_p = NULL;
      std::vector<std::string> local_residue_codes;
      local_residue_codes.push_back("ASN");
      local_residue_codes.push_back("TRP");

      for (unsigned int i=0; i<local_residue_codes.size(); i++) {
	 monomer_p = ccp4srs->getMonomer(local_residue_codes[i].c_str());
	 if (monomer_p) { 
	    std::cout << monomer_p->ID() << " " << monomer_p->chem_name() << std::endl;
	    for (int iat=0; iat<monomer_p->n_atoms(); iat++) {
	       ccp4srs::Atom *at = monomer_p->atom(iat);
	       std::cout << "    " << at->name() << " ("
			 << at->x() << "," << at->y() << ","
			 << at->z() <<")\n";
	    } 
	 } else {
	    std::cout << "WARNING:: structure " << local_residue_codes[i]
		      << " not found in CCP4SRS " << std::endl;
	 }
	 delete monomer_p;
	 monomer_p = NULL;
      }
   } else {
      std::cout << "WARNING:: CCP4SRS not initialised"  << std::endl;
   }
#endif    
}

// Here res_name is the tlc/comp_id.
// 
mmdb::Residue *
coot::protein_geometry::get_ccp4srs_residue(const std::string &res_name) const {

   mmdb::Residue *residue_p = NULL;

#ifdef HAVE_CCP4SRS   
   if (ccp4srs) {
      ccp4srs::Monomer *monomer_p = ccp4srs->getMonomer(res_name.c_str());
      if (monomer_p) {
	 residue_p = new mmdb::Residue;
	 for (int iat=0; iat<monomer_p->n_atoms(); iat++) {
	    ccp4srs::Atom *at = monomer_p->atom(iat);

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
	       mmdb::Atom *new_atom = new mmdb::Atom;
	       new_atom->SetCoordinates(at->x(), at->y(), at->z(), 1.0, 30.0);
	       new_atom->SetAtomName(new_atom_name.c_str());
	       new_atom->SetElementName(at->element());
	       new_atom->Het = true; 
	       residue_p->AddAtom(new_atom);
	    } else {
	       std::cout << "WARNING:: rejecting " << res_name << " ccp4srs atom :" << new_atom_name
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
   if (ccp4srs) {
      unsigned int n_entries = ccp4srs->n_entries();
      for (unsigned int i=1; i<n_entries; i++) {
	 ccp4srs::Monomer *Monomer = ccp4srs->getMonomer(i, NULL);
	 std::string id = Monomer->ID();
	 std::string chem_name = Monomer->chem_name();
	 // std::cout << "i " << i <<  " monomer id  " << id << std::endl;
	 if (chem_name.find(compound_name_frag) != std::string::npos) {
	    std::pair<std::string, std::string> p(id, chem_name);
	    v.push_back(p);
	 }
      }
   } else {
      std::cout << "WARNING:: Null CCP4SRS" << std::endl;
   } 
#endif // HAVE_CCP4SRS   
   return v;
}


// return mmdb ccp4srs return codes
//
// Try to use the MONOMER_DIR_STR, ie. COOT_CCP4SRS_DIR first, if that
// fails then use the fallback directory ccp4srs_monomer_dir_in
// 
int
coot::protein_geometry::init_ccp4srs(const std::string &ccp4srs_monomer_dir_in) {

#ifdef HAVE_CCP4SRS

   // if srs-dir is not given, then the default should be $CCP4_MASTER/share/ccp4srs

   bool debug = false;
   int RC = ccp4srs::CCP4SRS_FileNotFound; // initial status.
   // std::cout << "init_ccp4srs() with " << ccp4srs_monomer_dir_in << std::endl;
   std::string dir;
   const char *d1 = getenv(MONOMER_DIR_STR); // "COOT_CCP4SRS_DIR"
   const char *d2 = getenv("CCP4");

   if (d1) {
      if (file_exists(d1))
	 dir = d1;
   } else {
      if (d2) {
	 std::string dir_a = util::append_dir_dir(d2, "share");
	 std::string dir_b = util::append_dir_dir(dir_a, "ccp4srs");
	 if (file_exists(dir_b))
	    dir = dir_b;
      } else {
	 if (file_exists(ccp4srs_monomer_dir_in))
	    dir = ccp4srs_monomer_dir_in;
      }
   }
   
   if (dir.length()) {
      std::cout << "INFO:: CCP4SRS::loadIndex from dir: " << dir << std::endl;
      ccp4srs = new ccp4srs::Manager;
      RC = ccp4srs->loadIndex(dir.c_str());
      if (debug)
	 std::cout << "init_ccprsrs() ... loadIndex() returned " << RC << std::endl;
      if (RC != ccp4srs::CCP4SRS_Ok) {
         std::cout << "CCP4SRS init problem." << std::endl;
	 delete ccp4srs;
	 ccp4srs = NULL;
      }
   } else {
      std::cout << "WARNING:: init_ccp4srs() no dir" << std::endl;
   } 
   return RC;
#else
   return -1;
#endif   
}


// used to created data from ccp4srs to put into protein_geometry
// object.
//
// return a flag to signify success.
//
// This relies on CCP4SRS being setup before we get to make this
// call.
// 
bool
coot::protein_geometry::fill_using_ccp4srs(const std::string &monomer_type) {

   bool success = false;
#ifdef HAVE_CCP4SRS
   dictionary_residue_restraints_t rest;
   success = rest.fill_using_ccp4srs(ccp4srs, monomer_type);
#endif
   return false;
}

#ifdef HAVE_CCP4SRS
bool
coot::dictionary_residue_restraints_t::fill_using_ccp4srs(ccp4srs::Manager *srs_manager,
							  const std::string &monomer_type) {
   
   bool success = false;
   residue_info.comp_id = monomer_type;

   if (! srs_manager) {
     std::cout << "WARNING:: fill_using_ccp4srs() Null CCP4SRS database " << std::endl;
   } else { 
      ccp4srs::Monomer *monomer_p = srs_manager->getMonomer(monomer_type.c_str());
      if (! monomer_p) {
	std::cout << "WARNING:: Null monomer ccp4srs::getMonomer()" << std::endl;
      } else { 

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
	 residue_info = d;
	 

	 // atoms
	 
	 for (int iat=0; iat<monomer_p->n_atoms(); iat++) {
	    ccp4srs::Atom *at = monomer_p->atom(iat);
	    std::pair<bool, float> pc(0,0); // partial charge.
	    std::string type_symbol = at->element();
	    std::string type_energy = at->energy_type();
	    std::string atom_id = at->name();
	    std::string atom_id_4c = atom_id_mmdb_expand(atom_id, type_symbol);
	    char rcsb_chirality = at->rcsb_chirality();
	    char ccp4_chirality = at->ccp4_chirality();
	    coot::dict_atom dict_at(atom_id,
				    atom_id_4c,
				    type_symbol,
				    type_energy,
				    pc);

	    // we need to add a sanity check on the position of the
	    // atoms (e.g. S in SO4 is at (-1.7976e+308 -1.7976e+308
	    // -1.7976e+308)).
            // 
	    bool add_atom_coords_ok = false;
	    if (fabs(at->x()) < 2000) { 
	       if (fabs(at->y()) < 2000) { 
		  if (fabs(at->z()) < 2000) {
		     add_atom_coords_ok = true;
		  }
               }
            }
            if (add_atom_coords_ok) { 
               clipper::Coord_orth pos(at->x(), at->y(), at->z());
               std::pair<bool, clipper::Coord_orth> p(true, pos);
               // write_cif only looks for real model coords! This needs fixing there. 
               // not this hack.
               // dict_at.add_pos(dict_atom::IDEAL_MODEL_POS, p);
               dict_at.add_pos(dict_atom::REAL_MODEL_POS, p);
            }

	    atom_info.push_back(dict_at);
	 }

	 if (atom_info.size() > 1) {

	    // bonds
	 
	    for (int ib=0; ib<monomer_p->n_bonds(); ib++) {
	       ccp4srs::Bond *bond = monomer_p->bond(ib);
	       int ind_1 = bond->atom1();
	       int ind_2 = bond->atom2();
	       int order = bond->order();
	       double dist = bond->length();
	       double esd = bond->length_esd();
	       std::string type = "single";
	       switch (order) {
	       case ccp4srs::Bond::noOrder:
		  type = "single";
		  break;
	       case ccp4srs::Bond::Single:
		  type = "single";
		  break;
	       case ccp4srs::Bond::Aromatic:
		  type = "aromatic";
		  break;
	       case ccp4srs::Bond::Double:
		  type = "double";
		  break;
	       case ccp4srs::Bond::Triple:
		  type = "triple";
		  break;
	       case ccp4srs::Bond::Deloc:
		  type = "deloc";
		  break;
	       case ccp4srs::Bond::Covalent:
		  type = "covalent";
		  break;
	       case ccp4srs::Bond::Metal:
		  type = "metal";
		  break;
	       default:
		  type = "single";
	       } 


	       std::string atom_name_1 = monomer_p->atom(ind_1)->name();
	       std::string atom_name_2 = monomer_p->atom(ind_2)->name();

	       if (false) // debug
		  std::cout << "SRS monomer bond: " << monomer_type 
                            << " atom index " 
                            << std::setw(3) << ind_1 << " " << std::setw(3) << ind_2 
                            << " atom names: " << atom_name_1 << " and " << atom_name_2 
			    << " order " << order << " type: " << type << std::endl;

	       coot::dict_bond_restraint_t dict_bond(atom_name_1, atom_name_2, type, dist, esd);
	       success = true;
	       bond_restraint.push_back(dict_bond);
	    }

	    // angles
	 
	    for (int ia=0; ia<monomer_p->n_angles(); ia++) {
	       ccp4srs::Angle *angle = monomer_p->angle(ia);
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
	       angle_restraint.push_back(dict_angle);
	    }

	    // torsions
	    for (int ia=0; ia<monomer_p->n_torsions(); ia++) {
	       ccp4srs::Torsion *torsion = monomer_p->torsion(ia);
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
	       torsion_restraint.push_back(dict_torsion);
	    }

	    // chirals
	 
	    for (int ic=0; ic<monomer_p->n_chicenters(); ic++) {
	       ccp4srs::ChiCenter *chiral = monomer_p->chicenter(ic);
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
	       chiral_restraint.push_back(dict_chiral);
	    }

	    // planes

	    for (int ip=0; ip<monomer_p->n_planes(); ip++) {
	       ccp4srs::Plane *plane = monomer_p->plane(ip);
	       std::string id = plane->id();
	       std::vector<std::string> atom_names;
	       std::vector<double> esds;
	       for (int iat=0; iat<plane->size(); iat++) { 
		  std::string atom_name = monomer_p->atom(plane->atom(iat))->name();
		  atom_names.push_back(atom_name);
	       }
	       double esd = plane->esds()[0];
	       coot::dict_plane_restraint_t dict_plane(id, atom_names, esd);
	       plane_restraint.push_back(dict_plane);
	    }
	 }
      }
   }
   if (success) {
      if (false) // debugging
	 std::cout << "INFO:: adding restraint from SRS " << atom_info.size() << " atoms and "
		   << bond_restraint.size() << " bonds, " 
		   << angle_restraint.size() << " angles, " 
		   << torsion_restraint.size() << " torsions," 
		   << chiral_restraint.size() << " chirals " 
		   << plane_restraint.size() << std::endl;
   }
   return success;
}
#endif // HAVE_CCP4SRS


// A new pdb file has been read in (say).  The residue types
// have been compared to the dictionary.  These (comp_ids) are
// the types that are not in the dictionary.  Try to load an
// ccp4srs description at least so that we can draw their bonds
// correctly.  Use fill_using_ccp4srs().
//
// the passed type should be a set - so that we don't have to uniquify here.
bool
coot::protein_geometry::try_load_ccp4srs_description(const std::vector<std::string> &comp_ids_with_duplicates) {

#ifdef HAVE_CCP4SRS
   bool status = false; // none added initially.

   std::vector<std::string> uniques;
   for (unsigned int ic=0; ic<comp_ids_with_duplicates.size(); ic++) { 
      if (std::find(uniques.begin(), uniques.end(), comp_ids_with_duplicates[ic]) ==
	  uniques.end())
	 uniques.push_back(comp_ids_with_duplicates[ic]);
   }

   if (ccp4srs) {
      for (unsigned int i=0; i<uniques.size(); i++) {
	 std::cout << i << " " << uniques[i] << std::endl;
	 const std::string &comp_id = uniques[i];
	 if (is_non_auto_load_ligand(comp_id)) {
	    std::cout << "INFO:: ccp4srs-descriptions: comp-id: " << comp_id
		      << " is marked for non-autoloading - ignore " << std::endl;
	 } else { 
	    bool s = fill_using_ccp4srs(comp_id);
	    std::cout << "DEBUG:: ccp4srs dictionary for " << comp_id << " successfully loaded "
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
coot::protein_geometry::residue_from_best_match(mmdb::math::Graph &graph1, mmdb::math::Graph &graph2,
						mmdb::math::GraphMatch &match, int n_match,
						ccp4srs::Monomer *monomer_p) const {

   match_results_t r("", "", NULL);
   int best_match = -1;
   for (int imatch=0; imatch<n_match; imatch++) {
      r.success = 1;
      r.name = monomer_p->chem_name();
      r.comp_id = graph2.GetName();
      int n;
      mmdb::realtype p1, p2;
      mmdb::ivector FV1, FV2;
      match.GetMatch(imatch, FV1, FV2, n, p1, p2); // n p1 p2 set
      if (0)
	 std::cout << "   match " << imatch << " " << " set n pairs " << n << std::endl;
      int n_type_match = 0;
      for (int ipair=1; ipair<=n; ipair++) {
	 mmdb::math::Vertex *V1 = graph1.GetVertex ( FV1[ipair] );
	 mmdb::math::Vertex *V2 = graph2.GetVertex ( FV2[ipair] );
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
int
coot::protein_geometry::ccp4_srs_n_entries() const {

   if (! ccp4srs) {
      std::cout << "WARNING:: CCP4SRS is not initialized" << std::endl;
      return 0;
   } else {
      return ccp4srs->n_entries();
   }
}
#endif // HAVE_CCP4SRS


#ifdef HAVE_CCP4SRS

// if srs_idx_start is negative then run through all the
// monomers
// else do the range (eg. 1001->2000 (idxs are inclusive)).
// 
std::vector<coot::match_results_t>
coot::protein_geometry::compare_vs_ccp4srs(mmdb::math::Graph *graph_1,
					   float similarity, int n_vertices,
					   int srs_idx_start, int srs_idx_end,
					   bool fill_match_graphs) const {

   std::vector<coot::match_results_t> v;

   int minMatch = int(similarity * n_vertices);
   std::cout << "INFO:: match.MatchGraphs must match at least "
	     << minMatch << " atoms from " << similarity << " * " << n_vertices << std::endl;
   
   if (! ccp4srs) {
      std::cout << "WARNING:: CCP4SRS is not initialized" << std::endl;
   } else {
      int l = ccp4srs->n_entries();
      std::cout << "INFO:: compare_vs_ccp4srs(): found " << l << " entries in CCP4 SRS" << std::endl;
      mmdb::math::Graph  *graph_2 = NULL;
      int rc = 0;
      mmdb::io::File *fp = NULL;

      if (srs_idx_start < 1) {
	 srs_idx_start = 1;
	 srs_idx_end = l;
      } else {
	 if (srs_idx_start >= l)
	    srs_idx_start = 1;
	 if (srs_idx_end >= (l-1))
	    srs_idx_end = l-1;
      }
      
      for (int i=srs_idx_start; i<=srs_idx_end; i++) {
	 ccp4srs::Monomer *Monomer = ccp4srs->getMonomer(i, fp);
	 if (! fp) // the first time it is null
	    fp = ccp4srs->getStructFile();

	 if (! Monomer)  {
	    std::cout << "Null monomer " << i << std::endl;
	 } else {
	    std::string id = Monomer->ID();
	    // std::cout << "monomer index i " << i <<  " monomer id  " << id << std::endl;
	    if (id.length()) {
	       graph_2 = Monomer->getGraph(&rc);
	       if (graph_2) {
		  graph_2->Build(true);
		  graph_2->MakeSymmetryRelief(true); // necessary if graph_1 has symmetry relief?

		  if (rc < 10000) {
		     mmdb::math::GraphMatch match;
		     match.SetTimeLimit(2); // seconds

		     mmdb::math::VERTEX_EXT_TYPE vertex_ext=mmdb::math::EXTTYPE_Equal; // mmdb default
		     bool vertext_type = true;
		     match.SetMaxNofMatches(1, true); // only need find 1 match
		     match.MatchGraphs(graph_1, graph_2, minMatch, vertext_type, vertex_ext);
		     int n_match = match.GetNofMatches();
		     if (n_match > 0) {
			std::cout << "INFO:: " << id
				  << " match NumberofMatches (similar graphs): " << n_match << std::endl;

			bool really_match = false;
			std::vector<std::pair<int, int> > graph_match_atom_indices;
			for (int imatch=0; imatch<n_match; imatch++) {
			   int n;
			   mmdb::realtype p1, p2;
			   mmdb::ivector FV1, FV2;
			   match.GetMatch(imatch, FV1, FV2, n, p1, p2); // n p1 p2 set
			   std::cout << "   match " << imatch << " matched " << n << " atoms"
				     << std::endl;
			   if (n >= minMatch) really_match = true;
			   for (int ipair=1; ipair<=n; ipair++) {
			      mmdb::math::Vertex *V1 = graph_1->GetVertex(FV1[ipair]);
			      mmdb::math::Vertex *V2 = graph_2->GetVertex(FV2[ipair]);
			      if ((!V1) || (!V2))  {
				 std::cout << "Can't get vertices for match " << ipair << std::endl;
			      } else  {
				 std::cout << "   " << V1->GetUserID() << " " << V2->GetUserID()
					   << std::endl;
				 std::pair<int, int> p(V1->GetUserID(), V2->GetUserID());
				 graph_match_atom_indices.push_back(p);
			      }
			   }
			   if (really_match)
			      break;
			}

			if (really_match) {
			   mmdb::Residue *residue_p = NULL; // for now
			   std::string name = Monomer->chem_name();
			   match_results_t mr(id, name, residue_p);
			   mr.graph_match_atom_indices = graph_match_atom_indices;
			   v.push_back(mr);
			}
		     }
		  }
	       }
	       delete graph_2; // fixes memory leak?
	       graph_2 = NULL;
	    }
	 }
	 delete Monomer;
	 Monomer = NULL;
      }
   }
   return v;
}
#endif // HAVE_CCP4SRS

#ifdef HAVE_CCP4SRS
std::vector<std::string>
coot::protein_geometry::get_available_ligand_comp_id(const std::string &hoped_for_head,
						     unsigned int n_top) const {

   std::string r;
   std::vector<std::string> available_comp_ids;
   
   if (! ccp4srs) {
      std::cout << "WARNING:: CCP4SRS is not initialized" << std::endl;
   } else {
      std::string next_letter = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
      std::string::size_type l_head = hoped_for_head.length();
      std::vector<std::string> matching_codes;
      int l = ccp4srs->n_entries();

      // find comp_ids that match the letters of hoped_for_head
      // 
      for (int i=0; i<l ;i++)  {
	 ccp4srs::Monomer *Monomer = ccp4srs->getMonomer(i, NULL);
	 std::string id = Monomer->ID();
	 std::string::size_type l_monomer = id.length();
	 if (l_monomer >= l_head) {

	    bool matched = true;
	    for (std::string::size_type ii=0; ii<l_head; ii++) { 
	       if (hoped_for_head[ii] != id[ii]) {
		  matched = false;
		  break;
	       }
	    }
	    if (matched) {
	       matching_codes.push_back(id);
	    }
	 } 
      }

      // debug::
      if (false) {
	 std::cout << "dictionary matching codes:" << std::endl;
	 for (unsigned int i=0; i<matching_codes.size(); i++) {
	    if (i>0)
	       if (i%10 == 0) 
		  std::cout << std::endl;
	    std::cout << "    " << matching_codes[i];
	 }
	 std::cout << std::endl;
	 std::cout << "-----" << std::endl;
      }

      if (hoped_for_head.length() == 3) {
	 // simply check if hoped_for_head is in matching_codes.  If
	 // it is then the user can't have it.
	 std::vector<std::string>::const_iterator it;
	 it = std::find(matching_codes.begin(),
			matching_codes.end(),
			hoped_for_head);
	 if (it == matching_codes.end())
	    available_comp_ids.push_back(hoped_for_head);
      }


      //
      if (hoped_for_head.length() == 2) {
	 std::vector<std::string>::const_iterator it;
	 int nli = 0;
	 for (std::string::size_type idx=0; idx<next_letter.length(); idx++) { 
	    std::string test_code = hoped_for_head;
	    test_code += next_letter[idx];
	    it = std::find(matching_codes.begin(),
			   matching_codes.end(),
			   test_code);
	    if (it == matching_codes.end()) {
	       available_comp_ids.push_back(test_code);
	       if (available_comp_ids.size() >= n_top)
		  break;
	    } 
	 }
      }

      if (hoped_for_head.length() == 1) {
	 std::vector<std::string>::const_iterator it;
	 unsigned int nli = 0;
	 bool done = false;
	 for (std::string::size_type idx_1=0; idx_1<next_letter.length(); idx_1++) { 
	    for (std::string::size_type idx_2=0; idx_2<next_letter.length(); idx_2++) { 
	       std::string test_code = hoped_for_head;
	       test_code += next_letter[idx_1];
	       test_code += next_letter[idx_2];
	       it = std::find(matching_codes.begin(),
			      matching_codes.end(),
			      test_code);
	       if (it == matching_codes.end()) {
		  available_comp_ids.push_back(test_code);
		  if (available_comp_ids.size() >= n_top)
		     done = true;
	       }
	       if (done)
		  break;
	    }
	    if (done)
	       break;
	 }
      }

      if (hoped_for_head.length() == 0) {
	 std::vector<std::string>::const_iterator it;
	 unsigned int nli = 0;
	 bool done = false;
	 for (std::string::size_type idx_1=0; idx_1<next_letter.length(); idx_1++) { 
	    for (std::string::size_type idx_2=0; idx_2<next_letter.length(); idx_2++) { 
	       for (std::string::size_type idx_3=0; idx_3<next_letter.length(); idx_3++) { 
		  std::string test_code;
		  test_code += next_letter[idx_1];
		  test_code += next_letter[idx_2];
		  test_code += next_letter[idx_3];
		  it = std::find(matching_codes.begin(),
				 matching_codes.end(),
				 test_code);
		  if (it == matching_codes.end()) {
		     available_comp_ids.push_back(test_code);
		     if (available_comp_ids.size() >= n_top)
			done = true;
		  }
		  if (done)
		     break;
	       }
	       if (done)
		  break;
	    }
	    if (done)
	       break;
	 }
      }
   }
   return available_comp_ids;
}
#endif  // HAVE_CCP4SRS



#ifdef HAVE_CCP4SRS
// return 0 on strangeness, to pass in search.
// 
int coot::get_min_match(const int &n1, const float similarity) {
   int most_1 = int (similarity * float(n1));
   return most_1;
   // return n1;
}
#endif
