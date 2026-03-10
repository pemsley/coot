/* ideal/torsion-bonds.cc
 * 
 * Copyright 2015 by Medical Research Council
 * Author: Paul Emsley
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
#include <cstdlib>

#ifdef USE_BACKWARD
#include <utils/backward.hpp>
#endif

#include "simple-restraint.hh"
#include "coot-utils/coot-map-heavy.hh"

#include "coot-utils/coot-coord-extras.hh"
#include "coot-utils/contact-info.hh"
#include "coot-utils/atom-tree.hh"

#include "torsion-bonds.hh"

#include "utils/logging.hh"
extern logging logger;


// this can throw an exception
// 
std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> >
coot::torsionable_bonds(int imol, mmdb::Manager *mol, mmdb::PPAtom atom_selection,
			int n_selected_atoms,
			protein_geometry *geom_p) {

   std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > v;
   bool include_pyranose_ring_torsions_flag = false;

   std::vector<mmdb::Residue *> residues;
   std::map<mmdb::Residue *, std::vector<int> > atoms_in_residue;
   // fill residues and atoms_in_residue
   for (int i=0; i<n_selected_atoms; i++) {
      mmdb::Residue *r = atom_selection[i]->residue;
      if (std::find(residues.begin(), residues.end(), r) == residues.end())
	 residues.push_back(r);
      atoms_in_residue[r].push_back(i);
   }
   std::map<mmdb::Residue *, dictionary_residue_restraints_t> res_restraints;
   for (unsigned int ires=0; ires<residues.size(); ires++) {
      std::string rn = residues[ires]->GetResName();
      std::pair<bool, dictionary_residue_restraints_t> rest =
	 geom_p->get_monomer_restraints(rn, imol);
      if (! rest.first) {
	 std::string m = "Restraints not found for type ";
	 m += rn;
	 throw std::runtime_error(m);
      }
      res_restraints[residues[ires]] = rest.second;
   }

   for (unsigned int ires=0; ires<residues.size(); ires++) {
      // a coot-coord-extras function
      std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > v_inner =
	 torsionable_bonds_monomer_internal(residues[ires],
					    atom_selection,
					    n_selected_atoms,
					    include_pyranose_ring_torsions_flag,
					    geom_p);
      // std::cout << "found " << v_inner.size() << " monomer internal torsions for "
      // << residues[ires]->GetResName() << std::endl;
      for (unsigned int ip=0; ip<v_inner.size(); ip++)
	 v.push_back(v_inner[ip]);
   }

   std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > v_link =
      coot::torsionable_link_bonds(residues, mol, geom_p);
   for (unsigned int il=0; il<v_link.size(); il++)
      v.push_back(v_link[il]);

   if (0) // debug
      for (unsigned int ipair=0; ipair<v.size(); ipair++) {
	 std::cout << "   torsionable bond: "
		   << atom_spec_t(v[ipair].first) << "  "
		   << atom_spec_t(v[ipair].second)
		   << std::endl;
      }

   return v;
}

std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> >
coot::torsionable_link_bonds(std::vector<mmdb::Residue *> residues_in,
			     mmdb::Manager *mol, protein_geometry *geom_p) {

   std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > v;

   if (! mol)
      return v;

   std::vector<std::pair<bool, mmdb::Residue *> > residues(residues_in.size());
   for (unsigned int i=0; i<residues_in.size(); i++)
      residues[i] = std::pair<bool, mmdb::Residue *> (0, residues_in[i]);

   std::vector<atom_spec_t> dummy_fixed_atom_specs;
   std::vector<mmdb::Link> links;
   clipper::Xmap<float> dummy_xmap;
   coot::restraints_container_t restraints(residues, links, *geom_p, mol, dummy_fixed_atom_specs, &dummy_xmap);
   bonded_pair_container_t bpc = restraints.bonded_residues_from_res_vec(*geom_p);

   // add in the torsion: CB-CG-ND2-C1 (Psi-N)
   // add in the torsion: CG-ND2-C1-O5 (Phi-N)
   // t.geom.link_add_torsion("NAG-ASN", 1, 2, 2, 2, " C1 ", "ND2 ", " CG ", " CB ", 180, 40, 3, "Psi-N");
   // t.geom.link_add_torsion("NAG-ASN", 1, 1, 2, 2, " O5 ", " C1 ", "ND2 ", " CG ", 180, 40, 3, "Psi-N");
      
   // std::cout << "found LINKR linked pairs:\n   " <<  bpc;

   for (unsigned int i=0; i<bpc.bonded_residues.size(); i++) { 
      dictionary_residue_link_restraints_t link = geom_p->link(bpc[i].link_type);
      if (link.link_id != "") {
	 if (0)
	    std::cout << "   dictionary link found " << link.link_id << " with "
		      << link.link_bond_restraint.size() << " bond restraints and "
		      << link.link_torsion_restraint.size() << " link torsions " << std::endl;
	 
	 for (unsigned int ib=0; ib<link.link_bond_restraint.size(); ib++) { 
	    // we need to get the atoms and add them to "pairs".

	    if (0) 
	       std::cout << "   "
			 << link.link_bond_restraint[ib].atom_id_1_4c() << " "
			 << link.link_bond_restraint[ib].atom_1_comp_id << " to "
			 << link.link_bond_restraint[ib].atom_id_2_4c() << " " 
			 << link.link_bond_restraint[ib].atom_2_comp_id << " "
			 << " of "
			 << bpc[i].res_1->GetResName() << " to "
			 << bpc[i].res_2->GetResName()
			 << std::endl;
	    
	    mmdb::Atom *link_atom_1 = bpc[i].res_1->GetAtom(link.link_bond_restraint[ib].atom_id_1_4c().c_str());
	    mmdb::Atom *link_atom_2 = bpc[i].res_2->GetAtom(link.link_bond_restraint[ib].atom_id_2_4c().c_str());
	    if (link_atom_1 && link_atom_2) { 
	       std::pair<mmdb::Atom *, mmdb::Atom *> pair(link_atom_1, link_atom_2);
	       v.push_back(pair);
	    }
	 }


	 // Add in link torsion atoms if they were not already added
	 // above (because they were link bond restraints)
	 // 
	 for (unsigned int it=0; it<link.link_torsion_restraint.size(); it++) {
	    mmdb::Residue *res_for_at_2 = bpc[i].res_1;
	    mmdb::Residue *res_for_at_3 = bpc[i].res_1;
	    if (link.link_torsion_restraint[it].atom_2_comp_id == 1) res_for_at_2 = bpc[i].res_1;
	    if (link.link_torsion_restraint[it].atom_2_comp_id == 2) res_for_at_2 = bpc[i].res_2;
	    if (link.link_torsion_restraint[it].atom_3_comp_id == 1) res_for_at_3 = bpc[i].res_1;
	    if (link.link_torsion_restraint[it].atom_3_comp_id == 2) res_for_at_3 = bpc[i].res_2;

	    if (res_for_at_2 && res_for_at_3) {
	       mmdb::Atom *link_atom_1 = res_for_at_2->GetAtom(link.link_torsion_restraint[it].atom_id_2_4c().c_str());
	       mmdb::Atom *link_atom_2 = res_for_at_3->GetAtom(link.link_torsion_restraint[it].atom_id_3_4c().c_str());

	       if (link_atom_1 && link_atom_2) {
		  std::pair<mmdb::Atom *, mmdb::Atom *> pair(link_atom_1, link_atom_2);
		  if (std::find(v.begin(), v.end(), pair) == v.end())
		     v.push_back(pair);
	       }
	    }
	 }
      }
   }


   if (0) {
      std::cout << "---------------- torsionable_link_bonds() returns: " << std::endl;
      for (unsigned int i=0; i<v.size(); i++) { 
	 std::cout << "    " << i << " " << atom_spec_t(v[i].first) << " - "
		   << atom_spec_t(v[i].second) << std::endl;
      }
   }
   
   return v;
}


// And the atom_quad version of the above 2 functions (for setting link torsions)
// 
// this can throw an exception
// 
std::vector<coot::torsion_atom_quad>
coot::torsionable_quads(int imol, mmdb::Manager *mol, mmdb::PPAtom atom_selection,
			int n_selected_atoms,
			protein_geometry *geom_p) {

   bool pyranose_ring_torsion_flag = false; // no thanks
   std::vector<torsion_atom_quad> quads;
   std::vector<mmdb::Residue *> residues;
   for (int i=0; i<n_selected_atoms; i++) {
      mmdb::Residue *r = atom_selection[i]->residue;
      if (std::find(residues.begin(), residues.end(), r) == residues.end())
	 residues.push_back(r);
   }
   std::vector<torsion_atom_quad> link_quads =
      torsionable_link_quads(imol, residues, mol, geom_p);
   for (unsigned int iquad=0; iquad<link_quads.size(); iquad++)
      quads.push_back(link_quads[iquad]);
   for (unsigned int ires=0; ires<residues.size(); ires++) {
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      residues[ires]->GetAtomTable(residue_atoms, n_residue_atoms);
      std::vector<torsion_atom_quad> monomer_quads =
	 torsionable_bonds_monomer_internal_quads(residues[ires], residue_atoms,
							n_residue_atoms,
							pyranose_ring_torsion_flag, geom_p);
      for (unsigned int iquad=0; iquad<monomer_quads.size(); iquad++)
	 quads.push_back(monomer_quads[iquad]);
   }
   return quads;
}


// And the atom_quad version of that (for setting link torsions)
//
std::vector<coot::torsion_atom_quad>
coot::torsionable_link_quads(int imol,
			     std::vector<mmdb::Residue *> residues_in,
			     mmdb::Manager *mol, protein_geometry *geom_p) {

   if (false) {
      std::cout << "DEBUG:: torsionable_link_quads called with residues in size "
                << residues_in.size() << std::endl;
      for (unsigned int i=0; i<residues_in.size(); i++) {
         std::cout << "DEBUG::   " << coot::residue_spec_t(residues_in[i]) << std::endl;
      }
      logger.log(log_t::DEBUG, logging::ltw("torsionable_link_quads called with residues in size"),
                 logging::ltw(residues_in.size()));
      for (unsigned int i=0; i<residues_in.size(); i++)
         logger.log(log_t::DEBUG, residue_spec_t(residues_in[i]).format());
   }

   std::vector<torsion_atom_quad> quads;
   std::vector<std::pair<bool, mmdb::Residue *> > residues(residues_in.size());
   for (unsigned int i=0; i<residues_in.size(); i++)
      residues[i] = std::pair<bool, mmdb::Residue *> (0, residues_in[i]);

   // We want a quick way of getting to the restaints of an atom's
   // residue (link_atom_1 and link_atom_2 below).
   // So here we set up res_restraints map indexed by a mmdb::Residue *.
   //
   std::map<mmdb::Residue *, dictionary_residue_restraints_t> res_restraints;
   for (unsigned int ires=0; ires<residues_in.size(); ires++) {
      std::string rn = residues_in[ires]->GetResName();
      std::pair<bool, dictionary_residue_restraints_t> rest = geom_p->get_monomer_restraints(rn, imol);
      if (! rest.first) {
	 std::string m = "Restraints not found for type ";
	 m += rn;
	 throw std::runtime_error(m);
      }
      res_restraints[residues_in[ires]] = rest.second;
   }

   std::vector<atom_spec_t> dummy_fixed_atom_specs;
   std::vector<mmdb::Link> links;
   clipper::Xmap<float> dummy_xmap;
   restraints_container_t restraints(residues, links, *geom_p, mol, dummy_fixed_atom_specs, &dummy_xmap);
   bonded_pair_container_t bpc = restraints.bonded_residues_from_res_vec(*geom_p);

   for (unsigned int i=0; i<bpc.bonded_residues.size(); i++) {

      const dictionary_residue_link_restraints_t &link = geom_p->link(bpc[i].link_type);

      // In a NAG-ASN link the first residue should be the NAG and residue-2 should be the ASN.
      //
      // std::cout << "DEBUG:: link-idx: " << i << " link_id: \"" << link.link_id << "\"" << std::endl;

      if (link.link_id != "") {

	 // Don't use bonds - use link torsions if you can.
	 //
	 // std::cout << "DEBUG:: link " << link.link_id << " has " << link.link_torsion_restraint.size()
	 //           << " torsion restraints " << std::endl;

	 if (link.link_torsion_restraint.size() > 0) {

	    for (unsigned int il=0; il<link.link_torsion_restraint.size(); il++) {
	       // std::cout << "----------- link torsion restraint " << il << std::endl;
	       const dict_link_torsion_restraint_t &rest = link.link_torsion_restraint[il];
	       if (rest.is_pyranose_ring_torsion()) {
		  // pass
		  // std::cout << "   link # " << il << " is pyranose ring torsion # PASS" << std::endl;
                  logger.log(log_t::INFO, "   link #", il, "is pyranose ring torsion # PASS");

		  if (false) { // debug
		     mmdb::Residue *r_1 = bpc[i].res_1;
		     mmdb::Residue *r_2 = bpc[i].res_1;
		     mmdb::Residue *r_3 = bpc[i].res_1;
		     mmdb::Residue *r_4 = bpc[i].res_1;
		     if (rest.atom_1_comp_id == 2)
			r_1 = bpc[i].res_2;
		     if (rest.atom_2_comp_id == 2)
			r_2 = bpc[i].res_2;
		     if (rest.atom_3_comp_id == 2)
			r_3 = bpc[i].res_2;
		     if (rest.atom_4_comp_id == 2)
			r_4 = bpc[i].res_2;
		     mmdb::Atom *link_atom_1 = r_1->GetAtom(rest.atom_id_1_4c().c_str());
		     mmdb::Atom *link_atom_2 = r_2->GetAtom(rest.atom_id_2_4c().c_str());
		     mmdb::Atom *link_atom_3 = r_3->GetAtom(rest.atom_id_3_4c().c_str());
		     mmdb::Atom *link_atom_4 = r_4->GetAtom(rest.atom_id_4_4c().c_str());

		     // std::cout << "   link # " << il << " has link_atoms: "
		     //           << atom_spec_t(link_atom_1) << " " << atom_spec_t(link_atom_2) << " "
		     //           << atom_spec_t(link_atom_3) << " " << atom_spec_t(link_atom_4) << " "
		     //           << std::endl;
                     logger.log(log_t::INFO, {"   link #", il, "has link_atoms:",
                                              atom_spec_t(link_atom_1).format(),
                                              atom_spec_t(link_atom_2).format(),
                                              atom_spec_t(link_atom_3).format(),
                                              atom_spec_t(link_atom_4).format()});

		  }
	       } else {
		  mmdb::Residue *r_1 = bpc[i].res_1;
		  mmdb::Residue *r_2 = bpc[i].res_1;
		  mmdb::Residue *r_3 = bpc[i].res_1;
		  mmdb::Residue *r_4 = bpc[i].res_1;
		  if (rest.atom_1_comp_id == 2)
		     r_1 = bpc[i].res_2;
		  if (rest.atom_2_comp_id == 2)
		     r_2 = bpc[i].res_2;
		  if (rest.atom_3_comp_id == 2)
		     r_3 = bpc[i].res_2;
		  if (rest.atom_4_comp_id == 2)
		     r_4 = bpc[i].res_2;
		  mmdb::Atom *link_atom_1 = r_1->GetAtom(rest.atom_id_1_4c().c_str());
		  mmdb::Atom *link_atom_2 = r_2->GetAtom(rest.atom_id_2_4c().c_str());
		  mmdb::Atom *link_atom_3 = r_3->GetAtom(rest.atom_id_3_4c().c_str());
		  mmdb::Atom *link_atom_4 = r_4->GetAtom(rest.atom_id_4_4c().c_str());

		  // std::cout << "   link # " << il << " has residues    "
		  //           << residue_spec_t(r_1) << "         " << residue_spec_t(r_2)
		  //           << "         "
		  //           << residue_spec_t(r_3) << "         " << residue_spec_t(r_4)
		  //           << std::endl;
		  // std::cout << "   link # " << il << " has link_atoms: "
		  //           << atom_spec_t(link_atom_1) << " " << atom_spec_t(link_atom_2) << " "
		  //           << atom_spec_t(link_atom_3) << " " << atom_spec_t(link_atom_4) << " "
		  //           << std::endl;

                  logger.log(log_t::INFO, {"   link #", il, "has residues",
                                           residue_spec_t(r_1).format(),
                                           residue_spec_t(r_2).format(),
                                           residue_spec_t(r_3).format(),
                                           residue_spec_t(r_4).format()});
                  logger.log(log_t::INFO, {"   link #", il, "has link_atoms:",
                                           atom_spec_t(link_atom_1).format(),
                                           atom_spec_t(link_atom_2).format(),
                                           atom_spec_t(link_atom_3).format(),
                                           atom_spec_t(link_atom_4).format()});

		  if (link_atom_1 && link_atom_2 && link_atom_3 && link_atom_4) {
		     torsion_atom_quad q(link_atom_1, link_atom_2, link_atom_3, link_atom_4,
					 rest.angle(),
					 rest.angle_esd(),
					 rest.period());
		     q.name = rest.id();
		     quads.push_back(q);
		  }
	       }
	    }

	 } else {

	    // std::cout << "INFO:: link torsion generated from link bond restraint"  << std::endl;
	    logger.log(log_t::INFO, "link torsion generated from link bond restraint");

	    // bleugh... OK, no torsion restraints.
	    // So use a bond restaint to make one torsion (around the link bond).
	    //
	    for (unsigned int ib=0; ib<link.link_bond_restraint.size(); ib++) {
	       mmdb::Atom *link_atom_1 = bpc[i].res_1->GetAtom(link.link_bond_restraint[ib].atom_id_1_4c().c_str());
	       mmdb::Atom *link_atom_2 = bpc[i].res_2->GetAtom(link.link_bond_restraint[ib].atom_id_2_4c().c_str());
	       if (link_atom_1 && link_atom_2) {
		  // What are the neightbours of link_atom_1 (and link_atom_2)?
		  // Try to find a non-hydrogen atom to which it is bonded.
		  bool H_flag = false;
		  std::string atom_name_1 = link_atom_1->name;
		  std::string atom_name_2 = link_atom_2->name;
		  std::vector<std::string> n1;
		  std::vector<std::string> n2;
		  n1 = res_restraints[link_atom_1->residue].neighbours(atom_name_1, H_flag);
		  n2 = res_restraints[link_atom_2->residue].neighbours(atom_name_2, H_flag);
		  if (n1.size() && n2.size()) {
		     std::string neigbhour_1_name = n1[0];
		     std::string neigbhour_2_name = n2[0];
		     mmdb::Atom *n_at_1 = bpc[i].res_1->GetAtom(neigbhour_1_name.c_str());
		     mmdb::Atom *n_at_2 = bpc[i].res_2->GetAtom(neigbhour_2_name.c_str());
		     if (n_at_1 && n_at_2) {
 			torsion_atom_quad q(n_at_1, link_atom_1, link_atom_2, n_at_2,
					    -84, 40, 3); // synthetic values
			// If link is NAG-ASN then this torsion is phi (C1-ND2)
			q.name = "Bond-derived synthetic torsion-phi-" + util::int_to_string(ib);
			quads.push_back(q);

			// also we need psi, CB, CG, ND2, C1 (of NAG)
			//
			std::vector<std::string> n3;
			n3 = res_restraints[link_atom_2->residue].neighbours(atom_name_2, H_flag);
			if (n3.size()) {
			   mmdb::Atom *n_at_3 = bpc[i].res_2->GetAtom(n3[0].c_str());
			   if (n_at_3) {
			      torsion_atom_quad q(link_atom_1, link_atom_2, n_at_2, n_at_3,
						  185.1, 40, 3); // synthetic values
			      // If link is NAG-ASN then this torsion is phi (C1-ND2)
			      q.name = "Bond-derived synthetic torsion-psi-" + util::int_to_string(ib);
			      // quads.push_back(q);
			   }
			}
		     }
		  } 
	       } else {
		  std::cout << "WARNING:: oops missing link atoms " << std::endl;
		  if (! link_atom_1)
		     std::cout << "   " << link.link_bond_restraint[ib].atom_id_1_4c().c_str()
			       << " is missing from residue " << residue_spec_t(bpc[i].res_1) << std::endl;
		  if (! link_atom_2)
		     std::cout << "   " << link.link_bond_restraint[ib].atom_id_2_4c().c_str()
			       << " is missing from residue " << residue_spec_t(bpc[i].res_2) << std::endl;
	       }
	    }
	 }
      }
   }
   return quads;
}


// this can throw an exception.
//
// This presumes that mol is actually a mol-fragment, created something like:
//    mmdb::Manager *moving_mol =
//      coot::util::create_mmdbmanager_from_residue_specs(residue_specs, atom_sel.mol);
//
void
coot::multi_residue_torsion_fit_map(int imol,
				    mmdb::Manager *mol,
				    const clipper::Xmap<float> &xmap,
				    const std::vector<std::pair<bool, clipper::Coord_orth> > &avoid_these_atoms, // 20170613 flag is is-water
				    int n_trials,
				    protein_geometry *geom_p) {

   // First fill the atoms vector: all the atoms in the input mol
   // Get the torsionable quads atom names
   // Get the contact info and use that to make torsions
   // For n_trials,
   //     make a model using set_dihedral_multi
   //     get self_clash_score
   //     if self_clash_score < 6
   //        get density fit score
   //        if its better than the best score so far
   //           update best_quads and best_tree_dihedral_quads
   // update model so that it uses best_tree_dihedral_quads

   std::vector<std::pair<std::string, int> > atom_numbers = util::atomic_number_atom_list();

   try {
      mmdb::PPAtom atom_selection = 0;
      int n_selected_atoms;
      int selhnd = mol->NewSelection(); // d
      mol->SelectAtoms(selhnd, 0, "*",
		       mmdb::ANY_RES, "*",
		       mmdb::ANY_RES, "*",
		       "*", "*", "*", "*");
      mol->GetSelIndex(selhnd, atom_selection, n_selected_atoms);
      std::vector<std::pair<mmdb::Atom *, float> > atoms(n_selected_atoms); // for density fitting
      for (int iat=0; iat<n_selected_atoms; iat++) {
	 int atomic_number = util::atomic_number(atom_selection[iat]->element, atom_numbers);
	 float z = atomic_number;
	 if (atomic_number == -1)
	    z = 6.0f;
	 atoms[iat] = std::pair<mmdb::Atom *, float> (atom_selection[iat], z);
      }

      if (n_selected_atoms > 0) {
	 std::vector<torsion_atom_quad> quads =
	    torsionable_quads(imol, mol, atom_selection, n_selected_atoms, geom_p);

	 // FIXME for future, calculate link_angle_atom_triples, using something analoguous to
	 // torsionable_link_quads()

	 if (false) // debug
	    for (unsigned int iquad=0; iquad<quads.size(); iquad++)
	       std::cout << "DEBUG multi-residue-torsion-fit-map: tosion quads:  "
                         << iquad << " "
			 << atom_spec_t(quads[iquad].atom_1) << " "
			 << atom_spec_t(quads[iquad].atom_2) << " "
			 << atom_spec_t(quads[iquad].atom_3) << " "
			 << atom_spec_t(quads[iquad].atom_4) << " \""
			 << quads[iquad].name << "\" torsion: "
			 << quads[iquad].torsion()
			 << std::endl;

	 contact_info contacts(mol, imol, selhnd, quads, geom_p);
	 std::vector<std::vector<int> > contact_indices =
	    contacts.get_contact_indices_with_reverse_contacts();
	 atom_tree_t tree(contact_indices, 0, mol, selhnd);

	 bool reverse_flag = 1;
	 double pre_score = util::z_weighted_density_score_new(atoms, xmap);

	 double best_score = pre_score;
	 int n_quads = quads.size();
	 std::vector<double> best_quads(n_quads, -1);
	 std::vector<atom_tree_t::tree_dihedral_quad_info_t> best_tree_dihedral_quads;
	 map_index_t fixed_index(0);

	 // save the current
	 for (int iquad=0; iquad<n_quads; iquad++)
	    best_quads[iquad] = quads[iquad].torsion();

	 // Make the first few shifts small, because we could be close
	 // to the correct solution by initial placement. At a guess,
	 // the first 15% should be small.
	 //
	 int itrial_n_small_lim(0.15*n_trials);
	 //
	 for (int itrial=0; itrial<n_trials; itrial++) {

	    bool allow_conformer_switch = true;
	    bool small_torsion_changes = false;
	    if (itrial < itrial_n_small_lim) {
	       allow_conformer_switch = false;
	       small_torsion_changes = true;
	    }

	    if (false)
	       std::cout << "Round " << itrial << " of " << n_trials << " for " << n_quads << " quads "
			 << std::endl;

	    std::vector<atom_tree_t::tree_dihedral_quad_info_t> torsion_quads;

	    // debug, store angles in rand_angles: name current trial-value
	    std::vector<std::tuple<std::string, double, double> > rand_angles(n_quads);

	    for (int iquad=0; iquad<n_quads; iquad++) {
	       // quads[iquad] is passed for debugging
               double rand_angle = get_rand_angle(best_quads[iquad], quads[iquad], itrial,
						  n_trials, allow_conformer_switch, small_torsion_changes);

	       std::tuple<std::string, double, double> tup(quads[iquad].name, best_quads[iquad], rand_angle);
	       rand_angles[iquad] = tup;

	       atom_tree_t::tree_dihedral_quad_info_t tor(quads[iquad], rand_angle, fixed_index);
	       torsion_quads.push_back(tor);
	    }


	    if (false) { //debug
	       for (int iquad=0; iquad<n_quads; iquad++) {
		  std::cout << "debug: itrial " << itrial << " "
			    << "iquad " << iquad << " "
			    << std::get<0>(rand_angles[iquad]) << " "
			    << std::get<1>(rand_angles[iquad]) << " "
			    << std::get<2>(rand_angles[iquad]) << " ";
	       }
	       std::cout << std::endl;
	    }

            if (false) {
               std::string file_name = "A-trial-" + std::to_string(itrial) + ".pdb";
               mol->WritePDBASCII(file_name.c_str());
            }

	    tree.set_dihedral_multi(torsion_quads);

            if (false) {
               std::string file_name = "B-trial-" + std::to_string(itrial) + ".pdb";
               mol->WritePDBASCII(file_name.c_str());
            }

	    // FIXME for futures, also include link_angle_atom_triples (for excluding of bumps)
	    double self_clash_score = get_self_clash_score(mol, atom_selection, n_selected_atoms, quads);

	    double env_clash_score = get_environment_clash_score(mol, atom_selection, n_selected_atoms,
                                                                 avoid_these_atoms);

	    if (false) {
	       std::cout << "DEBUG:: self_clash_score: " << self_clash_score << std::endl;
	       std::cout << "DEBUG::  env_clash_score: " <<  env_clash_score << std::endl;
	    }

	    // self-clash scores have mean 7.5, median 3.3 and sd 14, IRQ 0.66
	    // Is this a good clash score lim?  Not clear, but 10.0 is better than 1.0
	    //
	    if ((self_clash_score > 6) || (env_clash_score > 30.0)) {

	       // crash and bangs into itself (between residues)
	       // or into its neighbours (the 1.0 might need tuning)

	    } else {

	       // happy path

	       double this_score = util::z_weighted_density_score_new(atoms, xmap);

	       // debugging of scores
	       if (false) {
		  std::cout << "debug trial " << itrial << " fit-score: " << this_score
			    << " self-clash-score " << self_clash_score
			    << " for quads ";
		  for (unsigned int iquad=0; iquad<quads.size(); iquad++)
		     std::cout << "   " << quads[iquad].torsion();
		  std::cout << std::endl;
	       }

	       if (this_score > best_score) {
                  // std::cout << "Round " << itrial << " improved! was " << best_score
                  // << " now " << this_score << std::endl;
                  logger.log(log_t::DEBUG, {std::string("Round"), itrial,
                                            std::string("improvement - was"), best_score,
                                            std::string("now"), this_score});
                  // util::debug_z_weighted_density_score_new(atoms,xmap);
		  // save best torsion angles
		  best_score = this_score;
		  for (int iquad=0; iquad<n_quads; iquad++)
		     best_quads[iquad] = quads[iquad].torsion();
		  best_tree_dihedral_quads = torsion_quads;
	       }

	       if (false) {  // debugging.

		  // set the b-factor of the atoms to the score
		  int imod = 1;
		  mmdb::Model *model_p = mol->GetModel(imod);
		  int n_chains = model_p->GetNumberOfChains();
		  for (int ichain=0; ichain<n_chains; ichain++) {
                     mmdb::Chain *chain_p = model_p->GetChain(ichain);
		     int nres = chain_p->GetNumberOfResidues();
		     for (int ires=0; ires<nres; ires++) {
                        mmdb::Residue *residue_p = chain_p->GetResidue(ires);
			int n_atoms = residue_p->GetNumberOfAtoms();
			for (int iat=0; iat<n_atoms; iat++) {
                           mmdb::Atom *at = residue_p->GetAtom(iat);
			   at->tempFactor = this_score * 0.4;
			   at->tempFactor = self_clash_score;
			}
		     }
		  }
		  std::string file_name = "trial-" + util::int_to_string(itrial) + ".pdb";
		  mol->WritePDBASCII(file_name.c_str());
	       }
	    }
	 }
         if (!best_tree_dihedral_quads.empty()) {
            // std::string fn = "pre-set-dihedrals-for-best.pdb";
            // mol->WritePDBASCII(fn.c_str());
            tree.set_dihedral_multi(best_tree_dihedral_quads);
            // fn = "post-set-dihedrals-for-best.pdb";
            //mol->WritePDBASCII(fn.c_str());
         }
	 mol->DeleteSelection(selhnd);
      }
   }
   catch (const std::runtime_error &rte) {
      std::cout << "WARNING:: " << rte.what() << std::endl;
   }
}


double
coot::get_rand_angle(double current_angle,
                     const torsion_atom_quad &quad,
                     int itrial, int n_trials,
		     bool allow_conformer_switch,
		     bool small_torsion_changes) {

   double r = current_angle;
   double minus_one_to_one = -1 + 2 * double(util::random())/double(RAND_MAX);
   // trial_factor goes from 0 (start) to 1 (end)
   double trial_factor = double(itrial)/double(n_trials);
   // double angle_scale_factor = 0.2 + 0.8*(1-trial_factor);
   double angle_scale_factor = 0.2 + 0.8 - trial_factor;

   if (small_torsion_changes) {
      r += 5.0 * minus_one_to_one;
   } else {
      r += 30 * minus_one_to_one * angle_scale_factor;
   }

   // allow gauche+/gauche-/trans
   if (allow_conformer_switch) {
      double rn = float(util::random())/float(RAND_MAX);
      double tf = 1 - trial_factor; // tf goes from 1 (start) to 0 (end)
      if (rn < (0.02 + 0.25 * tf)) {
	 double rn_2 = float(util::random())/float(RAND_MAX);
	 double step = floor(6 * rn_2) * 60.0;
	 // std::cout << "      step " << step << std::endl;
	 r += step;
      }
   }

   if (r > 360)
      r -= 360;
   
   return r; 
} 


double
coot::get_self_clash_score(mmdb::Manager *mol,
			   mmdb::PPAtom atom_selection,
			   int n_selected_atoms,
			   const std::vector<torsion_atom_quad> &quads) {

   // Score is
   // sum of (d-bump_max)^2 for atom pairs i,j where j<i where d < bump_max

   mmdb::realtype bump_max = 3.6; // find distances between atoms that are less than this.
   bump_max = 2.8; // 20170615 try this (for less self bumping)
   double clash_score = 0;

   // setup for SeekContacts():
   // 
   mmdb::Contact *pscontact = NULL;
   int n_contacts;
   long i_contact_group = 1;
   mmdb::mat44 my_matt;
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
         my_matt[i][j] = 0.0;      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

   mol->SeekContacts(atom_selection, n_selected_atoms,
		     atom_selection, n_selected_atoms,
		     0.001, bump_max,
		     0, // seqDist: 1 means in different residues,
		        // but if I set that, then no contacts are found
		        // even for atoms that are in different residues:
		        // mmdb bug.
		     pscontact, n_contacts,
		     0, &my_matt, i_contact_group);
   
   if (n_contacts > 0) {
      if (pscontact) {
         for (int i=0; i<n_contacts; i++) {
	    if (pscontact[i].id1 < pscontact[i].id2) { 
	       mmdb::Atom *at_1 = atom_selection[pscontact[i].id1];
	       mmdb::Atom *at_2 = atom_selection[pscontact[i].id2];
	       if (at_1->residue != at_2->residue) {
		  std::string e1 = at_1->element;
		  std::string e2 = at_2->element;

		  if ((e1 != " H") && (e2 != " H")) {  // PDB vs 3 FIXME
		     // ignore bumps to O5 (e.g. O4(prev)-O5(new)) on newly added residue

		     std::string atom_name_2 = at_2->name;
		     if (atom_name_2 != " O5 ") {
			double d_sqd =
			   (at_1->x-at_2->x) * (at_1->x-at_2->x) +
			   (at_1->y-at_2->y) * (at_1->y-at_2->y) + 
			   (at_1->z-at_2->z) * (at_1->z-at_2->z);

			// are they either in a bond, angle or torsion of any of quads?
			// 
			bool in_a_tors = both_in_a_torsion_p(at_1, at_2, quads);
			if (! in_a_tors) {
			   double delta = bump_max - sqrt(d_sqd);
			   clash_score += delta * delta;
			   if (false)
			      std::cout << "adding to clash_score " << delta * delta << " for dist " << sqrt(d_sqd)
					<< " between " << atom_spec_t(at_1) << " and " << atom_spec_t(at_2) << std::endl;
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   return clash_score;
}


bool
coot::both_in_a_torsion_p(mmdb::Atom *at_1,
			  mmdb::Atom *at_2,
			  const std::vector<torsion_atom_quad> &quads) {

   // this doesn't check angles and (non-dictionary) torsion that are
   // made as a results of the link but not in a torsion.

   bool in_a_tors = false;
   for (unsigned int i=0; i<quads.size(); i++) {
      bool found_at_1 = false;
      bool found_at_2 = false;
      const torsion_atom_quad &q = quads[i];
      if (q.atom_1 == at_1) found_at_1 = true;
      if (q.atom_2 == at_1) found_at_1 = true;
      if (q.atom_3 == at_1) found_at_1 = true;
      if (q.atom_4 == at_1) found_at_1 = true;
      if (q.atom_1 == at_2) found_at_2 = true;
      if (q.atom_2 == at_2) found_at_2 = true;
      if (q.atom_3 == at_2) found_at_2 = true;
      if (q.atom_4 == at_2) found_at_2 = true;
      if (found_at_1 && found_at_2) {
	 in_a_tors = true;
	 break;
      }
   }
   return in_a_tors;
}

// return a positive number for a clash - the bigger the number the worse the clash.
//
double
coot::get_environment_clash_score(mmdb::Manager *mol,
				  mmdb::PPAtom atom_selection,
				  int n_selected_atoms,
				  const std::vector<std::pair<bool, clipper::Coord_orth> > &avoid_these_atoms) {
   double cs = 0;
   double sf = 1.0;
   for (int iat=0; iat<n_selected_atoms; iat++) {
      // we expect that the ASN will be close to its polypeptide neighbours.  We don't want to
      // include such clashes
      std::string res_name = atom_selection[iat]->GetResName();
      // std::cout << "res_name is " << res_name << std::endl;
      if (res_name != "ASN") {
	 clipper::Coord_orth at_pt = co(atom_selection[iat]);
	 for (unsigned int jat=0; jat<avoid_these_atoms.size(); jat++) {
	    double close_lim = 3.3;
	    if (avoid_these_atoms[jat].first) close_lim = 2.5; // we can get close to waters without worry
	    double close_lim_sqrd = close_lim * close_lim;
	    double d_sqd = (at_pt - avoid_these_atoms[jat].second).lengthsq();
	    if (d_sqd < close_lim_sqrd) {
	       double diff = close_lim - sqrt(d_sqd);
	       cs += diff*diff*sf;
	       if (false)
		  std::cout << "DEBUG:: env clash: atom " << atom_spec_t(atom_selection[iat]) << " is close to "
			    << jat << " " << avoid_these_atoms[jat].second.format() << " " << sqrt(d_sqd) << std::endl;
	    }
	 }
      }
   }
   return cs;
}
