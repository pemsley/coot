/* coot-utils/coot-coord-extras.cc
 * 
 * Copyright 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2009 by The University of Oxford
 * Author: Paul Emsley
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


#include <stdexcept>
#include <sstream>
#include <queue>

#include "coot-utils.hh"
#include "coot-coord-utils.hh"
#include "coot-coord-extras.hh"

#include "coot-sysdep.h"


// Return 0 if any of the residues don't have a dictionary entry
// geom_p gets updated to include the residue restraints if necessary
// 
std::pair<int, std::vector<std::string> >
coot::util::check_dictionary_for_residues(PCResidue *SelResidues, int nSelResidues,
					  coot::protein_geometry *geom_p,
					  int read_number) {

   std::pair<int, std::vector<std::string> > r;

   int status;
   int fail = 0; // not fail initially.

   for (int ires=0; ires<nSelResidues; ires++) { 
      std::string resname(SelResidues[ires]->name);
      status = geom_p->have_dictionary_for_residue_type(resname, read_number);
      // This bit is redundant now that try_dynamic_add has been added
      // to have_dictionary_for_residue_type():
      if (status == 0) { 
	 status = geom_p->try_dynamic_add(resname, read_number);
	 if (status == 0) {
	    fail = 1; // we failed to find it then.
	    r.second.push_back(resname);
	 }
      }
   }

   if (fail)
      r.first = 0;
   return r;
} 





// We also now pass regular_residue_flag so that the indexing of the
// contacts is inverted in the case of not regular residue.  I don't
// know why this is necessary, but I have stared at it for hours, this
// is a quick (ugly hack) fix that works.  I suspect that there is
// some atom order dependency in mgtree that I don't understand.
// Please fix (remove the necessity of depending on
// regular_residue_flag) if you know how.
// 
std::vector<std::vector<int> >
coot::util::get_contact_indices_from_restraints(CResidue *residue,
						coot::protein_geometry *geom_p,
						short int regular_residue_flag) {

   int nResidueAtoms = residue->GetNumberOfAtoms(); 
   std::vector<std::vector<int> > contact_indices(nResidueAtoms);
   std::string restype(residue->name);
   CAtom *atom_p;

   int n_restr = geom_p->size();

   for (int icomp=0; icomp<n_restr; icomp++) {
      if ((*geom_p)[icomp].comp_id == restype) {
// 	 std::cout << "There are " << (*geom_p)[icomp].bond_restraint.size()
// 		   << " bond restraints " << "for " << restype << std::endl;
	 for (unsigned int ibr=0; ibr< (*geom_p)[icomp].bond_restraint.size(); ibr++) {
	    for (int iat=0; iat<nResidueAtoms; iat++) {
	       atom_p = residue->GetAtom(iat);
	       std::string at_name(atom_p->GetAtomName());
	       if ( (*geom_p)[icomp].bond_restraint[ibr].atom_id_1_4c() == at_name ) {
// 		  std::cout << "found a bond match "
// 			    << (*geom_p)[icomp].bond_restraint[ibr].atom_id_1_4c()
// 			    << " to "
// 			    << (*geom_p)[icomp].bond_restraint[ibr].atom_id_2_4c()
// 			    << std::endl;
		  int ibond_to = -1;  // initially unassigned.
		  std::string at_name_2;
		  for (int iat2=0; iat2<nResidueAtoms; iat2++) {
		     atom_p = residue->GetAtom(iat2);
		     at_name_2 = atom_p->GetAtomName();
		     if ( (*geom_p)[icomp].bond_restraint[ibr].atom_id_2_4c() == at_name_2 ) {
			ibond_to = iat2;
			break;
		     }
		  }
		  if (ibond_to > -1 ) { 
		     if (regular_residue_flag) {
			contact_indices[iat].push_back(ibond_to);  // for ALA etc
		     } else {
			contact_indices[ibond_to].push_back(iat);  // ligands
			// contact_indices[iat].push_back(ibond_to);  // ALA etc
		     }
		  } 
//		  else
		     // This spits out the names of Hydrogens, often.
//  		     std::cout << "failed to find bonded atom "
//  			       << (*geom_p)[icomp].bond_restraint[ibr].atom_id_2_4c()
//  			       << std::endl;
	       }
	    }
	 }
      }
   }
   return contact_indices;
}

// The atoms of residue_atoms are in the "right" order for not making
// a tree along the main chain.
// 
std::vector<std::vector<int> >
coot::util::get_contact_indices_for_PRO_residue(PPCAtom residue_atoms,
						int nResidueAtoms, 
						coot::protein_geometry *geom_p) { 

   std::vector<std::vector<int> > contact_indices(nResidueAtoms);
   CAtom *atom_p;
   int n_restr = geom_p->size();
   for (int icomp=0; icomp<n_restr; icomp++) {
      if ((*geom_p)[icomp].comp_id == "PRO") {
	 for (unsigned int ibr=0; ibr< (*geom_p)[icomp].bond_restraint.size(); ibr++) {
	    for (int iat=0; iat<nResidueAtoms; iat++) {
	       atom_p = residue_atoms[iat];
	       std::string at_name(atom_p->GetAtomName());
	       if ( (*geom_p)[icomp].bond_restraint[ibr].atom_id_1_4c() == at_name ) {
		  int ibond_to = -1;  // initially unassigned.
		  std::string at_name_2;
		  for (int iat2=0; iat2<nResidueAtoms; iat2++) {
		     atom_p = residue_atoms[iat2];
		     at_name_2 = atom_p->GetAtomName();
		     if ( (*geom_p)[icomp].bond_restraint[ibr].atom_id_2_4c() == at_name_2 ) {
			ibond_to = iat2;
			break;
		     }
		  }
		  if (ibond_to != -1)
		     contact_indices[iat].push_back(ibond_to);
	       }
	    }
	 }
      }
   }
   return contact_indices;
}


coot::util::dict_residue_atom_info_t::dict_residue_atom_info_t(const std::string &residue_name_in,
							       coot::protein_geometry *geom_p) {

   residue_name = residue_name_in;

   std::pair<short int, dictionary_residue_restraints_t> p = 
      geom_p->get_monomer_restraints(residue_name);

   if (p.first) {
      for (unsigned int iat=0; iat<p.second.atom_info.size(); iat++) {
	 std::string atom_name = p.second.atom_info[iat].atom_id_4c;
	 short int isHydrogen = 0;
	 if (p.second.atom_info[iat].type_symbol == "H" ||
	     p.second.atom_info[iat].type_symbol == "D") {
	    isHydrogen = 1;
	 }
	 atom_info.push_back(coot::util::dict_atom_info_t(atom_name, isHydrogen));
      }
   }

}

// This one we can do a dynamic add.
// 
short int
coot::util::is_nucleotide_by_dict_dynamic_add(CResidue *residue_p, coot::protein_geometry *geom_p) {

   short int is_nuc = 0;
   short int ifound = 0;
   std::string residue_name = residue_p->GetResName();

   int n_restr = geom_p->size();
   for (int icomp=0; icomp<n_restr; icomp++) {
      if ((*geom_p)[icomp].comp_id == residue_name) {
	 ifound = 1;
	 if ((*geom_p)[icomp].residue_info.group == "RNA" ||
	     (*geom_p)[icomp].residue_info.group == "DNA" ) {
	    is_nuc = 1;
	 }
	 break;
      }
   }

   int read_number = 40;
   if (ifound == 0) {
      int status = geom_p->try_dynamic_add(residue_name, read_number);
      if (status != 0) {
	 // we successfully added it, let's try to run this function
	 // again.  Or we could just test the last entry in
	 // geom_p->dict_res_restraints(), but it is not public, so
	 // it's messy.
	 // 
	 is_nuc = is_nucleotide_by_dict_dynamic_add(residue_p, geom_p);
      } 
   }

   return is_nuc;
}


// This one we can NOT do a dynamic add.
//
short int
coot::util::is_nucleotide_by_dict(CResidue *residue_p, const coot::protein_geometry &geom) {

   short int is_nuc = 0;
   std::string residue_name = residue_p->GetResName();

   int n_restr = geom.size();
   for (int icomp=0; icomp<n_restr; icomp++) {
      if (geom[icomp].comp_id == residue_name) {
	 if (geom[icomp].residue_info.group == "RNA" ||
	     geom[icomp].residue_info.group == "DNA" ) {
	    is_nuc = 1;
	 }
	 break;
      }
   }

   return is_nuc;
}



coot::atom_tree_t::atom_tree_t(const coot::dictionary_residue_restraints_t &rest, CResidue *res,
			       const std::string &altconf) {

   if (! res) {
      std::string mess = "Null residue in atom tree constructor";
      throw std::runtime_error(mess);
   }
   residue = res; // class data now, used in rotate_about() - so that
		  // we don't have to pass the CResidue * again.

   if (rest.tree.size() == 0) {
      std::string mess = "No tree in restraint";
      throw std::runtime_error(mess);
   } 

   // Fill bonds
   PPCAtom residue_atoms = 0;
   int n_residue_atoms;
   res->GetAtomTable(residue_atoms, n_residue_atoms);
   // fill the bonds
   for (unsigned int i=0; i<rest.bond_restraint.size(); i++) {
      int idx1 = -1;
      int idx2 = -1;
      for (int iat=0; iat<n_residue_atoms; iat++) {
	 std::string atom_name = residue_atoms[iat]->name;
	 std::string atom_altl = residue_atoms[iat]->altLoc;
	 if (atom_name == rest.bond_restraint[i].atom_id_1())
	    if (atom_altl == "" || atom_altl == altconf)
	       idx1 = iat;
	 if (atom_name == rest.bond_restraint[i].atom_id_2())
	    if (atom_altl == "" || atom_altl == altconf)
	       idx2 = iat;
	 // OK, we have found them, no need to go on searching.
	 if ((idx1 != -1) && (idx2 != -1))
	    break;
      }

      if ((idx1 != -1) && (idx2 != -1))
	 bonds.push_back(std::pair<int,int>(idx1, idx2));
   }

   fill_name_map(altconf);

   bool success_vertex = fill_atom_vertex_vec(rest, res, altconf);
   if (! success_vertex) {
      std::string mess = "Failed to fill atom vector from cif atom tree - bad tree?";
      throw std::runtime_error(mess);
   }

   bool success_torsion = fill_torsions(rest, res, altconf);

}

void
coot::atom_tree_t::fill_name_map(const std::string &altconf) {
   
   PPCAtom residue_atoms = 0;
   int n_residue_atoms;
   residue->GetAtomTable(residue_atoms, n_residue_atoms);

   // atom-name -> index, now class variable
   // std::map<std::string, int, std::less<std::string> > name_to_index;
   for (int iat=0; iat<n_residue_atoms; iat++) {
      std::string atom_name(residue_atoms[iat]->name);
      std::string atom_altl = residue_atoms[iat]->altLoc;
      if (atom_altl == "" || atom_altl == altconf)
	 name_to_index[atom_name] = iat;
   }
}


// the constructor, given a list of bonds and a base atom index.
// Used perhaps as the fallback when the above raises an
// exception.
coot::atom_tree_t::atom_tree_t(const std::vector<std::vector<int> > &contact_indices,
			       int base_atom_index, 
			       CResidue *res,
			       const std::string &altconf) {
   if (! res) {
      std::string mess = "null residue in alternate atom_tree_t constructor";
      throw std::runtime_error(mess);
   } else { 
      residue = res;
      fill_name_map(altconf);
      fill_atom_vertex_vec_using_contacts(contact_indices, base_atom_index);
   } 

}


bool
coot::atom_tree_t::fill_torsions(const coot::dictionary_residue_restraints_t &rest, CResidue *res,
				 const std::string &altconf) {


   bool r = 0;
   int n_torsions_inserted = 0;
   if (rest.torsion_restraint.size() > 0) {
      std::vector<coot::atom_index_quad> quads;
      for (unsigned int itr=0; itr<rest.torsion_restraint.size(); itr++) { 
	 try {
	    coot::dict_torsion_restraint_t tr = rest.torsion_restraint[itr];
	    coot::atom_index_quad quad = get_atom_index_quad(tr, res, altconf);
	    quads.push_back(quad);
	 }
	 catch (std::runtime_error rte) {
	    // We don't want to know that we can't fill hydrogen
	    // torsions when there are no hydrogens in the model.
	    // std::cout << rte.what() << std::endl;
	 }
      }

      // Now we have a set of atom index quads, put them in the
      // atom_vertex_vec, at the position of the second atom in the
      // torsion (the position of the first atom in the rotation
      // vector).

      // std::cout << " debug:: " << quads.size() << " quads" << std::endl;
      for (unsigned int iquad=0; iquad<quads.size(); iquad++) {
	 bool inserted = 0;
	 for (unsigned int iv=0; iv<atom_vertex_vec.size(); iv++) {
	    if (iv == quads[iquad].index2) {
	       for (unsigned int ifo=0; ifo<atom_vertex_vec[iv].forward.size(); ifo++) { 
		  if (atom_vertex_vec[iv].forward[ifo] == quads[iquad].index3) {

		     // now check that the forward atom of this
		     // forward atom is index4
		     int this_forward = atom_vertex_vec[iv].forward[ifo];
		     for (unsigned int ifo2=0; ifo2<atom_vertex_vec[this_forward].forward.size(); ifo2++) {
			if (atom_vertex_vec[this_forward].forward[ifo2] == quads[iquad].index4) {
			   atom_vertex_vec[iv].torsion_quad.first = 1;
			   atom_vertex_vec[iv].torsion_quad.second = quads[iquad];
			   if (0) {
			      std::cout << "            DEBUG:: inserting torsion " << iquad
					<< " into vertex " << iv << std::endl;
			   }
			   r = 1;
			   inserted = 1;
			   n_torsions_inserted++;
			}
		     } 
		  }
	       }
	    }
	    if (inserted)
	       break;
	 }
      } 
   }
   // std::cout << "DEBUG:: inserted " << n_torsions_inserted << " torsions" << std::endl;
   return r;
}

bool
coot::atom_tree_t::fill_atom_vertex_vec_using_contacts(const std::vector<std::vector<int> > &contact_indices,
						       int base_atom_index) {

   bool r=0;

   PPCAtom residue_atoms;
   int n_residue_atoms;
   residue->GetAtomTable(residue_atoms, n_residue_atoms);
   atom_vertex_vec.resize(n_residue_atoms);

   coot::atom_vertex av;
   av.connection_type = coot::atom_vertex::START;
   atom_vertex_vec[base_atom_index] = av;


   if (0) { 
      std::cout << " debug:: =========== contact indices in -======= " << std::endl;
      for (unsigned int ic1=0; ic1<contact_indices.size(); ic1++) {
	 std::cout << " index " << ic1 << " : ";
	 for (unsigned int ic2=0; ic2<contact_indices[ic1].size(); ic2++)
	    std::cout << contact_indices[ic1][ic2] << " ";
	 std::cout << std::endl;
      }
   }


   std::queue<int> q;
   q.push(base_atom_index);
   std::vector<int> done; // things that were in the queue that are
			  // not done (on poping). So that we don't
			  // loop endlessly doing atom indices that we
			  // have already considered.

   while (q.size()) {
      int this_base_atom = q.front();
      // now what are the forward atoms of av?
      std::vector<int> av_contacts = contact_indices[this_base_atom];
      
      for (unsigned int iav=0; iav<av_contacts.size(); iav++) {

	 int i_forward = av_contacts[iav];

	 // Check that this forward atom is not already in the forward
	 // atoms of this_base_atom
	 bool ifound_forward = 0;
	 for (unsigned int ifo=0; ifo<atom_vertex_vec[this_base_atom].forward.size(); ifo++) {
	    if (atom_vertex_vec[this_base_atom].forward[ifo] == av_contacts[iav]) { 
	       ifound_forward = 1;
	       break;
	    }
	 }

	 if (! ifound_forward) { 
	    // Add the contact as a forward atom of this_base_atom but
	    // only if the forward atom does not have a forward atom
	    // which is this_base_atom (and this keeps the tree going in
	    // up and not back again).
	    bool ifound_forward_forward = 0;
	    for (unsigned int ifo=0; ifo<atom_vertex_vec[i_forward].forward.size(); ifo++) {
	       if (atom_vertex_vec[i_forward].forward[ifo] == this_base_atom) {
		  ifound_forward_forward = 1;
		  break;
	       }
	    }
	    if (!ifound_forward_forward) {
	       atom_vertex_vec[this_base_atom].forward.push_back(av_contacts[iav]);
	    }
	 }

	 // add the forward atoms to the queue if they are not already
	 // in the done list.
	 bool in_done = 0;
	 for (unsigned int idone=0; idone<done.size(); idone++) {
	    if (done[idone] == av_contacts[iav]) {
	       in_done = 1;
	       break;
	    } 
	 }
	 if (!in_done)
	    q.push(av_contacts[iav]);
      }
      
      // if the forward atoms of this_atom_index do not have a back
      // atom, add one (but not, of course, if forward atom is the
      // start point of the tree).
      for (unsigned int iav=0; iav<av_contacts.size(); iav++) {
	 if (atom_vertex_vec[av_contacts[iav]].backward.size() == 0) {
	    if (atom_vertex_vec[av_contacts[iav]].connection_type != coot::atom_vertex::START) 
	       atom_vertex_vec[av_contacts[iav]].backward.push_back(this_base_atom);
	 }
      }
      q.pop();
      done.push_back(this_base_atom);
      r = 1; // return success
   }
   
   // print out the name_to_index map
   if (0) {
      std::cout << "==== atom indexes ===\n";
      for(std::map<std::string, coot::atom_tree_t::atom_tree_index_t>::const_iterator it = name_to_index.begin();
	  it != name_to_index.end(); ++it)
	 std::cout << "Atom :" << it->first << ": Index: " << it->second.index() << '\n';
   }

   // print out the atom tree
   if (0) {
      std::cout << "debug:: ==== atom_vertex_vec === from fill_atom_vertex_vec_using_contacts "
		<< std::endl;
      for (unsigned int iv=0; iv<atom_vertex_vec.size(); iv++) {
	 std::cout << "   atom_vertex_vec[" << iv << "] forward atom ";
	 for (unsigned int ifo=0; ifo<atom_vertex_vec[iv].forward.size(); ifo++) 
	    std::cout << atom_vertex_vec[iv].forward[ifo] << " ";
	 std::cout << std::endl;
      }
   }

   return r;
} 


// Throw an exception on not able to fill.
// 
coot::atom_index_quad
coot::atom_tree_t::get_atom_index_quad(const coot::dict_torsion_restraint_t &tr,
				       CResidue *res,
				       const std::string &altconf) const {
   coot::atom_index_quad quad(-1,-1,-1,-1);
   PPCAtom residue_atoms = 0;
   int n_residue_atoms;
   res->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      std::string atom_name(   residue_atoms[iat]->name);
      std::string atom_altconf(residue_atoms[iat]->altLoc);
      if (atom_name == tr.atom_id_1_4c())
	 if (atom_altconf == "" || atom_altconf == altconf)
	    quad.index1 = iat;
      if (atom_name == tr.atom_id_2_4c())
	 if (atom_altconf == "" || atom_altconf == altconf)
	    quad.index2 = iat;
      if (atom_name == tr.atom_id_3_4c())
	 if (atom_altconf == "" || atom_altconf == altconf)
	    quad.index3 = iat;
      if (atom_name == tr.atom_id_4_4c())
	 if (atom_altconf == "" || atom_altconf == altconf)
	    quad.index4 = iat;
   }
   if ((quad.index1 == -1) || (quad.index2 == -1) ||
       (quad.index3 == -1) || (quad.index4 == -1)) {
      std::string mess = "Can't fill atom_quad with indices ";
      std::ostringstream o;
      o << tr;
      mess += o.str();
      throw std::runtime_error(mess);
   } 
   return quad;
}


bool
coot::atom_tree_t::fill_atom_vertex_vec(const coot::dictionary_residue_restraints_t &rest,
					CResidue *res,
					const std::string &altconf) {

   // Note that we add an extra check on the forward and back atom
   // indices before adding them.  This is a sanity check - if we come
   // here with a big tree but a small residue, then the forward and
   // back atom indices could go beyond the limits of the size of the
   // atom_vertex_vec.  Which would be badness.  Hmm.. I am not sure
   // about this now. The indices are into the residue table, are they
   // not?  And the residue table corresponds here to residue_atoms
   // and n_residue_atoms. Hmm.. maybe that was not the problem.

   bool retval = 0; // fail initially.

   bool found_start = 0;
   int rest_tree_start_index = -1;
   for (unsigned int i=0; i<rest.tree.size(); i++) {
      if (rest.tree[i].connect_type == "START") { 
	 found_start = 1;
	 rest_tree_start_index = i;
	 break;
      }
   }

   // print out the name_to_index map
   if (0) {
      std::cout << "==== atom indexes ===\n";
      for(std::map<std::string, coot::atom_tree_t::atom_tree_index_t>::const_iterator it = name_to_index.begin();
	  it != name_to_index.end(); ++it)
	 std::cout << "Atom :" << it->first << ": Index: " << it->second.index() << '\n';
   }
 

   if (found_start) {

      // atom vertex is based on residue atoms (not dict).
      int n_residue_atoms = res->GetNumberOfAtoms();
      atom_vertex_vec.resize(n_residue_atoms);

      // fill atom_vertex_vec, converting dictionary atom names to
      // residue indices.
      // 
      for (unsigned int itree=0; itree<rest.tree.size(); itree++) {
	 coot::atom_tree_t::atom_tree_index_t atom_id_index = name_to_index[rest.tree[itree].atom_id];
	 if (atom_id_index.is_assigned()) {
	    retval = 1;
	    int idx = atom_id_index.index();
	    coot::atom_tree_t::atom_tree_index_t atom_back_index =
	       name_to_index[rest.tree[itree].atom_back];
	    if (rest.tree[itree].atom_back != "") { 
	       // if connect_type is START then the back atom is n/a -> UNASSIGNED
	       if (atom_back_index.is_assigned()) { 
		  if (atom_back_index.index() < n_residue_atoms) { 
		     // there should be only one back atom for each vertex
		     atom_vertex_vec[idx].backward.push_back(atom_back_index.index());
		     // 
		     // make a synthetic forward atom
		     add_unique_forward_atom(atom_back_index.index(), idx);
		  }
	       }
	    }
	       
	    coot::atom_tree_t::atom_tree_index_t atom_forward_index =
	       name_to_index[rest.tree[itree].atom_forward];

	    if (atom_forward_index.is_assigned()) {
	       if (atom_forward_index.index() < n_residue_atoms)
		  add_unique_forward_atom(idx, atom_forward_index.index());
	    }

	    if (0) { 
	       std::cout << "   itree " << itree << " :"
			 << rest.tree[itree].atom_id << ": :"
		      << rest.tree[itree].atom_back << ": :"
			 << rest.tree[itree].atom_forward << ": -> "
			 << atom_id_index.index() << " "
			 << atom_back_index.index() << " "
			 << atom_forward_index.index()
			 << std::endl;
	    }
	    
	    atom_vertex_vec[idx].connection_type = coot::atom_vertex::STANDARD;
	    if (rest.tree[itree].connect_type == "START")
	       atom_vertex_vec[idx].connection_type = coot::atom_vertex::START;
	    if (rest.tree[itree].connect_type == "END")
	       atom_vertex_vec[idx].connection_type = coot::atom_vertex::END;
	 }
      }
   }

   // print out the atom tree
   if (0) {
      std::cout << "debug:: ==== atom_vertex_vec === from fill_atom_vertex_vec " << std::endl;
      for (unsigned int iv=0; iv<atom_vertex_vec.size(); iv++) {
	 std::cout << "   atom_vertex_vec[" << iv << "] forward atom ";
	 for (unsigned int ifo=0; ifo<atom_vertex_vec[iv].forward.size(); ifo++) 
	    std::cout << atom_vertex_vec[iv].forward[ifo] << " ";
	 std::cout << std::endl;
      }
   }
   return retval;
}




// Add forward_atom_index as a forward atom of this_index - but only
// if forward_atom_index is not already a forward atom of this_index.
void
coot::atom_tree_t::add_unique_forward_atom(int this_index, int forward_atom_index) { 

   bool ifound = 0;
   for (unsigned int ifo=0; ifo<atom_vertex_vec[this_index].forward.size(); ifo++) {
      if (atom_vertex_vec[this_index].forward[ifo] == forward_atom_index) {
	 ifound = 1;
	 break;
      }
   }

   std::pair<int, std::vector<coot::atom_tree_t::atom_tree_index_t> >
      forward_atoms_of_forward_atom_index_pair = get_forward_atoms(forward_atom_index);

   if (0) 
      std::cout << " in add_unique_forward_atom get_forward_atoms called "
		<< forward_atoms_of_forward_atom_index_pair.first 
		<< " times - indices size: "
		<< forward_atoms_of_forward_atom_index_pair.second.size() << std::endl;
   
   if (0) { 
      std::cout << "debug:: forward atoms of " << forward_atom_index << ":"; 
      for (unsigned int i=0; i<forward_atoms_of_forward_atom_index_pair.second.size(); i++)
	 std::cout << " " << forward_atoms_of_forward_atom_index_pair.second[i].index();
      std::cout << std::endl;
   }

   for (unsigned int i=0; i<forward_atoms_of_forward_atom_index_pair.second.size(); i++)
      if (this_index == forward_atoms_of_forward_atom_index_pair.second[i].index()) {
	 ifound = 1;
	 if (0) { 
	    std::cout << " reject this attempt to add synthetic forward atom because "
		      << this_index << " is a forward atom of " << forward_atom_index
		      << std::endl;
	 }
      }

   if (ifound == 0) {
//       std::cout << "So, new forward: to this_index " << this_index << " added forward "
// 		<< forward_atom_index << std::endl;
      atom_vertex_vec[this_index].forward.push_back(forward_atom_index);
   } 

} 


// Throw exception on able to rotate atoms.
//
// Return the new torsion angle (use the embedded torsion on index2 if you can)
// Otherwise return -1.0;
// 
double
coot::atom_tree_t::rotate_about(const std::string &atom1, const std::string &atom2,
				double angle, bool reversed_flag) {

   double new_torsion = 0.0;
   
   coot::atom_tree_t::atom_tree_index_t index2 = name_to_index[atom1];
   coot::atom_tree_t::atom_tree_index_t index3 = name_to_index[atom2];

   if (index2.is_assigned()) { 
      if (index3.is_assigned()) {

	 // is index3 in the forward atoms of index2?
	 // if not, then swap and test again.
	 //
	 bool index3_is_forward = 0;
	 for (unsigned int ifo=0; ifo<atom_vertex_vec[index2.index()].forward.size(); ifo++) {
	    if (atom_vertex_vec[index2.index()].forward[ifo] == index3.index()) {
	       index3_is_forward = 1;
	       break;
	    }
	 }

	 if (! index3_is_forward) {
	    // perhaps index2 is the forward atom of index3?
	    bool index2_is_forward = 0;
	    for (unsigned int ifo=0; ifo<atom_vertex_vec[index3.index()].forward.size(); ifo++) {
	       if (atom_vertex_vec[index3.index()].forward[ifo] == index2.index()) {
		  index2_is_forward = 1;
		  break;
	       }
	    }

	    // if index2 *was* the forward atom of index3, then swap
	    // around index2 and 3 into "standard" order.
	    // 
	    if (index2_is_forward) {
	       std::swap(index2, index3);
	       index3_is_forward = 1;
	    }
	 }

	 // OK, try again
	 if (index3_is_forward) {

	    std::pair<int, std::vector<coot::atom_tree_t::atom_tree_index_t> > moving_atom_indices_pair =
	       get_forward_atoms(index3);
	    std::vector<coot::atom_tree_t::atom_tree_index_t> moving_atom_indices =
	       moving_atom_indices_pair.second;
// 	    std::cout << " in rotate_about get_forward_atoms called " << moving_atom_indices_pair.first
// 		      << " times - indices size: " << moving_atom_indices_pair.second.size() << std::endl;

	    // Maybe a synthetic forward atom was made, and later on
	    // in the dictionary, it was a real (normal) forward atom
	    // of a different atom.  In that case there can be 2
	    // copies of the atom in moving_atom_indices.  So let's
	    // filter them one.  Only the first copy.
	    // 
	    std::vector<coot::atom_tree_t::atom_tree_index_t> unique_moving_atom_indices =
	       uniquify_atom_indices(moving_atom_indices);

	    if (0) {
	       std::cout << "number of moving atoms based on atom index " << index3.index()
			 << " is " << moving_atom_indices.size() << std::endl;
	       for (unsigned int imov=0; imov<unique_moving_atom_indices.size(); imov++) {
		  std::cout << "now  moving atom[" << imov << "] is "
			    << unique_moving_atom_indices[imov].index() << std::endl;
	       }
	    }
	    

	    if (reversed_flag) {
	       unique_moving_atom_indices = complementary_indices(unique_moving_atom_indices,
								  index2, index3);
	    } 

	    // so now we have a set of moving atoms:
	    // set up the coordinates for the rotation, and rotate
	    PPCAtom residue_atoms = 0;
	    int n_residue_atoms;
	    residue->GetAtomTable(residue_atoms, n_residue_atoms);
	    CAtom *at2 = residue_atoms[index2.index()];
	    CAtom *at3 = residue_atoms[index3.index()];
	    clipper::Coord_orth base_atom_pos(at2->x, at2->y, at2->z);
	    clipper::Coord_orth    third_atom(at3->x, at3->y, at3->z);;
	    clipper::Coord_orth direction = third_atom - base_atom_pos;
	    rotate_internal(unique_moving_atom_indices, direction, base_atom_pos, angle);

	    // set the new_torsion (return value) if possible.
	    // 
	    if (atom_vertex_vec[index2.index()].torsion_quad.first) {
	       new_torsion = quad_to_torsion(index2);
	    } 
	 }
      }
   }
   return new_torsion;
}


double
coot::atom_tree_t::quad_to_torsion(const coot::atom_tree_t::atom_tree_index_t &index2) const {

   // std::cout << " this bond has a chi angle assigned" << std::endl;
   coot::atom_index_quad quad = atom_vertex_vec[index2.index()].torsion_quad.second;
   clipper::Coord_orth co[4];
   PPCAtom residue_atoms;
   int n_residue_atoms;
   residue->GetAtomTable(residue_atoms, n_residue_atoms);
   co[0] = clipper::Coord_orth(residue_atoms[quad.index1]->x,
			       residue_atoms[quad.index1]->y,
			       residue_atoms[quad.index1]->z);
   co[1] = clipper::Coord_orth(residue_atoms[quad.index2]->x,
			       residue_atoms[quad.index2]->y,
			       residue_atoms[quad.index2]->z);
   co[2] = clipper::Coord_orth(residue_atoms[quad.index3]->x,
			       residue_atoms[quad.index3]->y,
			       residue_atoms[quad.index3]->z);
   co[3] = clipper::Coord_orth(residue_atoms[quad.index4]->x,
			       residue_atoms[quad.index4]->y,
			       residue_atoms[quad.index4]->z);
   double ar = clipper::Coord_orth::torsion(co[0], co[1], co[2], co[3]);
   double new_torsion = clipper::Util::rad2d(ar);
   return new_torsion;
} 

// Back atoms, not very useful - is is only a single path.
// 
// some sort of recursion here
std::vector<coot::atom_tree_t::atom_tree_index_t>
coot::atom_tree_t::get_back_atoms(const coot::atom_tree_t::atom_tree_index_t &index) const {

   std::vector<coot::atom_tree_t::atom_tree_index_t> v;

   if (atom_vertex_vec[index.index()].connection_type == coot::atom_vertex::END)
      return v;
   
   // all of these
   for (unsigned int ib=0; ib<atom_vertex_vec[index.index()].backward.size(); ib++) {
      v.push_back(atom_vertex_vec[index.index()].backward[ib]);
   }

   // and the back atoms of back atoms.
   for (unsigned int ib=0; ib<atom_vertex_vec[index.index()].backward.size(); ib++) {
      int index_back = atom_vertex_vec[index.index()].backward[ib];
      coot::atom_tree_t::atom_tree_index_t back_index(index_back); // ho ho
      std::vector<coot::atom_tree_t::atom_tree_index_t> nv = 
	 coot::atom_tree_t::get_back_atoms(index_back);
      for (unsigned int inv=0; inv<nv.size(); inv++)
	 v.push_back(nv[inv]);
   }
   return v;
}


// with forward recursion
// 
std::pair<int, std::vector<coot::atom_tree_t::atom_tree_index_t> > 
coot::atom_tree_t::get_forward_atoms(const coot::atom_tree_t::atom_tree_index_t &index) const {

   std::vector<coot::atom_tree_t::atom_tree_index_t> v;
   int n_forward_count = 1; // this time at least

   // If atom_vertex_vec was not filled, then we should not index into
   // it with index.index():  Stops a crash, at least.
   // 
   if (index.index() >= atom_vertex_vec.size())
      return std::pair<int, std::vector<coot::atom_tree_t::atom_tree_index_t> > (n_forward_count, v);


   // how can this happen?
   if (atom_vertex_vec[index.index()].connection_type == coot::atom_vertex::START) {
      // std::cout << " in get_forward_atoms index " <<  index.index() << " at start" << std::endl;
      return std::pair<int, std::vector<coot::atom_tree_t::atom_tree_index_t> > (n_forward_count, v);
   } 

   // Bizarre atom trees in the refmac dictionary - the END atom is
   // marks the C atom - not the O.  How can that be sane?
   // 
//    if (atom_vertex_vec[index.index()].connection_type == coot::atom_vertex::END) { 
//       // std::cout << " in get_forward_atoms index " <<  index.index() << " at end" << std::endl;
//       return v;
//    } 
   
   // all of these
   // 
   if (0) {
      std::cout << "DEBUG::get_forward_atoms of " << index.index() << " has "
		<< atom_vertex_vec[index.index()].forward.size()
		<< " forward atoms ";
      for (unsigned int ifo=0; ifo<atom_vertex_vec[index.index()].forward.size(); ifo++)
	 std::cout << " " << atom_vertex_vec[index.index()].forward[ifo];
      std::cout << std::endl;
   }
   
   for (unsigned int ifo=0; ifo<atom_vertex_vec[index.index()].forward.size(); ifo++) { 
      v.push_back(atom_vertex_vec[index.index()].forward[ifo]);
      if (0) 
	 std::cout << " adding to forward atoms " << atom_vertex_vec[index.index()].forward[ifo]
		   << " which is atom_vertex_vec[" << index.index() << "].forward[" << ifo << "]"
		   << std::endl;
   }

   // and the forward atoms of the forward atoms
   for (unsigned int ifo=0; ifo<atom_vertex_vec[index.index()].forward.size(); ifo++) {
      int index_forward = atom_vertex_vec[index.index()].forward[ifo];
      coot::atom_tree_t::atom_tree_index_t forward_index(index_forward);
      std::pair<int, std::vector<coot::atom_tree_t::atom_tree_index_t> > nv_pair = 
	 coot::atom_tree_t::get_forward_atoms(forward_index);
      std::vector<coot::atom_tree_t::atom_tree_index_t> nv = nv_pair.second;
      n_forward_count += nv_pair.first;
      for (unsigned int inv=0; inv<nv.size(); inv++) {
	 if (0) 
	    std::cout << " adding item  " << inv << " which has index " << nv[inv].index()
		      << " as a forward atom of " << index.index() << std::endl;
	 v.push_back(nv[inv]);
      }
   }

   return std::pair<int, std::vector<coot::atom_tree_t::atom_tree_index_t> > (n_forward_count, uniquify_atom_indices(v));
}



std::vector<coot::atom_tree_t::atom_tree_index_t>
coot::atom_tree_t::uniquify_atom_indices(const std::vector<coot::atom_tree_t::atom_tree_index_t> &vin) const {

   std::vector<coot::atom_tree_t::atom_tree_index_t> v;

   for (unsigned int iv=0; iv<vin.size(); iv++) {
      // vin[iv] is in v?
      bool found = 0;

      for (unsigned int ii=0; ii<v.size(); ii++) {
	 if (vin[iv] == v[ii]) {
	    found = 1;
	    break;
	 } 
      }

      if (found == 0)
	 v.push_back(vin[iv]);
   } 
   return v; 
} 



// Return the complementary indices c.f. the moving atom set, but do
// not include index2 or index3 in the returned set (they do not move
// even with the reverse flag (of course)).
std::vector<coot::atom_tree_t::atom_tree_index_t> 
coot::atom_tree_t::complementary_indices(const std::vector<coot::atom_tree_t::atom_tree_index_t> &moving_atom_indices,
					 const coot::atom_tree_t::atom_tree_index_t &index2,
					 const coot::atom_tree_t::atom_tree_index_t &index3) const {

   std::vector<coot::atom_tree_t::atom_tree_index_t> v;
   for (unsigned int ivert=0; ivert<atom_vertex_vec.size(); ivert++) {
      bool ifound = 0;
      for (unsigned int im=0; im<moving_atom_indices.size(); im++) {
	 if (moving_atom_indices[im].index() == ivert) {
	    ifound = 1;
	    break;
	 }
      }
      if (ifound == 0)
	 if (index2.index() != ivert)
	    if (index3.index() != ivert)
	       v.push_back(ivert);
   }
   
   return v;
} 




// so now we have a set of moving and non-moving atoms:
//
// the angle is in radians. 
void
coot::atom_tree_t::rotate_internal(std::vector<coot::atom_tree_t::atom_tree_index_t> moving_atom_indices,
				   const clipper::Coord_orth &dir,
				   const clipper::Coord_orth &base_atom_pos,
				   double angle) {

//    std::cout << "in rotate_internal with " << moving_atom_indices.size() << " moving atoms "
// 	     << std::endl;
   PPCAtom residue_atoms = 0;
   int n_residue_atoms;
   residue->GetAtomTable(residue_atoms, n_residue_atoms);
   
   for (unsigned int im=0; im<moving_atom_indices.size(); im++) {
      int idx = moving_atom_indices[im].index();
      CAtom *at = residue_atoms[idx];
      clipper::Coord_orth po(at->x, at->y, at->z);
      clipper::Coord_orth pt = coot::util::rotate_round_vector(dir, po, base_atom_pos, angle);
      if (0)
	 std::cout << " rotate_internal() moving atom number " << im << " " << at->name
		   << " from " << at->x << "," << at->y << "," << at->z << " "
		   << pt.format() << std::endl;
      at->x = pt.x(); 
      at->y = pt.y(); 
      at->z = pt.z(); 
   } 
}

// can throw an exception
double
coot::atom_tree_t::set_dihedral(const std::string &atom1, const std::string &atom2,
				const std::string &atom3, const std::string &atom4,
				double angle) {

   double dihedral_angle = 0.0;
   coot::atom_tree_t::atom_tree_index_t i1 = name_to_index[atom1];
   coot::atom_tree_t::atom_tree_index_t i2 = name_to_index[atom2];
   coot::atom_tree_t::atom_tree_index_t i3 = name_to_index[atom3];
   coot::atom_tree_t::atom_tree_index_t i4 = name_to_index[atom4];

   if (i1.is_assigned() && i2.is_assigned() && i3.is_assigned() && i4.is_assigned()) {
      if (atom_vertex_vec[i2.index()].torsion_quad.first) { 
	 double current_dihedral_angle = quad_to_torsion(i2);
	 double diff = angle - current_dihedral_angle;
	 if (diff > 360.0)
	    diff -= 360.0;
	 if (diff < -360.0)
	    diff += 360.0;
	 rotate_about(atom2, atom3, clipper::Util::d2rad(diff), 0);
	 dihedral_angle = quad_to_torsion(i2);
      } else {
	 std::string mess = "Torsion for ";
	 mess += atom2;
	 mess += " to ";
	 mess += atom3;
	 mess += " not found ";
	 throw std::runtime_error(mess);
      } 
   } else {
      std::string mess = "Atom name(s) not found in residue ";
      throw std::runtime_error(mess);
   }
   return dihedral_angle;
}

