/* high-res/high-res.cc
 * 
 * Copyright 2003, 2004  by Paul Emsley, The University of York
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
#include <map>
#include "coords/mmdb-extras.h" // for atom_selection_container_t used in coot-close
#include "coords/mmdb.h" // for formatting of mmdb::Atom
#include "coords/coot-close.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-atom-graph.hh"
#include "mini-mol/mini-mol-utils.hh"
#include "high-res.hh"

coot::high_res::high_res(const coot::minimol::molecule &mol) {

   // first get the atom at the "middle" of the mol.  Let's get the
   // "most connected" atom (i.e. the atom with the most contacts with
   // a radius of 15A).
   // 
   // Also return the mmdb_manager that we used to find this.
   // 
   std::pair<clipper::Coord_orth, mmdb::Manager *> middle_pos = get_middle_pos(mol);
   clipper::Coord_orth target_pos = middle_pos.first;

   fill_globular_protein(mol, middle_pos.first, middle_pos.second);
   // clear up middle_pos.second here

   make_trees();
}

coot::high_res::high_res(const coot::minimol::molecule &mol,
			 int iflag) {

   std::pair<clipper::Coord_orth, mmdb::Manager *> middle_pos = get_middle_pos(mol);
   fill_globular_protein_by_fragments(mol, middle_pos.first, middle_pos.second);
   delete middle_pos.second;
}

coot::high_res::high_res(const coot::minimol::molecule &mol,
			 const clipper::Coord_orth &given_centre) {

   std::pair<clipper::Coord_orth, mmdb::Manager *> middle_pos = get_middle_pos(mol);
   fill_globular_protein_by_fragments(mol, given_centre, middle_pos.second);
   delete middle_pos.second;

}


void
coot::high_res::output_pdb(const std::string &filename) const { 
   globular_molecule.write_file(filename, 20.0);
}


std::pair<clipper::Coord_orth, mmdb::Manager *> 
coot::high_res::get_middle_pos(const coot::minimol::molecule &minimol_mol) const { 

   std::pair<clipper::Coord_orth, mmdb::Manager *> r;
   long i_contact_group = 1;

   mmdb::Manager *mol = minimol_mol.pcmmdbmanager();
   r.second = mol;
   mmdb::Contact *pscontact = NULL;
   int n_contacts = -1;
   float min_dist = 1.0;
   float max_dist = 15.0;
   mmdb::mat44 my_matt;
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
   atom_selection_container_t asc = make_asc(mol);
   
   int err = mol->GetTMatrix(my_matt, 0, 0, 0, 0);
   if (err != mmdb::SYMOP_Ok) {
      std::cout << "!! Warning:: No symmetry available for this molecule"
		<< std::endl;
   } else { 

      mol->SeekContacts(asc.atom_selection, asc.n_selected_atoms, 
			asc.atom_selection, asc.n_selected_atoms,
			min_dist, max_dist, // min, max distances
			0,        // seqDist 0 -> in same res also
			pscontact, n_contacts,
			0, &my_matt, i_contact_group);

      std::cout << "INFO:: There were " << n_contacts << " contacts. " << std::endl;

      std::vector<int> contact_atom(asc.n_selected_atoms);
      for (int i=0; i<asc.n_selected_atoms; i++)
	 contact_atom[i] = 0;
   
      for (int i=0; i<n_contacts; i++) { 
	 contact_atom[pscontact[i].id1]++;;
      } 
      int most_contacts_index = -1;
      int most_contacts = -1;
      for (int i=0; i<asc.n_selected_atoms; i++) { 
	 if (contact_atom[i] > most_contacts) { 
	    most_contacts_index = i;
	    most_contacts = contact_atom[i];
	 }
      }
      if (most_contacts >= 0) { 
	 mmdb::Atom *at = asc.atom_selection[most_contacts_index];
	 r.first = clipper::Coord_orth(at->x, at->y, at->z);
      } 
      delete [] pscontact;
      std::cout << "INFO:: get_middle_pos: returns " << r.first.format()  << " with " 
		<< asc.n_selected_atoms << " atoms " << std::endl;
   }
   return r;
}

// needed?
#include "coords/coot-close.hh"

void
coot::high_res::fill_globular_protein(const coot::minimol::molecule &mol, 
				      const clipper::Coord_orth &target_pos_in,
				      mmdb::Manager *mmdb_mol) { 

   clipper::Coord_orth sum_atoms = target_pos_in;
   clipper::Coord_orth target_pos = target_pos_in;
   // now a bit of jiggery pokery
   // Lets add 20 atoms at the target_pos;
   sum_atoms = clipper::Coord_orth(20.0*target_pos.x(), 
				   20.0*target_pos.y(), 
				   20.0*target_pos.z());
   double n_atoms = 20;

   clipper::Coord_orth t;
   globular_molecule.set_cell_symm(mol);
   for (unsigned int ifrag=0; ifrag<mol.fragments.size(); ifrag++) {
      int igfrag = globular_molecule.fragment_for_chain(mol[ifrag].fragment_id);
      coot::minimol::residue residue(1, // 1 is better for bonding 
  				        // rather than mol[ifrag][ires].seqnum 
				     "ALA");
      for (int ires=mol[ifrag].min_res_no(); ires<=mol[ifrag].max_residue_number();
	   ires++) {
	 for (unsigned int iat=0; iat<mol[ifrag][ires].n_atoms(); iat++) {

	    t = ::closest_approach(mol[ifrag][ires][iat].pos, 
				   target_pos, mmdb_mol);
	    
// 	    std::cout << "MOVING start pos " 
// 		      << mol[ifrag][ires][iat].pos.x() << " " 
// 		      << mol[ifrag][ires][iat].pos.y() << " " 
// 		      << mol[ifrag][ires][iat].pos.z() << "\n";

	    sum_atoms += t;
	    n_atoms += 1.0;

	    residue.addatom(" C  ", " C", t, "", 1.0, 30.0);
	 }
      }
      try { 
	 globular_molecule[igfrag].addresidue(residue, 0);
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "ERROR:: fill_globular_protein() " << rte.what() << std::endl;
      } 
   }
   std::cout << "DEBUG:: ##################### globular_molecule created"
	     << std::endl;

}


void
coot::high_res::fill_globular_protein_by_fragments(const coot::minimol::molecule &mol, 
						   const clipper::Coord_orth &target_pos_in,
						   mmdb::Manager *mmdb_mol) {
   
   // first we want a fragmented minimol molecule
   globular_molecule = mol.fragmentize();
   globular_molecule.set_cell_symm(mol);

   clipper::Coord_orth target_pos = target_pos_in;
   std::cout << "Target centre of protein " << target_pos.format() << std::endl;

   // Here we move fragments
   // 
   for (unsigned int ifrag=0; ifrag<globular_molecule.fragments.size(); ifrag++) {
      clipper::Coord_orth p = globular_molecule[ifrag].midpoint();
//       std::cout << "   (" << p.x() << " " << p.y() << " " << p.x()
// 		<< ") ; fragment mid point"<< std::endl;
      clipper::RTop_orth rtop = closest_approach_transformation(p, target_pos, mmdb_mol);
      globular_molecule[ifrag].transform(rtop);
   }

}



void
coot::high_res::make_trees() { 

   mmdb::Manager *mol = globular_molecule.pcmmdbmanager();
   atom_selection_container_t asc = make_asc(mol);

   mmdb::Contact *pscontact = NULL;
   int n_contacts;
   float min_dist = 1.0;
   float max_dist = 1.9;
   long i_contact_group = 1;
   std::vector<std::vector<node_info> > contact_indices(asc.n_selected_atoms);
   int nsymops = 0;

   if (mol->GetNumberOfSymOps() == 0) {
      std::cout << "WARNING:: no symmetry available! no connections!\n";
   } else { 
      mmdb::realtype cell[6];
      mmdb::realtype vol;
      int orthcode;
      mol->GetCell(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5], vol, orthcode);
      std::cout << "  Spacegroup: " << mol->GetSpaceGroup() << "\n";
      std::cout << "        Cell: " 
		<< cell[0] << " " << cell[1] << " " << cell[2] << " " 
		<< cell[3] << " " << cell[4] << " " << cell[5] << " " << "\n";
      std::cout << "Here there are " << mol->GetNumberOfSymOps()
		<< " symmetry operators\n";
      std::cout << "Finding contacts for " << asc.n_selected_atoms << " atoms\n";
      

      mmdb::mat44 my_matt;
      for (int i=0; i<4; i++) 
	 for (int j=0; j<4; j++) 
	    my_matt[i][j] = 0.0;
      for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

      for (int ix = -1; ix < 2; ix++) { 
	 for (int iy = -1; iy < 2; iy++) { 
	    for (int iz = -1; iz < 2; iz++) {

	       for (int isym = 0; isym < mol->GetNumberOfSymOps(); isym++) { 
	       
// 		  std::cout << "Round: " <<  isym << " " 
// 			    <<  ix << " " << iy << " " << iz << std::endl;
		  int err = mol->GetTMatrix(my_matt, isym, ix, iy, iz); 
		  nsymops++;
		  if (err != 0) { 
		     std::cout << "WARNING:: something BAD with mmdb::CMMDBCryst.GetTMatrix\n";
		  } else { 
		     pscontact = NULL;
		     mol->SeekContacts(asc.atom_selection, asc.n_selected_atoms, 
				       asc.atom_selection, asc.n_selected_atoms,
				       min_dist, max_dist, // min, max distances
				       0,        // seqDist 0 -> in same res also
				       pscontact, n_contacts,
				       0, &my_matt, i_contact_group);

// 		     for (int i=0; i<4; i++) 
// 			for (int j=0; j<4; j++) 
// 			   if (fabs(my_matt[i][j]) < 0.001)
// 			      my_matt[i][j] = 0.0;
// 		     clipper::Mat33<double> clipper_mat(my_matt[0][0], my_matt[0][1], my_matt[0][2],
// 							my_matt[1][0], my_matt[1][1], my_matt[1][2],
// 							my_matt[2][0], my_matt[2][1], my_matt[2][2]);
// 		     clipper::Coord_orth  cco(my_matt[0][3], my_matt[1][3], my_matt[2][3]);
// 		     clipper::RTop_orth rtop(clipper_mat, cco);
// 		     std::cout << "rtop: \n" << rtop.format() << "\n"; 
		      
		  
		     // Now we need to do tree bond indexing like we do for mgtree elsewhere.
		     // 
		     int in1, in2;
		     for (int i=0; i<n_contacts; i++) { 
			in1 = pscontact[i].id1;
			in2 = pscontact[i].id2;
			if (isym == 0 && ix == 0 && iy == 0 && iz == 0) {
			   // don't need symmetry:
			   contact_indices[in1].push_back(coot::node_info(in2));
			} else { 
// 			   std::cout << "got a contact: " << in1 << " to " 
// 				     << in2 << " via " 
// 				     << isym << " " <<  ix << " " << iy << " " << iz << "\n";
			   contact_indices[in1].push_back(coot::node_info(mol, in2, isym, ix, iy, iz));
			}
		     }
		     delete [] pscontact;
		  }
	       }
	    }
	 }
      }

      std::vector<clipper::Coord_orth> coords;
      for (int iat=0; iat<asc.n_selected_atoms; iat++) 
	 coords.push_back(clipper::Coord_orth(asc.atom_selection[iat]->x,
					      asc.atom_selection[iat]->y,
					      asc.atom_selection[iat]->z));
   
      //    std::cout << "There are " << contact_indices.size()
      // 	     << " atoms "<< std::endl;
      //    for(int i=0; i<contact_indices.size(); i++) {
      //       for(int j=0; j<contact_indices[i].size(); j++) {
      // 	 std::cout << contact_indices[i][j] << " ";
      //       }
      //       std::cout << "\n";
      //    }


      // the output coordinates of atom_graph should contain
      // symmetry/cell info (in the pdb file, I mean).
      coot::atom_graph ag(mol, contact_indices, coords);
      ag.sort();

   }
}


void
coot::high_res::buccafilter_neighbours() {

   // Tinker with globular_molecule
   //

   // First find the neighbours.  We need an mmdb molecule and an atom
   // selection to do that.
   //

   mmdb::Manager *mol = globular_molecule.pcmmdbmanager();
   atom_selection_container_t asc = make_asc(mol);
   
   mmdb::Contact *contact = NULL;
   int ncontacts;

   mmdb::realtype min_dist = 0.0;
   mmdb::realtype max_dist = 2.0;

   std::cout << "INFO asc has has " << asc.n_selected_atoms
	     << " selected atoms" << std::endl;

   asc.mol->SeekContacts(asc.atom_selection, asc.n_selected_atoms,
			 asc.atom_selection, asc.n_selected_atoms,
			 min_dist, max_dist, 0,  // seqDist
			 contact, ncontacts);

   for (int i=0; i<10; i++)
      std::cout << "  " << asc.atom_selection[i] << std::endl;
   
   std::cout << "found " << ncontacts << " contacts" << std::endl;
   if (ncontacts > 0) {

      // mark each atom with a group (cluster) id.  Initially set the
      // group id to -1 (unset).
      int uddhandle = asc.mol->RegisterUDInteger(mmdb::UDR_ATOM,
						 "buccaneer filter group");
      for (int i=0; i<asc.n_selected_atoms; i++) {
	 int istat = asc.atom_selection[i]->PutUDData(uddhandle, -1);
	 if (istat == mmdb::UDDATA_WrongUDRType) {
	    std::cout << "ERROR::  mmdb:UDDATA_WrongUDRType in "
		      << "buccafilter" << std::endl;
	 }
      }

      // OK now every atom has a group marker - currently unset.
      //
      // now something tricky.  We need to find neighbours of this
      // atom that has the same atom name as this atom... and
      // neighbours of that neighbours until there are no more
      // neighbours of neighbours that are unmarked.
      //

      // neighbours is a vector: for each atom we have a vector that
      // is a list of the atoms to which this atom has neighbours
      // 
      std::vector<std::vector<int> > neighbours(asc.n_selected_atoms);
      for (int ic=0; ic< ncontacts; ic++)
	 neighbours[ contact[ic].id1 ].push_back(contact[ic].id2);

      int igroup = 0;
      int ierr;
      int ic;
      for (int iat=0; iat<asc.n_selected_atoms; iat++) {
	 // Tinker with the UDD of the atom selection
	 ierr = asc.atom_selection[iat]->GetUDData(uddhandle, ic);
	 if (ierr == mmdb::UDDATA_Ok) {
	    if (ic == -1) {
	       char *at_name = asc.atom_selection[iat]->name;
	       std::string atom_name(at_name);
	       mark_neighbours(iat, igroup, at_name, neighbours,
			       asc.atom_selection, uddhandle);
	       igroup++;
	    }
	 }
      }
      std::cout << " Found " << igroup << " groups " << std::endl;
      // 
      //
      // So every atom is marked now.
      // Let's make the groups.
      //
      std::vector<std::vector<int> > groups(igroup);
      for (int iat=0; iat<asc.n_selected_atoms; iat++) {
	 asc.atom_selection[iat]->GetUDData(uddhandle, igroup);
	 groups[igroup].push_back(iat);
      }
      // replace globular_molecule with the result of the filtering
      globular_molecule = filter_on_groups(groups, asc.mol, 
 					   asc.atom_selection,
 					   asc.n_selected_atoms);
   } 

   delete mol; // give back the pcmmdbmanager memory
   std::cout << "Done buccafilter filtering" << std::endl;
}


void
coot::high_res::buccafilter() {

   // Change globular_molecule
   //
   // First Lets make a molecule that has fragments that are sorted by
   // length.

   std::vector<coot::minimol::fragment> fragments = globular_molecule.fragments;

   std::sort(fragments.begin(), fragments.end(), fragment_sorter);

//    DEBUGGING
//    for (int ifrag=0; ifrag<fragments.size(); ifrag++) {
//       int a_min_res_no = 0;
//       for (int ir=0; ir<fragments[ifrag].residues.size(); ir++) {
// 	 if (fragments[ifrag].residues[ir].atoms.size() > 0) {
// 	    a_min_res_no = ir;
// 	    break;
// 	 }
//       }
//       std::cout << "  fragment " << ifrag << " has "
// 		<< fragments[ifrag].residues.size()
// 		<< " residues with max_residue_number: "
// 		<< fragments[ifrag].max_residue_number()
// 		<< " range length " << fragments[ifrag].max_residue_number() - a_min_res_no
// 		<< std::endl;
//    }

   // Actually, we don't need to sort by length, we can use the order
   // that the fragments came in as the sorted order - Kevin says that
   // they are sorted by density fit (likelihood, perhaps).
   // 
   // globular_molecule.fragments = fragments;
   // globular_molecule.write_file("fragmented.pdb"); // debug

   mmdb::Manager *mol = globular_molecule.pcmmdbmanager();
   atom_selection_container_t asc = make_asc(mol);
   
   mmdb::Contact *contact = NULL;
   int ncontacts;

   mmdb::realtype min_dist = 0.0;
   mmdb::realtype max_dist = 2.0;

   std::cout << "INFO asc has " << asc.n_selected_atoms
	     << " selected atoms" << std::endl;

   asc.mol->SeekContacts(asc.atom_selection, asc.n_selected_atoms,
			 asc.atom_selection, asc.n_selected_atoms,
			 min_dist, max_dist, 0,  // seqDist
			 contact, ncontacts);

   std::cout << "found " << ncontacts << " contacts" << std::endl;

   if (ncontacts > 0) {

      for (int i=0; i<10; i++)
	 std::cout << "  " << asc.atom_selection[i] << std::endl;

      std::map<mmdb::Chain *, int> chain_numbers;
      int imod = 1;
      mmdb::Model *model_p = asc.mol->GetModel(imod);
      mmdb::Chain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      std::cout << "mol has " << nchains << " chains " << std::endl;
      for (int ich=0; ich<nchains; ich++) {
	 chain_p = model_p->GetChain(ich);
	 chain_numbers[chain_p] = ich;
      }
   
      int ich1;
      int ich2;
      int uddhandle = asc.mol->RegisterUDInteger(mmdb::UDR_ATOM, "delete me flag");
      for (int iat=0; iat<asc.n_selected_atoms; iat++)
	 asc.atom_selection[iat]->PutUDData(uddhandle, 0);
   
      for (int ic=0; ic< ncontacts; ic++) {
	 mmdb::Chain *chain1 = asc.atom_selection[contact[ic].id1]->GetChain();
	 mmdb::Chain *chain2 = asc.atom_selection[contact[ic].id2]->GetChain();
	 if (chain1 != chain2) {
	    ich1 = chain_numbers[chain1];
	    ich2 = chain_numbers[chain2];
	    if (ich2 > ich1) {
	       // delete
	       asc.atom_selection[contact[ic].id2]->PutUDData(uddhandle, 1);
	    }
	 }
      }
      int idelete = 0;
      for (int iat=0; iat<asc.n_selected_atoms; iat++) {
	 asc.atom_selection[iat]->GetUDData(uddhandle, idelete);
	 if (idelete == 1) {
	    delete asc.atom_selection[iat];
	    asc.atom_selection[iat] = 0;
	 }
      }
      asc.mol->FinishStructEdit();
      asc.mol->WritePDBASCII("asc.mol.pdb");
   }

   coot::minimol::molecule mmol(asc.mol);
   mmol.write_file("mmol.pdb", 20.0);
   globular_molecule = mmol;
   delete asc.mol; // give back the memory
}


void
coot::high_res::add_os() {

   for (unsigned int ifrag=0; ifrag<globular_molecule.fragments.size(); ifrag++) {
      for (int ires = globular_molecule[ifrag].min_res_no();
	   ires<=(globular_molecule[ifrag].max_residue_number()-1); ires++) {
	 if (globular_molecule[ifrag][ires].atoms.size() > 0) {
	    if (globular_molecule[ifrag][ires].atoms.size() > 0) {

	       // test for non-existance of O
	       short int o_exists_in_residue = 0;
	       for (unsigned int iat=0; iat<globular_molecule[ifrag][ires].atoms.size();
		    iat++) {
		  if (globular_molecule[ifrag][ires].atoms[iat].name == " O  ") {
		     o_exists_in_residue = 1;
		     break;
		  }
	       }

	       if (o_exists_in_residue == 0) { 
		  std::pair<short int, clipper::Coord_orth> p =
		     coot::o_position(globular_molecule[ifrag][ires],
				      globular_molecule[ifrag][ires+1]);
		  //  	       std::cout << "::::::::: o p returns " << p.first << " "
		  //  			 << p.second.format() << std::endl;
		  if (p.first) {
		     coot::minimol::atom o(" O  ", " O", p.second, "", 30.0);
		     globular_molecule[ifrag][ires].addatom(o);
		  }
	       }
	    }
	 }
      }
   }
}


void
coot::high_res::add_cbetas() {

   // Let's use a coot-utils function to find the cb position, since
   // we will add cbetas in more than "high-res" models.
   //

   for (unsigned int ifrag=0; ifrag<globular_molecule.fragments.size(); ifrag++) {
      for (int ires = globular_molecule[ifrag].min_res_no();
	   ires<=globular_molecule[ifrag].max_residue_number(); ires++) {
	 if (globular_molecule[ifrag][ires].atoms.size() > 0) {
	    
	    // test for non-existance of CBETA
	    short int cbeta_exists_in_residue = 0;
	    for (unsigned int iat=0; iat<globular_molecule[ifrag][ires].atoms.size();
		 iat++) {
	       if (globular_molecule[ifrag][ires].atoms[iat].name == " CB ") {
		  cbeta_exists_in_residue = 1;
		  break;
	       }
	    }

	    if (cbeta_exists_in_residue == 0) { 
	       std::pair<short int, clipper::Coord_orth> p =
		  coot::cbeta_position(globular_molecule[ifrag][ires]);
	       if (p.first == 1)
		  globular_molecule[ifrag][ires].addatom(" CB ", " C", p.second, "", 1.0, 30.0);
	    }
	 }
      }
   }

}



// static
bool
coot::high_res::fragment_sorter(const coot::minimol::fragment &a,
				const coot::minimol::fragment &b) {

   int a_min_res_no = 0;
   int b_min_res_no = 0;

   for (unsigned int ir=0; ir<a.residues.size(); ir++) {
      if (a.residues[ir].atoms.size() > 0) {
	 a_min_res_no = ir;
	 break;
      }
   }
   for (unsigned int ir=0; ir<b.residues.size(); ir++) {
      if (b.residues[ir].atoms.size() > 0) {
	 b_min_res_no = ir;
	 break;
      }
   }
      
   
   return (a.max_residue_number() - a_min_res_no) > (b.max_residue_number() - b_min_res_no);
}


// Tinker with the UDD of the atom selection
// 
void
coot::high_res::mark_neighbours(int iatom, int igroup,
				const std::string &atom_name,
				const std::vector<std::vector<int> > &neighbours,
				mmdb::PPAtom atom_selection, int uddhandle) {

   int ig;
   atom_selection[iatom]->GetUDData(uddhandle, ig);
   if (ig == -1) {
      std::string this_atom_name = atom_selection[iatom]->name;
      if (this_atom_name == atom_name) { 
	 atom_selection[iatom]->PutUDData(uddhandle, igroup);
	 std::vector<int> n = neighbours[iatom];
	 int ns = n.size();
	 for (int in=0; in<ns; in++)
	    mark_neighbours(n[in], igroup, atom_name,
			    neighbours, atom_selection, uddhandle);
//       } else {
// 	 std::cout << "not marking " << this_atom_name << " != " << atom_name << std::endl;
      }
   }
}

coot::minimol::molecule
coot::high_res::filter_on_groups(const std::vector<std::vector<int> > &groups,
				 mmdb::Manager *mol,
				 mmdb::PPAtom atom_selection,
				 int n_selected_atoms) const {

   coot::minimol::molecule m;
   for (unsigned int igroup=0; igroup<groups.size(); igroup++) {

      std::cout << "group " << igroup << " has " << groups[igroup].size()
		<< " members" << std::endl;

	 
      // Find the average position in the group
      //
      // Find all the atoms that are within 2A of the average
      // position.  Use these atoms to calculate new mean.  Then
      // average the points that are within 1A of this new mean and
      // reject (ignore) the others.
      //
      // Let the atom spec of this new position be the atom spec of
      // the first atom in group.

      clipper::Coord_orth sum(0.0, 0.0, 0.0);
      int n_grp_ats = groups[igroup].size();
      for (int iat=0; iat<n_grp_ats; iat++) {
	 mmdb::Atom *at = atom_selection[groups[igroup][iat]];
	 clipper::Coord_orth pt(at->x, at->y, at->z);
	 sum += pt;
      }
      double frac = 1.0/double(n_grp_ats);
      clipper::Coord_orth new_sum(0.0, 0.0, 0.0);
      clipper::Coord_orth mean_pt(sum.x()*frac, sum.y()*frac, sum.z()*frac);
      // double dist_crit = 0.5; // A
      int n_group_within_lim = 0;
      for (int iat=0; iat<n_grp_ats; iat++) {
	 mmdb::Atom *at = atom_selection[groups[igroup][iat]];
	 clipper::Coord_orth pt(at->x, at->y, at->z);
	 //      double d = clipper::Coord_orth::length(mean_pt, pt);
	 // 	 if (d < dist_crit) {
	    n_group_within_lim++;
	    new_sum += pt;
	    // }
      }
      if (groups[igroup].size() > 10) { 
	 for (unsigned int iat=0; iat<groups[igroup].size(); iat++) {
	    std::cout << "   " << atom_selection[groups[igroup][iat]] << std::endl;
	    std::cout << "   (set-rotation-centre "
		      << atom_selection[groups[igroup][iat]]->x << " "
		      << atom_selection[groups[igroup][iat]]->y << " "
		      << atom_selection[groups[igroup][iat]]->z << ")" << std::endl;
	 }
      }
      
      if (n_group_within_lim > 0) { 
	 frac = 1.0/double(n_group_within_lim);
	 clipper::Coord_orth av_pt(new_sum.x()*frac,
				   new_sum.y()*frac,
				   new_sum.z()*frac);
	 mmdb::Atom *speced_at = atom_selection[groups[igroup][0]];
	 std::string atom_name(speced_at->name);
	 std::string atom_element(speced_at->element);
	 int resno = speced_at->GetSeqNum();
	 std::string chain_id(speced_at->GetChainID());
	 int ifrag = m.fragment_for_chain(chain_id);
	 coot::minimol::atom atom(atom_name, atom_element, av_pt, "", 30.0);
	 // std::cout << "DEBUG:: resno is " << resno << std::endl;
	 m[ifrag][resno].name = "ALA";
	 m[ifrag][resno].seqnum = resno;
	 m[ifrag][resno].addatom(atom);
      } else {
	 std::cout << "OOps! No unfiltered atoms left for "
		   << atom_selection[groups[igroup][0]] << std::endl;
      }
   }
   return m;
}



// Tom does not discuss
// 
// How to convert from (asymmetric unit?) map peaks to nodes
// (potential atoms).  i.e. how do we apply symmetry to put the nodes
// next to each other rather than being in the wrong asymmetric unit.
// 
// Multiple conformations
// 
// Longest path 
// 
// stitching together fragments
// 
// stitching together wrong asymmetric unit fragments
//
// sequence assignement
// 
// I do like the distortion score though.
//

