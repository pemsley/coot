/* pyrogen/restraints-boost.cc
 * 
 * Copyright 2011 by the University of Oxford
 * Copyright 2014, 2015 by Medical Research Council
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

#include "Python.h"

#include "compat/coot-sysdep.h"
#include <GraphMol/GraphMol.h>

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <boost/python.hpp>

using namespace boost::python;

#include <ideal/simple-restraint.hh>
#include <coot-utils/coot-coord-utils.hh>
#include <lidia-core/rdkit-interface.hh>

#include "py-restraints.hh"
#include "restraints-private.hh" // for bond-order conversion

#include "mmff-restraints.hh"

#include "cairo-molecule.hh"

namespace coot {

   RDKit::ROMol *regularize(RDKit::ROMol &r);
   RDKit::ROMol *regularize_with_dict(RDKit::ROMol &r,
                                      PyObject *py_restraints,
                                      const std::string &comp_id);
   // This tries to get a residue from the dictionary using the
   // model_Cartn of the dict_atoms.  If that is not available, return
   // a molecule with atoms that do not have positions.
   RDKit::ROMol *rdkit_mol_chem_comp_pdbx(const std::string &chem_comp_dict_file_name,
                                          const std::string &comp_id);
   RDKit::ROMol *hydrogen_transformations(const RDKit::ROMol &r);
   RDKit::ROMol *mogulify(const RDKit::ROMol &r);

   // fiddle with mol
   void delocalize_guanidinos(RDKit::RWMol *mol);

   // compiling/linking problems - give up for now.
   // PyObject *convert_rdkit_mol_to_pyobject(RDKit::ROMol *mol);

}



BOOST_PYTHON_MODULE(pyrogen_boost) {

   // def("convert_rdkit_mol_to_pyobject", coot::convert_rdkit_mol_to_pyobject,
   // return_value_policy<manage_new_object>());

   def("regularize",               coot::regularize,               return_value_policy<manage_new_object>());
   def("regularize_with_dict",     coot::regularize_with_dict,     return_value_policy<manage_new_object>());
   def("rdkit_mol_chem_comp_pdbx", coot::rdkit_mol_chem_comp_pdbx, return_value_policy<manage_new_object>());
   def("hydrogen_transformations", coot::hydrogen_transformations, return_value_policy<manage_new_object>());
   def("mogulify",                 coot::mogulify,                 return_value_policy<manage_new_object>());
   def("mmff_bonds_and_angles",    coot::mmff_bonds_and_angles,    return_value_policy<manage_new_object>());
   // rdkit-like function name
   def("MolFromPDBXr", coot::rdkit_mol_chem_comp_pdbx, return_value_policy<manage_new_object>());
   def("cairo_png_depict_from_mmcif", coot::cairo_png_depict_from_mmcif);
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
   def("cairo_png_depict_to_string",  coot::cairo_png_string_from_mol);
   def("cairo_svg_depict_to_string",  coot::cairo_svg_string_from_mol);
#endif

   class_<coot::mmff_bond_restraint_info_t>("mmff_bond_restraint_info_t")
      .def("get_idx_1",         &coot::mmff_bond_restraint_info_t::get_idx_1)
      .def("get_idx_2",         &coot::mmff_bond_restraint_info_t::get_idx_2)
      .def("get_type",          &coot::mmff_bond_restraint_info_t::get_type)
      .def("get_resting_bond_length", &coot::mmff_bond_restraint_info_t::get_resting_bond_length)
      .def("get_sigma",               &coot::mmff_bond_restraint_info_t::get_sigma)
      ;

   class_<coot::mmff_angle_restraint_info_t>("mmff_angle_restraint_info_t")
      .def("get_idx_1",    &coot::mmff_angle_restraint_info_t::get_idx_1)
      .def("get_idx_2",    &coot::mmff_angle_restraint_info_t::get_idx_2)
      .def("get_idx_3",    &coot::mmff_angle_restraint_info_t::get_idx_3)
      .def("get_resting_angle",  &coot::mmff_angle_restraint_info_t::get_resting_angle)
      .def("get_sigma",          &coot::mmff_angle_restraint_info_t::get_sigma)
      ;


   // established coot class, works with atom names though - so not useful ATM.
   class_<coot::dict_bond_restraint_t>("dict_bond_restraint_t")
      .def("atom_id_1",  &coot::dict_bond_restraint_t::atom_id_1)
      .def("atom_id_2",  &coot::dict_bond_restraint_t::atom_id_2)
      .def("type",       &coot::dict_bond_restraint_t::type)
      .def("value_dist", &coot::dict_bond_restraint_t::value_dist)
      .def("value_esd",  &coot::dict_bond_restraint_t::value_esd)
      ;

   class_<coot::mmff_b_a_restraints_container_t>("mmff_b_a_restraints_container_t")
      .def("bonds_size",  &coot::mmff_b_a_restraints_container_t::bonds_size)
      .def("angles_size", &coot::mmff_b_a_restraints_container_t::angles_size)
      .def("get_bond",    &coot::mmff_b_a_restraints_container_t::get_bond)
      .def("get_angle",   &coot::mmff_b_a_restraints_container_t::get_angle)
      ;
}


RDKit::ROMol*
coot::mogulify(const RDKit::ROMol &mol) {

   RDKit::RWMol rw(mol);
   coot::mogulify_mol(rw);
   RDKit::ROMol *ro = new RDKit::ROMol(rw);
   return ro;
}


RDKit::ROMol *
coot::regularize(RDKit::ROMol &mol_in) {

   RDKit::ROMol *m = new RDKit::ROMol(mol_in);
   return m;

}

RDKit::ROMol *
coot::regularize_with_dict(RDKit::ROMol &mol_in, PyObject *restraints_py, const std::string &comp_id) {

   coot::dictionary_residue_restraints_t dict_restraints = 
      monomer_restraints_from_python(restraints_py);

   RDKit::RWMol *m = new RDKit::RWMol(mol_in);
   mmdb::Residue *residue_p = make_residue(mol_in, 0, comp_id);
   if (! residue_p) {
      std::cout << "WARNING:: bad residue " << std::endl;
   } else {
      // deep copy residue and add to new molecule.
      mmdb::Manager *cmmdbmanager = util::create_mmdbmanager_from_residue(residue_p);
      mmdb::Residue *new_residue_p = coot::util::get_first_residue(cmmdbmanager);
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      new_residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      std::cout << "------------------ simple_refine() called from "
		<< "restraints-boost.cc:regularize_with_dict()" << std::endl;
      simple_refine(new_residue_p, cmmdbmanager, dict_restraints);
      std::cout << "------------------ simple_refine() finished" << std::endl;
      int iconf = 0;
      update_coords(m, iconf, new_residue_p);
   }
   return static_cast<RDKit::ROMol *>(m);
}


RDKit::ROMol *
coot::rdkit_mol_chem_comp_pdbx(const std::string &chem_comp_dict_file_name,
			       const std::string &comp_id) {

   // This function should not get an mmdb::Residue, it should just use
   // rdkit_mol(const dictionary_residue_restraints_t &restraints)
   
   RDKit::ROMol *mol = new RDKit::ROMol;
   coot::protein_geometry geom;
   geom.set_verbose(false);
   int read_number = 0;
   geom.init_refmac_mon_lib(chem_comp_dict_file_name, read_number);
   bool idealized = false;
   idealized = true; // 20150622 - so that we pick up the coords of OXT in 01Y
   bool try_autoload_if_needed = false;
   int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;

   // We use this (mmdb::Restraints-using) interface to rdkit_mol()
   // because it allows a code path fall-back when there are no
   // coordinates for the atoms of the molecule
   
   mmdb::Residue *residue_p = geom.get_residue(comp_id, imol_enc, idealized, try_autoload_if_needed);

   int imol = 0; // dummy
   std::pair<bool, dictionary_residue_restraints_t> rest =
      geom.get_monomer_restraints(comp_id, imol);

   if (rest.first) {

      if (residue_p) {

	 // makes a 3d conformer

	 // 20140618: This was false - now try true so that it processes 110 (i.e. no more
	 //           valence error).
	 
	 // 20141115: But if this is true, then 313 fails to find C-OC single or
	 //           double bonds (it should be deloc - it's a carboxylate)
	 //

	 try {
           
            // experimental value - for Friday.
	    // bool undelocalize = false;
            // 
	    bool undelocalize_flag = true;
	    int iconf = 0;

	    RDKit::RWMol mol_rw = coot::rdkit_mol(residue_p, rest.second, "", undelocalize_flag);
	    RDKit::ROMol *m = new RDKit::ROMol(mol_rw);

	    // Let's skip this step here if we can - because the chiral centres should have been
	    // assigned in rdkit_mol() and we don't want to wipe them out here.
	    // 
	    bool clean = false;
	    bool force = true;
	    bool flag_stereo = true;
	    // RDKit::MolOps::assignStereochemistry(*m, clean, force, flag_stereo);

	    // Here test that the propety mmcif_chiral_volume_sign has
	    // been set on atoms of mol_rw and m
	    //
	    if (false) {
	       for (unsigned int iat=0; iat<m->getNumAtoms(); iat++) {

		  std::string name;
		  // std::cout << "DEBUG:: rdkit_mol_chem_comp_pdbx(): testing atom "
		  // << iat << std::endl;

		  const RDKit::Atom *at_p = (*m)[iat];
		  try {
		     at_p->getProp("name", name);
		     std::cout << "   name " << name << std::endl;
		  }
		  catch (const KeyErrorException &err) {
		     std::cout << "   no name " << std::endl;
		  }
		  try {
		     std::string cv;
		     at_p->getProp("mmcif_chiral_volume_sign", cv);
		     std::cout << "m: name " << name << " " << cv << std::endl;
		  }
		  catch (const KeyErrorException &err) {
		     // std::cout << "m: no name or cv " << std::endl;
		  }
		  try {
		     std::string cip;
		     at_p->getProp("_CIPCode", cip);
		     std::cout << name << " cip-code " << cip << std::endl;
		  }
		  catch (const KeyErrorException &err) {
		     std::cout << name << " cip-code - " << std::endl;
		  }
		  try {
		     int cip_rank;
		     at_p->getProp("_CIPRank", cip_rank);
		     std::cout << "m: CIP-rank " << cip_rank << std::endl;
		  }
		  catch (const KeyErrorException &err) {
		     std::cout << "no-error: no _CIPRank " << err.what() << std::endl;
		  }
		  catch (const boost::bad_any_cast &bac) {
		     // Goodness knows why this is thrown... 
		     std::cout << "strange - caught bad_any_cast on _CIPRank get" << std::endl;
		  } 
	       }
	    }

	    // debug.  OK, so the bond orders are undelocalized here.
	    // debug_rdkit_molecule(&mol_rw);

	    // Happy Path return

	    return m;
	 }

	 catch (const std::runtime_error &rte) {
	    std::cout << "ERROR:: " << rte.what() << std::endl;
	 }

      } else {

	 // Are you here unexpectedly?  That's because the input cif dictionary doesn't have
	 // _chem_comp_atom.x{yz} or model_Cartn_x{yz} or pdbx_model_Cartn_x{yz}_ideal for
	 // the given comp_id.
	 // 
	 std::cout << "INFO:: No 3d coords in dictionary : using 2d from dictionary bonds and angles."
		   << std::endl;
	 // makes a 2d conformer

	 // Bleugh.  We deal with mol vs m badly here.  
	 // 
	 RDKit::RWMol mol_rw = coot::rdkit_mol(rest.second);
	 RDKit::MolOps::sanitizeMol(mol_rw);
	 RDKit::ROMol *m = new RDKit::ROMol(mol_rw);
	 bool canon_orient = false;
	 bool clear_confs  = false;
	 RDDepict::compute2DCoords(*m, NULL, canon_orient, clear_confs, 10, 20);
	 delete mol;
	 return m;
      }
   }

   // this is not the normal (happy path) return
   return mol;
}
     
RDKit::ROMol *
coot::hydrogen_transformations(const RDKit::ROMol &mol) {


   // debug_rdkit_molecule(&mol);
   RDKit::RWMol *r = new RDKit::RWMol(mol);

   RDKit::ROMol *query_cooh = RDKit::SmartsToMol("[C^2](=O)O[H]");
   RDKit::ROMol *query_n    = RDKit::SmartsToMol("[N^3;H2]");
   std::vector<RDKit::MatchVectType>  matches_cooh;
   std::vector<RDKit::MatchVectType>  matches_n;
   bool recursionPossible = true;
   bool useChirality = true;
   bool uniquify = true;

   // Note that we get the atom indexing using mol, but we change the
   // molecule r (I hope that the atoms are equivalently indexed).

   // 20210923-PE FIXME
   // int matched_cooh = RDKit::SubstructMatch(mol,*query_cooh,matches_cooh,uniquify,recursionPossible, useChirality);
   // int matched_n    = RDKit::SubstructMatch(mol,*query_n,   matches_n,   uniquify,recursionPossible, useChirality);

   // delete atoms after they have all been identified otherwise the indexing goes haywire.
   // 
   std::vector<const RDKit::Atom *> atoms_to_be_deleted;

   if (true) // this is useful info (at the moment at least)
      std::cout << "Hydrogen_transformations:"
		<< "\n    number of COOH matches: " << matches_cooh.size()
		<< "\n    number of NH2  matches: " << matches_n.size()
		<< "\n";

   if (true) {  // debugging
      for (unsigned int imatch_cooh=0; imatch_cooh<matches_cooh.size(); imatch_cooh++) {
	 std::cout << "INFO:: Removable hydrogen COOH matches: ";
	 for (unsigned int i=0; i<matches_cooh[imatch_cooh].size(); i++) { 
	    std::cout << " [" << matches_cooh[imatch_cooh][i].first 
		      << ": "  << matches_cooh[imatch_cooh][i].second
		      << "]";
	 }
	 std::cout << std::endl;
      }
   }
	 
   for (unsigned int imatch_cooh=0; imatch_cooh<matches_cooh.size(); imatch_cooh++) {
      
      const RDKit::Atom *at_c  = (*r)[matches_cooh[imatch_cooh][0].second];
      const RDKit::Atom *at_o1 = (*r)[matches_cooh[imatch_cooh][1].second];
      const RDKit::Atom *at_o2 = (*r)[matches_cooh[imatch_cooh][2].second];
      const RDKit::Atom *at_h  = (*r)[matches_cooh[imatch_cooh][3].second];

      std::string at_c_name, at_o1_name, at_o2_name, at_h_name;

      at_c->setProp( "type_energy", "C");
      at_o1->setProp("type_energy", "OC");
      at_o2->setProp("type_energy", "OC");

      RDKit::Bond *bond_1 = r->getBondBetweenAtoms(at_c->getIdx(), at_o1->getIdx());
      RDKit::Bond *bond_2 = r->getBondBetweenAtoms(at_c->getIdx(), at_o2->getIdx());
      RDKit::Bond *bond_3 = r->getBondBetweenAtoms(at_h->getIdx(), at_o2->getIdx());

      if (bond_1) {
	 if (bond_2) {
	    
	    // We can't make the bonds deloc in combination with a
	    // setting of the formal charge on the at_o2 because
	    // molecule then fails the Oxygen atom valence test in
	    // sanitization.
	    // 
	    bond_1->setBondType(RDKit::Bond::ONEANDAHALF);
	    bond_2->setBondType(RDKit::Bond::ONEANDAHALF);
	 }
      }

      // at_o2->setFormalCharge(-1); // not if the bonds are deloc/ONEANDAHALF.
      if (bond_3)
	 r->removeBond(at_o2->getIdx(), at_h->getIdx());
      else
	 std::cout << "DEBUG:: hydrogen_transformations(): no bond_3 (O-H bond) found!" << std::endl;
      atoms_to_be_deleted.push_back(at_h);
   }

   for (unsigned int imatch_n=0; imatch_n<matches_n.size(); imatch_n++) {
      unsigned int n_idx = matches_n[imatch_n][0].second;
      RDKit::Atom *at_n  = (*r)[n_idx];
      unsigned int degree = at_n->getDegree();
      at_n->setFormalCharge(+1);

      if (false) {  // debugging
         std::string name;
         at_n->getProp("name", name);
         std::cout << "debug:: N-atom idx " << n_idx << " " << name << " has degree " << degree
                   << std::endl;
      }

      if (degree == 4) {
         // it has its 2 hydrogens already
         at_n->setProp("type_energy", "NT2"); // also set to NT3 by SMARTS match in pyrogen.py
      }

      if (degree == 3) {
         at_n->setProp("type_energy", "NT3"); // also set to NT3 by SMARTS match in pyrogen.py
         // add a hydrogen atom and a bond to the nitrogen.
         //
         RDKit::Atom *new_h_at = new RDKit::Atom(1);
         // we want to find the idx of this added H, so we do that by
         // keeping hold of the pointer (otherwise the atom gets copied
         // on addAtom() and we lose the handle on the new H atom in the
         // molecule).
         bool updateLabel=true;
         bool takeOwnership=true;
         r->addAtom(new_h_at, updateLabel, takeOwnership);
         unsigned int h_idx = new_h_at->getIdx();
         if (h_idx != n_idx) {
	    r->addBond(n_idx, h_idx, RDKit::Bond::SINGLE);
         } else {
	    std::cout << "OOOPs: bad indexing on adding an amine H " << h_idx << std::endl;
         }
      }
   }

#if 0
   for(unsigned int idel=0; idel<atoms_to_be_deleted.size(); idel++)
      r->removeAtom(atoms_to_be_deleted[idel]);
#endif

   std::vector<RDKit::Atom *> non_const_atoms_to_be_deleted(atoms_to_be_deleted.size(), nullptr);
   for(unsigned int idel=0; idel<atoms_to_be_deleted.size(); idel++)
      non_const_atoms_to_be_deleted[idel] = const_cast<RDKit::Atom *>(atoms_to_be_deleted[idel]);
   for(unsigned int idel=0; idel<atoms_to_be_deleted.size(); idel++)
      r->removeAtom(non_const_atoms_to_be_deleted[idel]);

   // Here do_deloc_bonds should be false when we deal with PDBe input
   //
   bool do_deloc_bonds = false;
   remove_phosphate_hydrogens(r, do_deloc_bonds);
   remove_sulphate_hydrogens (r, do_deloc_bonds);

   // debug
   if (false)
       std::cout << "DEBUG:: hydrogen_transformations calling sanitizeMol() " 
                 << std::endl;
   // do we neet to sanitize? Yes, we do because we go on to minimize this molecule
   RDKit::MolOps::sanitizeMol(*r);
   if (false)
       std::cout << "DEBUG:: hydrogen_transformations back from sanitizeMol() " 
                 << std::endl;

   // delocalize_guanidinos(r); // not yet.
   
   RDKit::ROMol *ro_mol = new RDKit::ROMol(*r);
   if (0)
      std::cout << "hydrogen_transformations returns mol: " << RDKit::MolToSmiles(*ro_mol) 
                << std::endl;
   delete r;
   return ro_mol;

}

void
coot::delocalize_guanidinos(RDKit::RWMol *mol) {

   RDKit::ROMol *query = RDKit::SmartsToMol("N[C^2](=N)N");
   std::vector<RDKit::MatchVectType>  matches;
   bool recursionPossible = true;
   bool useChirality = true;
   bool uniquify = true;

   // Note that we get the atom indexing using mol, but we change the
   // molecule r (I hope that the atoms are equivalently indexed).

   // 20210923-PE FIXME
   // int matched = RDKit::SubstructMatch(*mol,*query,matches,uniquify,recursionPossible, useChirality);
   int matched = false;

   if (true) // this is useful info (at the moment at least)
      std::cout << "   delocalize guanidinos matches: " << matches.size() << "\n";

   for (unsigned int imatch=0; imatch<matches.size(); imatch++) {
      std::cout << "INFO:: guanidino hydrogen exchanges matches: ";
      for (unsigned int i=0; i<matches[imatch].size(); i++) { 
	 std::cout << " [" <<  matches[imatch][i].first 
		   << ": "  << matches[imatch][i].second
		   << "]";
      }
      std::cout << std::endl;
   }

   
   for (unsigned int imatch=0; imatch<matches.size(); imatch++) {

      RDKit::Atom *at_c   = (*mol)[matches[imatch][1].second];
      RDKit::Atom *at_n1  = (*mol)[matches[imatch][0].second];
      const RDKit::Atom *at_n2  = (*mol)[matches[imatch][2].second];
      const RDKit::Atom *at_n3  = (*mol)[matches[imatch][3].second];

      std::cout << "INFO:: C hybridisation " << at_c->getHybridization() << std::endl;

      RDKit::Bond *bond_1 = mol->getBondBetweenAtoms(at_c->getIdx(), at_n1->getIdx());
      RDKit::Bond *bond_2 = mol->getBondBetweenAtoms(at_c->getIdx(), at_n2->getIdx());
      RDKit::Bond *bond_3 = mol->getBondBetweenAtoms(at_c->getIdx(), at_n3->getIdx());
      
      std::cout << "INFO:: C hybridisation " << at_c->getHybridization() << std::endl;
      
      if (bond_2) {
	 if (bond_3) {
	    at_n1->setHybridization(RDKit::Atom::SP2);
	    at_c->setHybridization(RDKit::Atom::SP2);
	    bond_2->setBondType(RDKit::Bond::ONEANDAHALF);
	    bond_3->setBondType(RDKit::Bond::ONEANDAHALF);
	    std::cout << "deloced bonds 2 and 3" << std::endl;
	    //
	    if (bond_1) {
	       std::cout << "deloced bond 1 " << std::endl;
	       bond_1->setBondType(RDKit::Bond::ONEANDAHALF);
	    }

	 }
      }
   }

   
}
